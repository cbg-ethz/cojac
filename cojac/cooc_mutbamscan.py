#!/usr/bin/env python3
import numpy as np
import pandas as pd
import sys
import os
import re
import csv
import json
import yaml
import gzip

# import pysam # HACK pysam isn't available on bioconda aarch64, yet. But loading it here cause every other function of cojac to fail, too.

import click

from .mut_parser import mut_decode, filter_decode_vartiant


def test_read(read, mut_dict):
    """
    test if mutations listed in mut_dict are present in the pysam read

    returns a list with:
            found_site:	list of site present in the read (no matter content)
            found_mut:	list of those position which have the mutations variant
    """

    # WARNING pysam is 0-based! (see here: https://pysam.readthedocs.io/en/latest/faq.html#pysam-coordinates-are-wrong )

    # 1. check which mutation' sites are in range of that read
    found_site = [
        p
        for p, m in mut_dict.items()
        if read.reference_start <= (p - 1) <= (read.reference_end - len(m))
    ]
    if not len(found_site):
        return (None, None)  # sites aren't present no point checking variants

    # 2. of those sites, check which content mutations' variants
    read_dict = dict(
        zip(read.get_reference_positions(full_length=True), read.query_sequence)
    )

    found_mut = []
    for p in found_site:
        l = len(mut_dict[p])  # lenght
        # base present ?
        srch = [(p - 1 + i) in read_dict for i in range(l)]
        if all(srch):  # all positions found!
            # check if it's the expected mutation(s)
            if all([read_dict[p - 1 + i] == mut_dict[p][i] for i in range(l)]):
                found_mut.append(p)
        elif not any(srch):  # none position found! (entire deletion)
            # check if we're hunting for a string of deletions (-)
            if "-" * l == mut_dict[p]:
                found_mut.append(p)
        # TODO give some thoughs about partial deletions

    if len(found_mut):  # found mutation as sites
        return (found_site, found_mut)
    else:  # sites present, but no mutation found
        return (found_site, None)


# scan an amplicon for a specific set of mutations
def scanamplicon(read_iter, mut_dict):
    # TODO inefficient, could be streamed (with an accumulator for pairs)
    reads = dict()
    for read in read_iter:
        name = str(read.query_name)
        R = "R1" if read.is_read1 else "R2"
        if name in reads:
            reads[name][R] = read
        else:
            reads[name] = {R: read}

    print("amplion:", len(reads))

    # tally the mutation sites and the presence of variant in there accross all reads
    # TODO merge with above
    all_muts = []
    all_sites = []
    for rds in reads:
        val = reads[rds]
        site_out = []
        mut_out = []
        for s in val:
            R = val[s]
            (t_pos, t_read) = test_read(R, mut_dict)
            if t_pos is not None:
                site_out.extend(t_pos)
            if t_read is not None:
                mut_out.extend(t_read)
        if len(site_out):
            all_sites.append(site_out)
            if len(mut_out):
                all_muts.append(mut_out)

    sites_cnt = np.unique([len(set(mut)) for mut in all_sites], return_counts=True)
    print("sites:", len(all_sites), sites_cnt)
    muts_cnt = np.unique([len(set(sit)) for sit in all_muts], return_counts=True)
    print("muts:", len(all_muts), muts_cnt)

    # look at last column only
    return {
        "sites": (
            dict(zip(sites_cnt[0].tolist(), sites_cnt[1].tolist()))
            if len(all_sites)
            else {}
        ),
        "muts": (
            dict(zip(muts_cnt[0].tolist(), muts_cnt[1].tolist()))
            if len(all_muts)
            else {}
        ),
    }


# TODO better naming
def findbam(prefix, batch, sample):
    """function to find bam from prefix, batch and sample"""
    alnfname = os.path.join(prefix, sample, batch, "alignments/REF_aln.bam")
    assert os.path.isfile(alnfname), f"cannot find alginment file {alnfname}"
    return alnfname


def scanbam(alnfname, amplicons, rq_chr):
    """scan a bamfile found at alnfname"""
    import pysam  # HACK pysam isn't available on bioconda aarch64, yet. So hot-load it only in the function that requires it. This lets all other parts of cojac working without it.

    amp_results = {}
    with pysam.AlignmentFile(alnfname, "rb") as alnfile:
        if rq_chr == None:
            # autoguess reference from alignment
            # HACK only using the first one
            # TODO handle multiple fragments (in the request itself)
            rq_chr = alnfile.references[0]
            print(f"autodecting reference as {rq_chr}")
        for amp_name, amp in amplicons.items():
            (start, stop, rq_b, rq_e, mut_dict) = amp

            # we need at least 2 to compute co-occurrence
            if len(mut_dict) < 1:  # HACK 2:
                continue

            print(f"amplicon_{amp_name}", rq_b, rq_e, mut_dict, sep="\t", end="\t")

            amplicon_iter = alnfile.fetch(rq_chr, rq_b, rq_e)
            amp_results[amp_name] = scanamplicon(amplicon_iter, mut_dict)
    return amp_results


def make_amplicons_dict(amp_bed, mut_dict, voc_name="", cooc=2):
    """function to make a dictionnary of places to look for coocurences of mutations.
    Input:
            amp_bed: pd.DataFrame with columns "qstart" and "qstop" of the query start and stop positions of the amplicons.
            mut_df: dict with keys "position" and value "mutation" for each mutation to look for
            voc_name: name given to the variant of concern
    optional:
            cooc: amplicons must have at lest that many variants
    Returns:
            amplicons_dict: dict() with entries of the form:
                    *amplicon_number*_*voc_name* : [start, stop, {position1 : mutation1, position2 : mutation2 , ...}]
                    for each amplicon where cooccurrences could be found
    """
    mut_df = pd.DataFrame({"position": mut_dict.keys(), "mutation": mut_dict.values()})
    amplicons_dict = {}

    mincoocthresh = (
        cooc if len(mut_dict) > 1 else 1
    )  # HACK exception for early B.1 variant defined by a single mutation (A23403G : D614G aka Doug)
    for i in range(amp_bed.shape[0]):
        start, end, qstart, qstop = amp_bed.iloc[i].loc[
            ["start", "stop", "qstart", "qstop"]
        ]
        t_df = mut_df[(start <= mut_df["position"]) & (mut_df["position"] <= end)]
        tmp_mut_dict = dict()
        for j in range(t_df.shape[0]):
            tmp_mut_dict[t_df.iloc[j]["position"]] = t_df.iloc[j]["mutation"]

        if len(tmp_mut_dict) >= mincoocthresh:
            # sort so that definitions with same mutations in different categories all yield identic result
            # (e.g.: YAML definitions with either "{mut:{123:A,567:C}}" or "{mut:{567:C},extra:{123:A}}"  or "{mut:{123:A},extra:{567:C}}"
            # will all always consistently yield a tmp_mut_dict with "{123:A,567:C}" )
            amplicons_dict["{}_{}".format(i + 1, voc_name)] = [
                start,
                end,
                qstart,
                qstop,
                dict(sorted(tmp_mut_dict.items())),
            ]

    return amplicons_dict


def bed_load(bedfile, sort=True):
    """function to load a bedfile and compute the pysam request parameters

    In samtools' libhts and in pysam if you request between the start and end of a window,
    the library will return all the read-pairs in that partially overlap that window.
    In our case, if we submit a request spanning the whole amplicon's window, we will get
    all the read-pairs of that amplicon, BUT because the multiplex PCR are tiled, we would
    get also the read-pairs from neighbouring amplicons.

    e.g.: (where b: some base and p: base of a primer)

    bbbbbbbbbbppp----------------ppppbbbbbbbbb : amp 119 and 121
    ----pppbbbbbbbbbbbbbbbbbbbbbbbbbbbbppp---- : amp 120

           |--------------------------| : request

    With such a request we would get all the (whole) reads-pairs of all 3 amplicons.
    We would need to then check for each read received if it actually
    covers the mutation position.

    Whereas using the neighbouring amplicon start/stop and skipping -30/+30 to skip
    the primers:

    bbbbbbbbbbppp----------------ppppbbbbbbbbb : amp 119 and 121
    ----pppbbbbbbbbbbbbbbbbbbbbbbbbbbbbppp---- : amp 120
             | +30 >         < 30- | : skip
                   |---------| : request

    The only read-pairs that overlap with the request are the ones from amplicon 120
    We will get only (whole) read-pairs of that amplicon, so we don't need to
    filter-out reads that absolutely don't cover the positions.

    But in some protocols, e.g. Nimagen, the amplicons are so short that our usual
    -30/+30 overshoots and crosses.
    In that case we're doing +5/-5 from the centre.

    bbbbbbbbbbppp----------------ppppbbbbbbbbb : amp 119 and 121
    ----pppbbbbbbbbbbbbbbbbbbbbbbbbbbbbppp---- : amp 120
                   < 5- | +5 > : skip
                   |---------| : request
    """

    amp_bed_raw = pd.read_table(
        bedfile,
        index_col=False,
        names=["ref", "start", "stop", "amp_num", "pool", "strand"],
        dtype={
            "ref": "str",
            "start": "uint32",
            "stop": "uint32",
            "amp_num": "str",
        },
    )
    amp_bed = amp_bed_raw.sort_values(by=["ref", "start"]) if sort else amp_bed_raw

    # make query start and query end for each amplicon
    amp_bed["qstart"] = (
        pd.concat([amp_bed["start"][0:1], amp_bed["stop"][:-1]]).reset_index(drop=True)
        + 30
    )
    amp_bed["qstop"] = (
        pd.concat([amp_bed["start"][1:], amp_bed["stop"][-1:]]).reset_index(drop=True)
        - 30
    )
    # if they cross use centre +/- 5
    tmp_filter = amp_bed["qstart"] >= amp_bed["qstop"]
    amp_bed.loc[tmp_filter, "qstart"] = (
        (amp_bed["start"] + amp_bed["stop"]) / 2 - 5
    ).astype("int")
    amp_bed.loc[tmp_filter, "qstop"] = (
        (amp_bed["start"] + amp_bed["stop"]) / 2 + 5
    ).astype("int")

    return amp_bed


def make_all_amplicons(amp_bed, vocs, revert=False, n_cooc=2, subset_fix=True):
    """given a .BED file and a directory with YAMLs decribing VOCs,
    generates all amplicons to search that have at least n_cooc mutations
    """

    # load all voc yamls
    loaded_yamls = []
    for full_path in vocs:
        with open(full_path, "r") as yf:
            loaded_yaml = yaml.load(yf, Loader=yaml.FullLoader)
        loaded_yamls.append(
            filter_decode_vartiant(
                loaded_yaml,
                categories=(
                    ["revert", "mut", "extra", "shared", "subset"]
                    if revert
                    else ["mut", "extra", "shared", "subset"]
                ),
            )
        )

    # make amplicon dict for each voc
    amplicons = {}
    for yam in loaded_yamls:
        amp_dict = make_amplicons_dict(amp_bed, yam["mut"], yam["name"], n_cooc)

        # merge identical amplicons from different variants
        # e.g.: P.1 and 501Y.V2 share amplicon 76 with following mutations
        #  - G23012A  E484K
        #  - A23063T  N501Y
        # both on Spike's RBD
        for k1, d1 in list(amp_dict.items()):  # WARNING inefficient but functional
            ampname = k1.split("_")[0]
            for k2, d2 in list(amplicons.items()):
                npart = k2.split("_")
                # search same amplicon number
                if npart[0] != ampname:
                    continue

                # compare amplicon definitions
                if d1 != d2:
                    # if we're only interested in matches
                    continue

                print(
                    f"{k1} is identical to {k2}: {  ','.join([ f'{p}{b}' for p, b in d1[4].items() ])}"
                )
                # sort the variant names part of the list, to be consistent between calls, no matter the filesystem's on-disk order
                # (i.e.: no 76_AY42_IN2 vs 76_IN2_AY42)
                newk2 = "_".join([npart[0]] + sorted(npart[1:] + [yam["name"]]))
                print(f"\t{k2} => {newk2}")
                amplicons[newk2] = amplicons.pop(k2)
                del amp_dict[k1]
                break

        # merge dicts
        if len(amp_dict):
            amplicons.update(amp_dict)

    # now check for subsets situations
    if subset_fix:
        # build look-up map
        # mutations sets back to amplicons
        rlum = {}
        for k, d in amplicons.items():
            # nothing can be a subset of a set of one
            if len(d[4]) == 1:
                continue

            mutset = {f"{p}{b}" for p, b in d[4].items()}
            varset = set(k.split("_")[1:])
            ampname = k.split("_")[0]

            if ampname not in rlum:
                rlum[ampname] = []

            rlum[ampname] += [(k, mutset, varset)]

        #  now examine each amplicon in succession in search of subsets
        for k1, d in list(amplicons.items()):
            origk1 = k1
            ampmutset = {f"{p}{b}" for p, b in d[4].items()}
            ampname = k1.split("_")[0]

            # are there *any* longer than 1 set of which we could be subsets?
            if ampname not in rlum:
                continue

            ampmutnum = len(ampmutset)
            # compare against the whole set (inefficient, but well)
            for k2, mutset, varset in rlum[ampname]:
                # skip self
                if k2 == origk1:
                    # sanity check
                    assert (
                        mutset == ampmutset
                    ), f"something went wrong while building index: wrong self '{mutset}' vs '{ampmutset}'"
                    continue

                # it must be a larger set if we have any chance of being a *sub*set
                if len(mutset) <= ampmutnum:
                    continue

                if 0 == len(ampmutset - mutset):
                    print(f"{origk1}'s '{ampmutset}' is a subset of {k2}'s '{mutset}'")

                    ampvset = set(k1.split("_")[1:])

                    newk1 = "_".join([ampname] + sorted(list(ampvset.union(varset))))
                    # only update the label
                    k1 = newk1
                    continue

            # rename the entry using the final label
            if k1 != origk1:
                print(f"\t{origk1} => {k1}")
                amplicons[k1] = amplicons.pop(origk1)

    return amplicons


def load_all_amplicons(inamp):
    with open(inamp, "rt") as yf:
        # type: force convert into numpy
        return {
            a: [np.uint32(p) for p in q[:4]]
            + [{np.uint32(p): m for p, m in q[4].items()}]
            for a, q in yaml.safe_load(yf).items()
        }


def add_name_comments(in_str, amp_bed):
    out_str = ""
    rxnum = re.compile(r"^(?:\? +)?(?P<num>\d+)(?=_)")

    n = None
    for l in in_str.splitlines():
        find_n = rxnum.match(l)
        if find_n:
            n = int(find_n.groupdict().get("num", 0))
        if l[-1] == "]" and n:
            l += f" #\t{amp_bed['amp_num'][n - 1]}"
            n = None
        out_str += f"{l}\n"

    return out_str


def write_all_amplicons(amplicons, outamp, amp_bed=None):
    with open(outamp, "wt") as yf:
        # wrapper that force either flow style or block style
        class blockmap(dict):
            pass

        def blockmap_rep(dumper, data):
            return dumper.represent_mapping(
                "tag:yaml.org,2002:map", data, flow_style=False
            )

        class flowseq(list):
            pass

        def flowseq_rep(dumper, data):
            return dumper.represent_sequence(
                "tag:yaml.org,2002:seq", data, flow_style=True
            )

        yaml.add_representer(blockmap, blockmap_rep)
        yaml.add_representer(flowseq, flowseq_rep)

        # type: force convert numpy numeric into standard YAML and python integers
        yaml_str = yaml.dump(
            blockmap(
                {
                    a: flowseq(
                        [int(p) for p in q[:4]] + [{int(p): m for p, m in q[4].items()}]
                    )
                    for a, q in amplicons.items()
                }
            ),
            sort_keys=False,
        )

        if amp_bed is not None:
            yaml_str = add_name_comments(yaml_str, amp_bed)

        print(
            yaml_str,
            file=yf,
        )


@click.command(
    help="Scan amplicon (covered by long read pairs) for mutation cooccurrence",
    epilog="@listfile can be used to pass a long list of parameters (e.g.: a large number of BAMs) in a file instead of command line",
)
@click.option(
    "-a",
    "--alignments",
    metavar="BAM/CRAM",
    multiple=True,
    help="alignment files",
)
@click.option(
    "-n",
    "--name",
    metavar="NAME",
    multiple=True,
    default=None,
    help="when using alignment files, name to use for the output",
)
@click.option(
    "-s",
    "--samples",
    metavar="TSV",
    type=str,
    help="V-pipe samples list tsv",
)
@click.option(
    "--batchname",
    metavar="SEP",
    required=False,
    default=None,
    type=str,
    help="concatenate samplename/batchname from samples tsv",
)
@click.option(
    "-p",
    "--prefix",
    metavar="PATH",
    required=False,
    default="working/samples",
    type=str,
    help="V-pipe work directory prefix for where to look at align files when using TSV samples list",
)
@click.option(
    "-r",
    "--reference",
    "rq_chr",
    metavar="REFID",
    required=False,
    default=None,  # default='NC_045512.2',
    type=str,
    help="reference to look for in alignment files",
)
@click.option(
    "-m",
    "--vocdir",
    metavar="DIR",
    required=False,
    default=None,
    type=str,
    help="directory containing the yamls defining the variant of concerns",
)
@click.option(
    "-V",
    "--voc",
    metavar="VOC",
    multiple=True,
    default=None,
    type=str,
    help="individual yamls defining the variant of concerns",
)
@click.option(
    "--rev/--no-rev",
    "--with-revert/--without-revert",
    "revert",
    default=False,
    help="also include reverts when compiling amplicons (requires VOC YAML files with revert category)",
)
@click.option(
    "-b",
    "--bedfile",
    metavar="BED",
    required=False,
    default="./nCoV-2019.insert.V3.bed",
    type=str,
    help="bedfile defining the amplicons, with format: ref\\tstart\\tstop\\tamp_num\\tpool\\tstrand",
)
@click.option(
    "--sort/--no-sort",
    "sort",
    default=True,
    help="sort the bedfile by 'reference name' and 'start position' (default: sorted)",
)
@click.option(
    "-#",
    "--cooc",
    metavar="COOC",
    required=False,
    default=2,
    type=int,
    help="minimum number of cooccurrences to search for",
)
@click.option(
    "--fix-subset/--no-fix-subset",
    "--fs/--no-fs",
    "subset_fix",
    default=False,
    help="Fix variants attribution when cooccurrence are subset/superset of other variants",
)
# TODO: use mutually exclusive groups
@click.option(
    "-Q",
    "--amplicons",
    "--in-amp",
    "--in-amplicons",
    "inamp",
    metavar="YAML",
    required=False,
    default=None,
    type=str,
    help="use the supplied YAML file to query amplicons instead of building it from BED + voc's DIR",
)
@click.option(
    "-A",
    "--out-amp",
    "--out-amplicons",
    "outamp",
    metavar="YAML",
    required=False,
    default=None,
    type=str,
    help="output amplicon query in a YAML file",
)
@click.option(
    "--comment/--no-comment",
    "comment",
    default=True,
    help="add comments in the out amplicon YAML with names from BED file (default: comment the YAML)",
)
@click.option(
    "-j",
    "--json",
    "json_fname",
    metavar="JSON",
    required=False,
    default=None,
    type=str,
    help="output results to a JSON file",
)
@click.option(
    "-y",
    "--yaml",
    "yaml_fname",
    metavar="YAML",
    required=False,
    default=None,
    type=str,
    help="output results to a yaml file",
)
@click.option(
    "-t",
    "--tsv",
    metavar="TSV",
    required=False,
    default=None,
    type=str,
    help="output results to a (raw) tsv file",
)
@click.option(
    "-d",
    "--dump",
    is_flag=True,
    help="dump the python object to the terminal",
)
def cooc_mutbamscan(
    samples,
    alignments,
    name,
    batchname,
    prefix,
    rq_chr,
    vocdir,
    voc,
    revert,
    bedfile,
    sort,
    cooc,
    subset_fix,
    inamp,
    outamp,
    comment,
    json_fname,
    yaml_fname,
    tsv,
    dump,
):
    # amplicons that will be searched
    amplicons = {}
    if inamp is not None:
        # load pre-computed amplicons
        amplicons = load_all_amplicons(inamp)
    else:
        if vocdir:
            if not voc:
                voc = []
            voc += [
                os.path.join(vocdir, path)
                for path in [p for p in os.listdir(vocdir) if not p.startswith(".")]
            ]
        # compute amplicons
        assert (
            len(voc) > 0
        ), "Neither --voc nor --vocdir provided. Please provide variants of concern definition yamls."
        amp_bed = bed_load(bedfile, sort)
        amplicons = make_all_amplicons(
            amp_bed, voc, revert=revert, n_cooc=cooc, subset_fix=subset_fix
        )
        # and save them for future reference
        if outamp:
            write_all_amplicons(amplicons, outamp, amp_bed if comment else None)

    rq_chr = rq_chr  # e.g.: 'NC_045512.2'

    # loop for if samples are given through the TSV list
    table = {}
    if samples is not None:
        i = 0
        with open(
            samples, "rt", encoding="utf-8", newline=""
        ) as tf:  # this file has the same content as the original experiment
            for r in csv.reader(tf, delimiter="\t"):  # dialect='excel-tab'):
                sample, batch = r[:2]
                print(sample)
                alnfname = findbam(prefix, batch, sample)
                table[f"{sample}{batchname}{batch}" if batchname else sample] = scanbam(
                    alnfname, amplicons, rq_chr
                )

    # loop for if samples are given through -a option
    # this option can also de used to dispatch per sample jobx on the cluster
    elif alignments is not None:
        if name:
            assert len(alignments) == len(
                name
            ), f"Error: the number of BAMs/CRAMs files given to the -a/--alignments parameter and the number of NAMEs given to -n/--name must mach.\n{len(alignments)} BAM(s)/CRAM(s) given vs {len(name)} NAMEs"
        for alnfname, sample in zip(alignments, name if name else alignments):
            # HACK use the whole BAM file as sample name if no names provided
            table[sample] = scanbam(alnfname, amplicons, rq_chr)
    else:
        # we only wrote out the outamp and have nothing else to do.
        return

    #
    # dumps, for being able to take it from here
    #
    if (dump) or (not (json_fname or yaml_fname)):
        print(table)
    if json_fname:
        with open(json_fname, "wt") as jf:
            json.dump(obj=table, fp=jf)
    if yaml_fname:
        with open(yaml_fname, "wt") as yf:
            print(yaml.dump(table, sort_keys=False), file=yf)

    # raw data into pands.DataFrame
    if tsv:
        raw_table_df = pd.DataFrame.from_dict(
            data={
                i: {(j, k): table[i][j][k] for j in table[i] for k in table[i][j]}
                for i in table
            },
            orient="index",
        )
        # with pd.option_context('display.max_rows', None): #, 'display.max_columns', None):
        # print(raw_table_df)
        raw_table_df.to_csv(tsv, sep="\t", compression={"method": "infer"})


if __name__ == "__main__":
    cooc_mutbamscan()
