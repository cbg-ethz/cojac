#!/usr/bin/env python3
import numpy as np
import pandas as pd
import sys
import os

# import re
import argparse
import csv
import json
import yaml
import gzip
import pysam
import re


def main():
    # parse command line
    argparser = argparse.ArgumentParser(
        fromfile_prefix_chars="@",  # support pass BAM list in file instead of parameters
        description="scan amplicon (covered by long read pairs) for mutation cooccurrence",
        epilog="@listfile can be used to pass a long list of parameters (e.g.: a large number of BAMs) in a file instead of command line",
    )
    inputgroup = argparser.add_mutually_exclusive_group(required=True)
    inputgroup.add_argument(
        "-s",
        "--samples",
        metavar="TSV",
        type=str,
        dest="samples",
        help="V-pipe samples list tsv",
    )
    inputgroup.add_argument(
        "-a",
        "--alignments",
        metavar="BAM/CRAM",
        nargs="+",
        dest="alignments",
        help="alignment files",
    )
    argparser.add_argument(
        "-/",
        "--batchname",
        nargs="?",
        const="/",
        default=None,
        dest="batchname",
        help="concatenate samplename/batchname from samples tsv",
    )
    argparser.add_argument(
        "-p",
        "--prefix",
        metavar="PATH",
        required=False,
        default="working/samples",
        type=str,
        dest="prefix",
        help="V-pipe work directory prefix for where to look at align files when using TSV samples list",
    )
    argparser.add_argument(
        "-r",
        "--reference",
        metavar="REFID",
        required=False,
        default=None,  # default='NC_045512.2',
        type=str,
        dest="rq_chr",
        help="reference to look for in alignment files",
    )
    argparser.add_argument(
        "-m",
        "--vocdir",
        metavar="DIR",
        required=False,
        default="./voc",
        type=str,
        dest="vocdir",
        help="directory containing the yamls defining the variant of concerns",
    )
    argparser.add_argument(
        "-b",
        "--bedfile",
        metavar="BED",
        required=False,
        default="./nCoV-2019.insert.bed",
        type=str,
        dest="bedfile",
        help="bedfile defining the amplicons, with format: ref\\tstart\\tstop\\tamp_num\\tpool\\tstrand",
    )
    argparser.add_argument(
        "-#",
        "--cooc",
        metavar="COOC",
        required=False,
        default=2,
        type=int,
        dest="cooc",
        help="minimum number of cooccurences to search for",
    )
    argparser.add_argument(
        "-j",
        "--json",
        metavar="JSON",
        required=False,
        default=None,
        type=str,
        dest="json",
        help="output results to as JSON file",
    )
    argparser.add_argument(
        "-y",
        "--yaml",
        metavar="YAML",
        required=False,
        default=None,
        type=str,
        dest="yaml",
        help="output results to as yaml file",
    )
    argparser.add_argument(
        "-t",
        "--tsv",
        metavar="TSV",
        required=False,
        default=None,
        type=str,
        dest="tsv",
        help="output results to as (raw) tsv file",
    )
    argparser.add_argument(
        "-d",
        "--dump",
        action="store_true",
        dest="dump",
        help="dump the python object to the terminal",
    )
    args = argparser.parse_args()

    if args.samples is not None:
        assert os.path.isfile(args.samples), f"cannot find sample file {args.samples}"
    else:
        for align in args.alignments:
            assert os.path.isfile(align), f"cannot find alginment file {align}"

    def test_read(read, mut_dict):
        """
        test if mutations listed in mut_dict are present in the pysam read

        returns a list with:
                found_site: list of site present in the read (no matter content)
                found_mut: list of those position which have the mutations variant
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
            "sites": dict(zip(sites_cnt[0].tolist(), sites_cnt[1].tolist()))
            if len(all_sites)
            else {},
            "muts": dict(zip(muts_cnt[0].tolist(), muts_cnt[1].tolist()))
            if len(all_muts)
            else {},
        }

    # TODO better naming
    def findbam(prefix, batch, sample):
        """function to find bam from prefix, batch and sample"""
        alnfname = os.path.join(prefix, sample, batch, "alignments/REF_aln.bam")
        assert os.path.isfile(alnfname), f"cannot find alginment file {alnfname}"
        return alnfname

    def scanbam(alnfname, amplicons, rq_chr):
        """scan a bamfile found at alnfname"""
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
                        for each amplicon where cooccurences could be found
        """
        mut_df = pd.DataFrame(
            {"position": mut_dict.keys(), "mutation": mut_dict.values()}
        )
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
                amplicons_dict["{}_{}".format(i + 1, voc_name)] = [
                    start,
                    end,
                    qstart,
                    qstop,
                    tmp_mut_dict,
                ]

        return amplicons_dict

    rxmutdec = re.compile(
        "^(?:(?:(?:(?P<ref>[ATCG]+)\>)?(?P<mut>[ATCG]+))|(?P<del>[\-]+)|(?:[\+](?P<ins>[ATGC]+)))$"
    )

    def mut_decode(mutstr):
        """function to decode the strings of a mutation
        mutation can specifie:
         - a single base (e.g.: 'C') or a run of bases ('CAT')
         - a substitution (e.g.: 'G>C'), can also be longer ('GTC>CAT')
         - a run of deletions ('---')
         - an insertions ('+TATA')

        (Basically, this function strips the reference part)

        Input:
                a string with any of the above
        Returns:
                a string with bases or deletions (either single or runs)
        """

        res = rxmutdec.match(mutstr)

        if res:
            match = res.groupdict()
            if match["mut"]:
                # print(f"{mutstr} : mutation {match['mut']}")
                return match["mut"]
            if match["del"]:
                # print(f"{mutstr} : deletion {match['del']}")
                return match["del"]
            if match["ins"]:
                print(f"insertions not supported (yet): {mutstr} : {match['ins']}")
                return None
        print(f"cannot parse mutation f{mutstr}")
        sys.exit(1)

    def filter_decode_vartiant(yam, categories=["mut", "extra", "shared", "subset"]):
        """function to filter mutations and decode the strings
        mutations can be classified as:
         - 'mut': defining mut
         - 'extra': officially not counted as defining but still specific
         - 'shared': present in this variant, but also occuring outside of it (e.g.: common ancestror)
         - 'subset': only present on some variants (i.e.: 'not fixed')
        mutation can be specifie:
         - a single base (e.g.: 'C') or a run of bases ('CAT')
         - a substitution (e.g.: 'G>C'), can also be longer ('GTC>CAT')
         - a run of deletions ('---')
         - an insertions ('+TATA')

        Input:
                yam: dictionnary loaded from a voc YAML
                categories: a list of categories to filter-in
        Returns:
                voc: dict() with
                        'name': short name of the variant
                        'mut': a an ordereddict() of 'position': 'bases or deletions'
        """

        return {
            "name": yam["variant"]["short"],
            "mut": {
                pos: mut_decode(mut)
                for c in categories
                if c in yam
                for (pos, mut) in yam[c].items()
                if mut_decode(mut)
            },
        }

    # load bedfile
    bedfile = args.bedfile
    amp_bed = pd.read_table(
        bedfile, names=["ref", "start", "stop", "amp_num", "pool", "strand"]
    )
    # make query start and query end for each amplicon
    amp_bed["qstart"] = (
        pd.concat([amp_bed["start"][0:1], amp_bed["stop"][:-1]]).reset_index(drop=True)
        + 30
    )
    amp_bed["qstop"] = (
        pd.concat([amp_bed["start"][1:], amp_bed["stop"][-1:]]).reset_index(drop=True)
        - 30
    )
    tmp_filter = amp_bed["qstart"] >= amp_bed["qstop"]
    amp_bed.loc[tmp_filter, "qstart"] = (
        (amp_bed["start"] + amp_bed["stop"]) / 2 - 5
    ).astype("int")
    amp_bed.loc[tmp_filter, "qstop"] = (
        (amp_bed["start"] + amp_bed["stop"]) / 2 + 5
    ).astype("int")

    # load all voc yamls
    vocdir = args.vocdir
    loaded_yamls = []
    for path in [p for p in os.listdir(vocdir) if not p.startswith(".")]:
        full_path = os.path.join(vocdir, path)
        with open(full_path, "r") as file:
            loaded_yaml = yaml.load(file, Loader=yaml.FullLoader)
        loaded_yamls.append(filter_decode_vartiant(loaded_yaml))

    # make amplicon dict for each voc
    amplicons = {}
    numbers = {}
    for yam in loaded_yamls:
        amp_dict = make_amplicons_dict(amp_bed, yam["mut"], yam["name"], args.cooc)

        # merge identical amplicons from different variants
        # e.g.: P.1 and 501Y.V2 share amplicon 76 with following mutations
        #  - G23012A  E484K
        #  - A23063T  N501Y
        # both on Spike's RBD
        for (k1, d) in list(amp_dict.items()):  # WARNING inefficient but functional
            if d in amplicons.values():
                k2 = list(amplicons.keys())[list(amplicons.values()).index(d)]
                print(f"{k1} is identical to {k2}")
                amplicons["_".join([k2, yam["name"]])] = amplicons.pop(k2)
                del amp_dict[k1]
                break

        # merge dicts
        if len(amp_dict):
            amplicons.update(amp_dict)

    rq_chr = args.rq_chr  # e.g.: 'NC_045512.2'

    # loop for if samples are given through the TSV list
    table = {}
    if args.samples is not None:
        i = 0
        with open(
            args.samples, "rt", encoding="utf-8", newline=""
        ) as tf:  # this file has the same content as the original experiment
            for r in csv.reader(tf, delimiter="\t"):  # dialect='excel-tab'):
                sample, batch = r[:2]
                print(sample)
                alnfname = findbam(args.prefix, batch, sample)
                table[
                    f"{sample}{args.batchname}{batch}" if args.batchname else sample
                ] = scanbam(alnfname, amplicons, rq_chr)

    # loop for if samples are given through -a option
    # this option can also de used to dispatch per sample jobx on the cluster
    else:
        for alnfname in args.alignments:
            sample = alnfname  # HACK use the whole BAM file as sample name
            table[sample] = scanbam(alnfname, amplicons, rq_chr)

    #
    # dumps, for being able to take it from here
    #
    if (args.dump) or (not (args.json or args.yaml)):
        print(table)
    if args.json:
        with open(args.json, "wt") as jf:
            json.dump(obj=table, fp=jf)
    if args.yaml:
        with open(args.yaml, "wt") as yf:
            print(yaml.dump(table, sort_keys=False), file=yf)

    # raw data into pands.DataFrame
    if args.tsv:
        raw_table_df = pd.DataFrame.from_dict(
            data={
                i: {(j, k): table[i][j][k] for j in table[i] for k in table[i][j]}
                for i in table
            },
            orient="index",
        )
        # with pd.option_context('display.max_rows', None): #, 'display.max_columns', None):
        # print(raw_table_df)
        raw_table_df.to_csv(args.tsv, sep="\t", compression={"method": "infer"})


if __name__ == "__main__":
    main()
