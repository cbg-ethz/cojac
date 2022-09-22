#!/usr/bin/env python3

import typing
import re

import click
import yaml
import csv

from .cooc_curate import listfilteredmutations


# regex
parsenuc = re.compile(
    "^(?P<orig>[ATGC])(?P<pos>\d+)(?P<mut>[\-ATCG])$", flags=re.IGNORECASE
)


@click.command(
    help="Helps generating a list of mutations frequently found in a variant by querying CoV-Spectrum",
    epilog="This tool queries LAPIS, see https://lapis.cov-spectrum.org/swagger/ and https://lapis.cov-spectrum.org/",
)
@click.option(
    "--var",
    "--variant",
    "var",
    metavar="PANGO",
    required=True,
    default=None,
    type=str,
    help="Pangolineage of the root variant to list",
)
@click.option(
    "--extras",
    metavar="LAPIS",
    default=None,
    type=str,
    help="Additional LAPIS query arguments passed as a YAML flow, e.g.: '"
    '{dateFrom: "2022-02-01", variantQuery: "[6-of: S:147E, S:152R, S:157L, S:210V, S:257S, S:339H, S:446S, S:460K, ORF1a:1221L, ORF1a:1640S, ORF1a:4060S]"}'
    "'. For more information about LAPIS, see: https://lapis-docs.readthedocs.io/en/latest/",
)
@click.option(
    "-f",
    "--minfreq",
    "minfreq",
    metavar="FREQ",
    default=0.8,
    type=float,
    help="Minimum frequency for inclusion in list",
)
@click.option(
    "-d",
    "--mindelfreq",
    "mindelfreq",
    metavar="FREQ",
    default=None,
    type=float,
    help="Use a different minimum frequency for deletions (useful early on when there are few sequences and some of those were produced by pipeline that don't handle deletions)",
)
@click.option(
    "-s",
    "--minseqs",
    "minseqs",
    metavar="NUM",
    default=100,
    type=int,
    help="Minimum number of sequence supporting for inclusion in list",
)
@click.option(
    "--covariants",
    "covariants",
    metavar="TSV",
    default=None,
    type=str,
    help="import from a covariants.org TSV file instead of covSpectrum. (See: https://github.com/hodcroftlab/covariants/blob/master/defining_mutations/)",
)
@click.option(
    "--debug/--no-debug",
    "debug",
    default=False,
)
def sig_generate(var, minfreq, mindelfreq, minseqs, extras, covariants, debug):
    if debug:
        import sys

        print(
            "extra:\n", yaml.load(extras, Loader=yaml.FullLoader), "\n", file=sys.stderr
        )

    # get initial list
    mut_request = []
    if not covariants:
        # fetch it from covSpectrum
        mut_request = listfilteredmutations(
            var,
            minfreq=minfreq,
            mindelfreq=mindelfreq,
            minseqs=minseqs,
            extras=yaml.load(extras, Loader=yaml.FullLoader) if extras else {},
        )
    else:
        # decode it from a TSV file from covariants.org
        with open(covariants, newline="") as csvfile:
            reader = csv.DictReader(csvfile, dialect="excel-tab")
            mut_request = [row["nuc_change"] for row in reader]
            # TODO handle all the additional metadata in other columns, like in tsv2cojac.py prototype

    mut_record = typing.NamedTuple("mut_record", [("mut", str), ("orig", str)])

    # load the mutation found by the request
    mutlist = {}
    for mut in mut_request:
        # parse single nucleotide mutation
        try:
            m = parsenuc.match(mut).groupdict()
        except:
            print(f"Cannot parse <{mut}>\n")
            next

        mutlist[int(m["pos"])] = mut_record(mut=m["mut"], orig=m["orig"])

    # sort and merge/extend consecutive
    mutmerged = {}
    for p, m in sorted(mutlist.items(), reverse=True):
        # consecutive positions ?
        nxt = p + 1
        if nxt in mutmerged:
            n = mutmerged[nxt]
            # are we extending a deletion? or nucleotides?
            if (m.mut == "-" and n.mut[0] == "-") or (m.mut != "-" and n.mut[0] != "-"):
                del mutmerged[nxt]
                mutmerged[p] = mut_record(mut=m.mut + n.mut, orig=m.orig + n.orig)
                continue

        # new mutation, not extending
        mutmerged[p] = m

    # TODO check neighbors for repeat patterns (i.e.: alternate deletion position span)

    # display
    for p, m in sorted(mutmerged.items()):
        if m.mut == "-" * len(m.mut):
            # deletion
            print(f"  {p}: '{m.mut}'")
        else:
            print(f"  {p}: '{m.orig}>{m.mut}'")
