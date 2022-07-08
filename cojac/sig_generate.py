#!/usr/bin/env python3

import click
import yaml

from .cooc_curate import listfilteredmutations


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
    help="Additional LAPIS query arguments passed as a YAML flow, e.g.: 6-of: S:147E, S:152R, S:157L, S:210V, S:257S, S:339H, S:446S, S:460K, ORF1a:1221L, ORF1a:1640S, ORF1a:4060S }",
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
    "-f",
    "--minseqs",
    "minseqs",
    metavar="NUM",
    default=100,
    type=int,
    help="Minimum number of sequence supporting for inclusion in list",
)
@click.option(
    "--debug/--no-debug",
    "debug",
    default=False,
)
# @click.argument("var", nargs=-1)


def sig_generate(var, minfreq, minseqs, extras, debug):
    if debug:
        import sys

        print(
            "extra:\n", yaml.load(extras, Loader=yaml.FullLoader), "\n", file=sys.stderr
        )
    for m in listfilteredmutations(
        var,
        minfreq=minfreq,
        minseqs=minseqs,
        extras=yaml.load(extras, Loader=yaml.FullLoader) if extras else {},
    ):
        print(m)
