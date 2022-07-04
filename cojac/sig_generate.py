#!/usr/bin/env python3

import click

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
#@click.argument("var", nargs=-1)

def sig_generate(var):
    for m in listfilteredmutations(var):
        print(m)

