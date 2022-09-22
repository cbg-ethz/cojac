#!/usr/bin/env python3
# import numpy as np
# import pandas as pd
import os
import sys
import re
import requests
from urllib.parse import urlparse

# import csv
import json
import yaml

import click

from .mut_parser import mut_decode, filter_decode_vartiant


#
# CoV-Spectrum pico-API
#


server = "https://lapis.cov-spectrum.org/gisaid/v1"


def getAccessKey():
    # TODO proper global handling

    # envrionment, e.g.: test workflow uses that
    if "COVSPECTRUM_ACCESSKEY" in os.environ:
        return os.environ["COVSPECTRUM_ACCESSKEY"]

    try:
        # we use ~/.netrc to obtain credentials
        import netrc

        return netrc.netrc().authenticators(urlparse(server).netloc)[2]
    except:
        print(
            f"no access key found for cov-spectrum.ethz.ch in ~/.netrc\n"
            f"your requests are likely to be not allowed by the server\n"
            f"please add the following entry in your ~/.netrc:\n"
            f"```\n"
            f"machine {urlparse(server).netloc}\n"
            f"password <YOUR_ACCESS_KEY_HERE>\n"
            f"```\n"
            f"or use the COVSPECTRUM_ACCESSKEY envrionment variable\n",
            file=sys.stderr,
        )
        # sys.exit(1)


def checkerror(reply):
    if reply.get("data", None) is None:
        # TODO prety-printing would help
        print("Error from server:", json.dumps(reply, indent=2))
        # TODO replace with proper exception throwing in the future
        sys.exit(1)

    return reply


def nucmutations(**kwargs):
    """use the 'nuc-mutation' LAPIS call to get
    - optionnally filtered (filter parameters such as 'pangoLineage')
    - list of mutation, frequencies and sample counts (for all frequencies > 0.05 )
    """
    kwargs["accessKey"] = getAccessKey()
    reply = json.loads(
        requests.get(f"{server}/sample/nuc-mutations", params=kwargs).text
    )

    checkerror(reply)

    return reply["data"]


def listmutations(lineage, extras={}):
    """list all mutation for a specific lineage (either the variant it self or one of its sub-variants)"""
    if "variantQuery" in extras.keys():
        # NOTE according to covspectrum: "Please specify the variant either by using the fields pangoLineage, nextstrainClade, gisaidClade, aaMutations and nucMutations, or by using variantQuery - don't use both at the same time."
        return nucmutations(**extras)
    return nucmutations(pangoLineage=(lineage + "*"), **extras)


def aggregated(**kwargs):
    """use the 'aggregated' LAPIS call to get
    - optionnally categorized ('fields' parameters)
    - count of samples that fit criteria (the other parameters)
    """
    kwargs["accessKey"] = getAccessKey()
    reply = json.loads(requests.get(f"{server}/sample/aggregated", params=kwargs).text)

    checkerror(reply)

    return reply["data"]


def listsublineages(lineage):
    """list all lineages which are known to be either the variant it self or one of its sub-variants"""
    return set(
        s["pangoLineage"]
        for s in aggregated(pangoLineage=(lineage + "*"), fields="pangoLineage")
    )


def listalllineages():
    """structure listing known lineages and their respective total count of samples"""
    return {s["pangoLineage"]: s["count"] for s in aggregated(fields="pangoLineage")}


def mutsinlineages(*mutations):
    """looks for presence of specific mutations,
    return lineages and count of samples carrying them"""
    return {
        s["pangoLineage"]: s["count"]
        for s in aggregated(nucMutations=",".join(mutations), fields="pangoLineage")
    }


#
# Compute useful stuff
#


def listfilteredmutations(
    lineage, minfreq=0.8, minseqs=100, mindelfreq=None, extras={}
):
    """select mutations to make a signature
    - optionnally uses an alternate minimal frequency for deletions
      (useful early on when there are few sequences, some of them produced
      by pipelines which don't handle deletions.)
    """
    return set(
        m["mutation"]
        for m in listmutations(lineage, extras)
        if m["proportion"]
        > (minfreq if mindelfreq is None or m["mutation"][-1] != "-" else mindelfreq)
        and m["count"] > minseqs
    )


def collapse_sublineages(vardict, sublineages, combinedname=None):
    """collapse multiple sublineages into a combined one
    e.g. BA.1 + BA.1.1 + BA.2 + BA.3 ... => B.1.529/omicron
    """
    if combinedname is None:
        combinedname = sublineages[0]

    combined = {v: c for v, c in vardict.items() if v not in sublineages}
    combined[combinedname] = sum([vardict.get(v, 0) for v in sublineages])

    return combined


def freqperlineage(mutations, lineages):
    """convert count of samples to fraction of samples by dividing mutations-carrying by total count
    results are sorted by frequency"""
    return dict(
        sorted(
            [(v, (c / lineages[v])) for v, c in mutations.items()],
            key=lambda var: var[1],
            reverse=True,
        )
    )


def rank(sublineages, freqs, limit=0.05):
    """find ranks of occurence of:
    - last sublineages still part of the seeked lineages
    - first other lineages which is not part of the seeked one
    - last member still above the limit
    Note: ranks are 1-based
    """
    lastsub = None
    firstother = None
    inrange = None

    s = set(sublineages)
    i = 0
    for v, f in freqs.items():
        i += 1
        if v in s:
            lastsub = i
            s.remove(v)
        else:
            if firstother is None:
                firstother = i

        if f >= limit:
            inrange = i

        if len(s) == 0 and firstother is not None and f < limit:
            break

    return (lastsub, firstother, inrange)


def curate_muts(
    alllin, sublineages, label, mutargs, low=0.05, high=0.80, ansi=False, combined=None
):
    """
    make a currated list:
    - label (printed on screen)
    - mutargs: list of mutations to search for
    - sublineages: lineage that they must be found in
    - alllin: table with all lineage and their counts (used to compute fractions)
    - combined:
      - None: use mutation prevalage list as-is
      - string: combine all counts of sublineages into a single entry called this way
    """
    marker = "\x1b[1;3;4m*" if ansi else "*"
    endl = "\x1b[0m" if ansi else ""

    muts = mutsinlineages(*mutargs)
    if combined is not None:
        # group the counts of sublineages together into the main lineage
        # e.g. combine counts of BA.1,BA.1.1,BA.2,BA.3,etc. => B.1.1.529*
        muts = collapse_sublineages(muts, sublineages, combined)
        su = {combined}
    else:
        su = set(sublineages)
    # print(muts)
    freqs = freqperlineage(muts, alllin)
    # print(freqs)
    (s, o, r) = rank(su, freqs, limit=low)
    # print(s,o,r)

    fi = [i for i in freqs.items()]
    # print({p:fi[p-1] for p in [s,o,r]})

    isspecific = (o is None) or (
        (s < o) and (fi[s - 1][1] >= high) and (fi[o - 1][1] <= low)
    )
    print(f"{marker if isspecific else ''}{label}: ", end=endl)
    end = min(
        15,
        len(fi),
        max(s + 1, o + 1 if o is not None else 0),
        (r + 1 if r is not None else 9999),
    )
    for i in range(0, end):
        (v, f) = fi[i]
        print(f"{', ' if i else ''}{marker if v in su else ''}{v}={f:.2f}", end=endl)

    # way after max
    if s > end:
        v, f = fi[s - 1]
        print(f" {chr(0x2026)} [{s}]:{marker}{v}={f:.2f}", end=endl)

    print()


@click.command(
    help="Helps determining specific mutations and cooccurences by querying CoV-Spectrum",
    epilog="This tool queries LAPIS, see https://lapis.cov-spectrum.org/swagger/ and https://lapis.cov-spectrum.org/",
)
@click.option(
    "-a",
    "--amplicons",
    "amp",
    metavar="YAML",
    required=False,
    default=None,
    type=str,
    help="use the YAML file generated by mutbamscan to query amplicons instead of mutations",
)
@click.option(
    "-m",
    "--mutations",
    "domuts",
    required=False,
    is_flag=True,
    help="always do mutations (even if amplicons YAML provided)",
)
@click.option(
    "-H",
    "--high",
    required=False,
    default=0.80,
    type=float,
    help="Fraction above which a mutation must be found among seeked lineages",
)
@click.option(
    "-l",
    "--low",
    required=False,
    default=0.20,
    type=float,
    help="Fraction under which a mutation must be found among other lineages",
)
@click.option(
    "--collapse/--no-collapse",
    required=False,
    default=True,
    help="combine counts of all sublineages together and consider a signle value that corresponds to a lineages family (e.g.: count all B.1.612.2* together). This is especially useful for assessing signature of old variants that have branched out by now.",
)
@click.option(
    "--colour/--no-colour",
    required=False,
    default=True,
    help="use coloured output",
)
@click.argument("voc", nargs=-1)
def cooc_curate(amp, domuts, high, low, collapse, colour, voc):
    amplicons = None
    if amp:
        with open(amp, "r") as yf:
            amplicons = yaml.load(yf, Loader=yaml.FullLoader)

    alllin = listalllineages()
    const_kwargs = {
        "alllin": alllin,
        "low": low,
        "high": high,
        "ansi": colour,
    }

    for voc_fname in voc:
        with open(voc_fname, "r") as yf:
            voc = yaml.load(yf, Loader=yaml.FullLoader)
        var = filter_decode_vartiant(
            voc, categories=["mut", "extra"]
        )  # shared and subset are known to not map perfectly anyway

        lineage = voc["variant"]["pangolin"]  # 'BA.1'

        # find all pangolin lineages which are known to be part of this lineage
        # e.g.: B.1.1.529 (omikron) is at the time of this writing BA.1,BA.1.1,BA.2 and BA.3
        sublineages = listsublineages(lineage)
        print(f"{lineage}: {','.join(sublineages)}")

        common_kwargs = {
            "sublineages": sublineages,
        }
        common_kwargs.update(const_kwargs)
        if len(sublineages) > 1 and collapse:
            # do not use the common list of all lineage, instead collapse those who are part of the current considered variant into a single category
            # e.g. combine counts of BA.1,BA.1.1,BA.2,BA.3,etc. => B.1.1.529*
            combinedname = f"{lineage}*"
            common_kwargs["alllin"] = collapse_sublineages(
                alllin, sublineages, combinedname
            )
            common_kwargs["combined"] = f"{lineage}*"
        else:
            common_kwargs["combined"] = None

        # mutations
        if domuts or amplicons is None:
            print("\nmutations:")
            for p, v in var["mut"].items():
                # on-screen: compact, e.g. 28881AAC
                label = f"{p}{v}"
                # requests: individual, e.g. 28881A,28882A,28883C
                mutargs = [f"{p2}{v2}" for p2, v2 in zip(range(p, p + len(v)), v)]

                curate_muts(label=label, mutargs=mutargs, **common_kwargs)

        # amplicons
        if amplicons is not None:
            print("\namplicons:")
            for aname, a in amplicons.items():
                # only amplicons at least listed for that variant
                if var["name"] not in aname.split("_")[1:]:
                    continue

                # on-screen: name + all compact mutation
                label = f"{aname}[{ ','.join(list(str(p)+v for p,v in a[4].items())) }]"
                # build full list of individual mutations
                mutargs = [
                    f"{p2}{v2}"
                    for p, v in a[4].items()
                    for p2, v2 in zip(range(p, p + len(v)), v)
                ]

                curate_muts(label=label, mutargs=mutargs, **common_kwargs)

        if len(voc) > 1:
            print()


if __name__ == "__main__":
    cooc_curate()
