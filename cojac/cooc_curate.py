#!/usr/bin/env python3
# import numpy as np
# import pandas as pd
import os
import sys
import re
import requests

# import csv
import json
import yaml

import click

# TODO proper library instead of copy-pasta
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
        # sort so that definitions with same mutations in different categories all yield identic result
        # (e.g.: YAML definitions with either "{mut:{123:A,567:C}}" or "{mut:{567:C},extra:{123:A}}"  or "{mut:{123:A},extra:{567:C}}"
        # will all always consistently yield a tmp_mut_dict with "{123:A,567:C}" )
        "mut": dict(
            sorted(
                (pos, mut_decode(mut))
                for c in categories
                if c in yam
                for (pos, mut) in yam[c].items()
                if mut_decode(mut)
            )
        ),
    }


# END of mutbamscan copy-pasta


#
# CoV-Spectrum pico-API
#

server = "https://cov-spectrum.ethz.ch/gisaid/api/v1"


def aggregated(**kwargs):
    """use the 'aggregated' LAPIS call to get
    - optionnally categorized ('fields' parameters)
    - count of samples that fit criteria (the other parameters)
    """
    # TODO proper error handling
    return json.loads(requests.get(f"{server}/sample/aggregated", params=kwargs).text)[
        "data"
    ]


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


def collapse(vardict, sublineages, combinedname=None):
    """collapse multiple sublineages into a combined one
    e.g. BA.1 + BA.1.1 + BA.2 + BA.3 ... => B.1.529/omicron
    """
    if combinedname is None:
        combinedname = sublineages[0]

    combined = {v: c for v, c in vardict.items() if v not in sublineages}
    combined[combinedname] = sum([vardict[v] for v in sublineages])

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
        muts = collapse(muts, sublineages, combined)
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
    help="helps determining specific mutations and cooccurences by querying CoV-Spectrum. This tool queries LAPIS, see https://lapis.cov-spectrum.org/swagger/ and https://lapis.cov-spectrum.org/"
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
            common_kwargs["alllin"] = collapse(alllin, sublineages, combinedname)
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
