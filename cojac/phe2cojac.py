#!/usr/bin/env python3

import sys
import io
import re
import yaml
import strictyaml
import argparse

import click


@click.command()
@click.option(
    "-s",
    "--shortname",
    metavar="SHRT",
    required=False,
    default=None,
    type=str,
    help="shortname to use (otherwise auto-build one based on phe-genomic's unique id)",
)
@click.option(
    "-y",
    "--yaml",
    "outname",
    metavar="OUT_YAML",
    required=False,
    default=None,
    type=str,
    help="write cojac variant to a YAML file instead of printing (if empty, build filename from shortname)",
)
@click.argument(
    "fname",
    metavar="IN_YAML",
    type=str,
)
def phe2cojac(shortname, outname, fname):
    # phe-genomics input
    with open(fname, "r") as yf:
        yam = strictyaml.dirty_load(yf.read(), allow_flow_style=True).data

    # cojac yaml output

    # general metadata
    outy = {
        "variant": {"voc": yam["phe-label"], "pheuid": yam["unique-id"]},
        "source": yam["information-sources"],
        "threshold": int(yam["calling-definition"]["probable"]["mutations-required"]),
    }

    rxshortify = re.compile("^([a-z]{2})[^-]*-([a-z]{2})")

    outy["variant"]["short"] = (
        shortname
        if shortname
        else "".join(rxshortify.search(yam["unique-id"]).group(1, 2))
    )

    for bl in yam["belongs-to-lineage"]:
        for (s, l) in bl.items():
            if s == "PANGO":
                outy["variant"]["pangolin"] = l
            elif s == "nextstrain":
                outy["variant"]["nextstrain"] = l
            else:
                print(f"unkown {s}:{l}")

    # parse mutations
    defcat = "mut"
    indelscat = (
        "extra"
        if 0 == int(yam["calling-definition"]["probable"]["indels-required"])
        else defcat
    )
    gene = {defcat: None, indelscat: None}
    variants = {defcat: "", indelscat: ""}
    for v in yam["variants"]:
        # TODO depending on phe-genomics issue #9 also handle 'shared' and 'subset'

        cat = defcat
        os = ""

        # actual mutation
        if v["type"] == "SNP" or v["type"] == "MNP":
            os += f"  {v['one-based-reference-position']}: '{v['reference-base']}>{v['variant-base']}'"
        elif v["type"] == "deletion":
            cat = indelscat
            os += f"  {1+int(v['one-based-reference-position'])}: '{ '-' * (len(v['reference-base'])-len(v['variant-base'])) }'"
        elif v["type"] == "insertion":
            cat = indelscat
            os += f"  {1+int(v['one-based-reference-position'])}: '+{ v['variant-base'][(len(v['reference-base'])):] }'"

        # comments
        if "amino-acid-change" in v:
            os += f" # {(v['protein']+':') if 'protein' in v else ''}{v['amino-acid-change']}\n"
        elif "predicted-effect" in v and v["predicted-effect"] == "synonymous":
            os += f" # syn{(' '+v['protein']) if 'protein' in v else ''}\n"
        else:
            os += "\n"

        # gene section header
        if gene[cat] != v["gene"]:
            gene[cat] = v["gene"]
            variants[cat] += f"  # { v['gene'] }\n"

        variants[cat] += os

    # auto derive output name from shortname
    if outname is not None and outname == "":
        outname = f"{outy['variant']['short'].lower()}_mutations.yaml"
    # HACK None, '', string tri-state is a hackish bad idea

    # start writing
    with open(outname, "wt") if outname is not None else sys.stdout as yf:
        print(yaml.dump(outy, sort_keys=False), end="", file=yf)
        for (t, s) in variants.items():
            if 0 == len(s):
                continue
            print(f"{t}:", file=yf)
            print(s, end="", file=yf)

    sys.exit(0)

    """
    unique-id: slinky-antennae
    phe-label: VOC-21JAN-02
    alternate-names:
    - VOC202101/02
    - Brazil 2
    belongs-to-lineage:
    - PANGO: P.1
    description: This variant was first identified in Japan in travellers from Brazil and is associated with Manaus in the Amazonas region assoicated with a severe second wave of COVID-19
    information-sources:
    - https://virological.org/t/genomic-characterisation-of-an-emergent-sars-cov-2-lineage-in-manaus-preliminary-findings/586
    - https://cov-lineages.org/global_report_P.1.html
    variants:
    - gene: S
      one-based-reference-position: 21764
      protein: surface glycoprotein
      reference-base: ATACATG
      type: deletion
      variant-base: A
    """
    """
    variant:
      short: 'BR'
      pangolin: 'P.1'
      nextstrain: '20J/501Y.V3'
      voc: 'VOC-202101/02'
    source:
      - https://virological.org/t/genomic-characterisation-of-an-emergent-sars-cov-2-lineage-in-manaus-preliminary-findings/586
      - doi:10.1101/2021.02.26.21252554
      - https://www.gov.uk/government/publications/investigation-of-novel-sars-cov-2-variant-variant-of-concern-20201201#attachment_4970549
    threshold: 5

    """


if __name__ == "__main__":
    phe2cojac()
