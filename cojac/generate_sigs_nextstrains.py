#!/usr/bin/env python3

import os
import json
import requests
from datetime import datetime
import yaml
import click


def load_dataset(url):
    response = requests.get(url)
    full_data = json.loads(response.content)

    return full_data


def prepare_header(data, url):
    return {
        "variant": {
            "nextstrain": data["nextstrainClade"],
            "pangolin": data["lineage"],
            "short": data["lineage"].lower().replace(".", "_"),
            "reference": {
                "address": url,
                "accessed_at": datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
            },
        }
    }


def prepare_mut(data):
    mut_dict = {}
    if data["nucSubstitutions"] != [""]:
        for item in data["nucSubstitutions"]:
            start_letter = item[0]
            end_letter = item[-1]
            number = int(item[1:-1])
            mut_dict[number] = f"{start_letter}>{end_letter}"
        return {"mut": mut_dict}
    else:
        return ""


def prepare_del(data):
    del_dict = {}
    if data["nucDeletions"] != [""]:
        for item in data["nucDeletions"]:
            if "-" in item:
                start, end = map(int, item.split("-"))
                key = start
                value = str((end - start) * "-")
            else:
                key = int(item)
                value = "-"
            del_dict[key] = value
        return {"del": del_dict}
    else:
        return ""


def prepare_yaml(data, url):
    header = prepare_header(data, url)
    mut = prepare_mut(data)
    del_section = prepare_del(data)

    output_data = header

    if mut != "":
        output_data["mut"] = mut["mut"]

    if del_section != "":
        output_data["del"] = del_section["del"]

    return output_data


def process_dataset(full_data, outdir, url):
    if not os.path.exists(outdir):
        os.makedirs(outdir)
    for data in full_data.values():
        output_data = prepare_yaml(data, url)
        with open(f"{outdir}/{output_data['variant']['short']}.yaml", "w") as yaml_file:
            yaml.dump(output_data, yaml_file, sort_keys=False)


@click.command(
    help="Generating a list of variants from nextstrain",
    epilog="This tool fetchs a JSON from Github",
)
@click.option(
    "-o",
    "--outdir",
    required=False,
    default="voc_nextstrain",
    type=str,
    help="The output directory for the YAML files",
)
@click.option(
    "-u",
    "--url",
    metavar="URL",
    required=False,
    default="https://raw.githubusercontent.com/corneliusroemer/pango-sequences/main/data/pango-consensus-sequences_summary.json",
    type=str,
    help="url to fetch the JSON from",
)
def generate_sigs_nextstrains(outdir, url):
    full_data = load_dataset(url)
    process_dataset(full_data, outdir, url)
