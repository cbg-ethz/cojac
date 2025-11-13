from .mut_parser import to_single
import click
import yaml


@click.command(
    help="display a voc YAML as a list of single mutations, e.g. for querying covSPECTRUM, V-pipe scout, etc.",
)
@click.option(
    "--deletions/--no-deletions",
    "--del/--no-del",
    "-d/-D",
    "deletions",
    default=False,
    help="also includes deletions",
)
@click.argument("voc", nargs=1)
def parser_voc2single(voc, deletions):

    with open(voc, "rt") as y:
        mutlist = yaml.load(y, Loader=yaml.Loader)["mut"]

    print(", ".join(to_single(dict(sorted(mutlist.items())), nodels=not deletions)))
