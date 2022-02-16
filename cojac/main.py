import click

from .cooc_colourmut import main as cooc_colourmut_script
from .cooc_curate import cooc_curate
from .cooc_mutbamscan import main as cooc_mutbamscan_script
from .cooc_pubmut import main as cooc_pubmut_script
from .cooc_tabmut import main as cooc_tabmut_script
from .phe2cojac import phe2cojac


@click.group()
def cli():
    pass


# TODO: make all subcommands also use click
@cli.command()
def cooc_colourmut():
    cooc_colourmut_script()


@cli.command()
def cooc_mutbamscan():
    cooc_mutbamscan_script()


@cli.command()
def cooc_pubmut():
    cooc_pubmut_script()


@cli.command()
def cooc_tabmut():
    cooc_tabmut_script()


cli.add_command(cooc_curate)
cli.add_command(phe2cojac)
