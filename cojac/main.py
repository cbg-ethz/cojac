import click

from .cooc_colourmut import cooc_colourmut
from .cooc_curate import cooc_curate
from .cooc_mutbamscan import cooc_mutbamscan
from .cooc_pubmut import cooc_pubmut
from .cooc_tabmut import cooc_tabmut
from .phe2cojac import phe2cojac


@click.group()
def cli():
    pass


cli.add_command(cooc_colourmut)
cli.add_command(cooc_curate)
cli.add_command(cooc_mutbamscan)
cli.add_command(cooc_pubmut)
cli.add_command(cooc_tabmut)
cli.add_command(phe2cojac)
