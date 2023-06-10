import click
from ._version import __version__

CONTEXT_SETTINGS = dict(help_option_names=["-h", "--help"])

from .cooc_colourmut import cooc_colourmut
from .cooc_curate import cooc_curate
from .cooc_mutbamscan import cooc_mutbamscan
from .cooc_pubmut import cooc_pubmut
from .cooc_tabmut import cooc_tabmut
from .phe2cojac import phe2cojac
from .sig_generate import sig_generate
from .generate_sigs_nextstrains import generate_sigs_nextstrains


@click.group(context_settings=CONTEXT_SETTINGS)
@click.version_option(__version__)
def cli():
    #  auto_envvar_prefix='...' for fetching options (such as passwords) from environment.
    pass


cli.add_command(cooc_colourmut)
cli.add_command(cooc_curate)
cli.add_command(cooc_mutbamscan)
cli.add_command(cooc_pubmut)
cli.add_command(cooc_tabmut)
cli.add_command(phe2cojac)
cli.add_command(sig_generate)
cli.add_command(generate_sigs_nextstrains)
