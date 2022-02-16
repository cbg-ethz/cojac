from .main import cli
from .cooc_colourmut import main as cooc_colourmut
from .cooc_curate import main as cooc_curate
from .cooc_mutbamscan import main as cooc_mutbamscan
from .cooc_pubmut import main as cooc_pubmut
from .cooc_tabmut import main as cooc_tabmut
from .phe2cojac import main as phe2cojac


__all__ = [
    "cli",
    "cooc_colourmut",
    "cooc_curate",
    "cooc_mutbamscan",
    "cooc_pubmut",
    "cooc_tabmut",
    "phe2cojac",
]
