"""Python wrapper for FPBase API."""

from importlib.metadata import PackageNotFoundError, version

try:
    __version__ = version("fpbasepy")
except PackageNotFoundError:
    __version__ = "uninstalled"
__author__ = "Talley Lambert"
__email__ = "talley.lambert@gmail.com"

from . import models
from ._fetch import (
    FPbaseClient,
    get_filter,
    get_fluorophore,
    get_microscope,
    get_protein,
    list_dyes,
    list_filters,
    list_fluorophores,
    list_microscopes,
    list_proteins,
)

__all__ = [
    "FPbaseClient",
    "get_filter",
    "get_fluorophore",
    "get_microscope",
    "get_protein",
    "list_dyes",
    "list_filters",
    "list_fluorophores",
    "list_microscopes",
    "list_proteins",
    "models",
]
