"""Python wrapper for FPBase API."""

from importlib.metadata import PackageNotFoundError, version

try:
    __version__ = version("fpbasepy")
except PackageNotFoundError:  # pragma: no cover
    __version__ = "uninstalled"
__author__ = "Talley Lambert"
__email__ = "talley.lambert@gmail.com"

from . import models
from ._fetch import (
    FPbaseClient,
    get_camera,
    get_filter,
    get_fluorophore,
    get_light_source,
    get_microscope,
    get_protein,
    graphql_query,
    list_cameras,
    list_dyes,
    list_filters,
    list_fluorophores,
    list_light_sources,
    list_microscopes,
    list_proteins,
)

__all__ = [
    "FPbaseClient",
    "get_camera",
    "get_filter",
    "get_fluorophore",
    "get_light_source",
    "get_microscope",
    "get_protein",
    "graphql_query",
    "list_cameras",
    "list_dyes",
    "list_filters",
    "list_fluorophores",
    "list_light_sources",
    "list_microscopes",
    "list_proteins",
    "models",
]
