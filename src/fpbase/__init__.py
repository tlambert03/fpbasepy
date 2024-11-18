"""Python wrapper for FPBase API."""

from importlib.metadata import PackageNotFoundError, version

try:
    __version__ = version("fpbasepy")
except PackageNotFoundError:
    __version__ = "uninstalled"
__author__ = "Talley Lambert"
__email__ = "talley.lambert@gmail.com"

from ._fetch import FPbaseClient, get_filter, get_fluorophore, get_microscope

__all__ = [
    "FPbaseClient",
    "get_filter",
    "get_fluorophore",
    "get_microscope",
]
