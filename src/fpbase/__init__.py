"""Python wrapper for FPBase API"""

from importlib.metadata import PackageNotFoundError, version

try:
    __version__ = version("fpbasepy")
except PackageNotFoundError:
    __version__ = "uninstalled"
__author__ = "Talley Lambert"
__email__ = "talley.lambert@gmail.com"
