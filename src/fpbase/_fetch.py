"""Main fetching logic."""

from __future__ import annotations

import hashlib
import json
import threading
from difflib import get_close_matches
from functools import cached_property
from typing import TYPE_CHECKING, Any, Final

import requests

from ._graphql import DYE_QUERY, MICROSCOPE_QUERY, PROTEIN_QUERY, SPECTRUM_QUERY
from .models import (
    Camera,
    DyeResponse,
    Filter,
    Fluorophore,
    LightSource,
    Microscope,
    MicroscopeResponse,
    Protein,
    ProteinResponse,
    Spectrum,
    SpectrumResponse,
)

if TYPE_CHECKING:
    from collections.abc import Mapping

FPBASE_URL: Final = "https://www.fpbase.org/graphql/"
_HEADERS = {"Content-Type": "application/json", "User-Agent": "fpbase-py"}


class FPbaseClient:
    __instance: FPbaseClient | None = None
    __lock: threading.Lock = threading.Lock()

    @classmethod
    def instance(cls) -> FPbaseClient:
        if cls.__instance is None:
            with cls.__lock:
                if cls.__instance is None:  # Double-checked locking
                    cls.__instance = cls()
        return cls.__instance

    def __init__(self, base_url: str = FPBASE_URL):
        self.base_url = base_url
        self.session = requests.Session()
        self.session.headers.update(_HEADERS)
        self._cache: dict[str, bytes] = {}

    def get_microscope(self, id: str = "i6WL2W") -> Microscope:
        """Get microscope by ID.

        Examples
        --------
        >>> get_microscope("i6WL2W")
        """
        resp = self._send_query(MICROSCOPE_QUERY, {"id": id})
        return MicroscopeResponse.model_validate_json(resp).data.microscope

    def get_fluorophore(self, name: str) -> Fluorophore:
        """Fetch fluorophore by name, slug, or ID.

        Examples
        --------
        >>> get_fluorophore("mTurquoise2")
        >>> get_fluorophore("Alexa Fluor 488")
        """
        fluor_info = _get_or_raise_suggestion(
            name, self._fluorophore_ids, "Fluorophore"
        )

        if fluor_info["type"] == "d":
            return self._get_dye_by_id(fluor_info["id"])
        elif fluor_info["type"] == "p":
            return self._get_protein_by_id(fluor_info["id"])
        raise ValueError(  # pragma: no cover
            f"Invalid fluorophore type {fluor_info['type']!r}"
        )

    def get_protein(self, name: str) -> Protein:
        """Fetch protein by name.

        Examples
        --------
        >>> get_protein("EGFP")
        """
        fluor_info = _get_or_raise_suggestion(name, self._fluorophore_ids, "Protein")
        if fluor_info["type"] != "p":  # pragma: no cover
            raise ValueError(f"Protein {name!r} not found.")
        return self._get_protein_by_id(fluor_info["id"])

    def list_proteins(self) -> list[str]:
        """List all available proteins."""
        return sorted(
            {
                info["name"]
                for info in self._fluorophore_ids.values()
                if info["type"] == "p"
            }
        )

    def list_dyes(self) -> list[str]:
        """List all available dyes."""
        return sorted(
            {
                info["name"]
                for info in self._fluorophore_ids.values()
                if info["type"] == "d"
            }
        )

    def list_fluorophores(self) -> list[str]:
        """List all available fluorophores."""
        return sorted({info["name"] for info in self._fluorophore_ids.values()})

    def list_microscopes(self) -> list[str]:
        """List all available microscopes."""
        resp = self._send_query("{ microscopes { id name } }")
        return [item["id"] for item in json.loads(resp)["data"]["microscopes"]]

    def list_filters(self) -> list[str]:
        """List all available filters."""
        return sorted(self._filter_spectrum_ids.keys())

    def list_cameras(self) -> list[str]:
        """List all available cameras."""
        return sorted(self._camera_spectrum_ids.keys())

    def list_light_sources(self) -> list[str]:
        """List all available lights."""
        return sorted(self._light_spectrum_ids.keys())

    def get_filter(self, name: str) -> Filter:
        """Fetch filter by name."""
        spectrum = self._get_spectrum(name, "Filter")
        if spectrum.owner_filter is None:  # pragma: no cover
            raise ValueError(f"Filter {name!r} not found.")
        return spectrum.owner_filter

    def get_camera(self, name: str) -> Camera:
        """Fetch camera spectrum by name."""
        spectrum = self._get_spectrum(name, "Camera")
        if spectrum.owner_camera is None:  # pragma: no cover
            raise ValueError(f"Camera {name!r} not found.")
        return spectrum.owner_camera

    def get_light_source(self, name: str) -> LightSource:
        """Fetch light spectrum by name."""
        spectrum = self._get_spectrum(name, "Light")
        if spectrum.owner_light is None:  # pragma: no cover
            raise ValueError(f"Light {name!r} not found.")
        return spectrum.owner_light

    def _get_spectrum(self, name: str, type_: str) -> Spectrum:
        possibilities: Mapping[str, int] = {
            "Filter": self._filter_spectrum_ids,
            "Light": self._light_spectrum_ids,
            "Camera": self._camera_spectrum_ids,
        }[type_]
        normed = _norm_name(name)
        filter_id = _get_or_raise_suggestion(normed, possibilities, type_)

        resp = self._send_query(SPECTRUM_QUERY, {"id": int(filter_id)})
        return SpectrumResponse.model_validate_json(resp).data.spectrum

    # -----------------------------------------------------------

    def _send_query(self, query: str, variables: dict | None = None) -> bytes:
        # Create a hash
        if (key := _hashargs(self.base_url, query, variables)) not in self._cache:
            payload = {"query": query, "variables": variables or {}}
            data = json.dumps(payload).encode("utf-8")
            response = self.session.post(self.base_url, data=data)
            response.raise_for_status()
            self._cache[key] = response.content
        return self._cache[key]

    @cached_property
    def _fluorophore_ids(self) -> dict[str, dict[str, str]]:
        """Return a lookup table of fluorophore {name: {id: ..., type: ...}}."""
        resp = self._send_query("{ dyes { id name slug } proteins { id name slug } }")
        data: dict[str, list[dict[str, str]]] = json.loads(resp)["data"]
        lookup: dict[str, dict[str, str]] = {}
        for key in ["dyes", "proteins"]:
            for item in data[key]:
                lookup[item["name"].lower()] = {**item, "type": key[0]}
                lookup[item["slug"]] = {**item, "type": key[0]}
                if key == "proteins":
                    lookup[item["id"].lower()] = {**item, "type": key[0]}
        return lookup

    @cached_property
    def _filter_spectrum_ids(self) -> Mapping[str, int]:
        return self._get_spectrum_ids("F")

    @cached_property
    def _light_spectrum_ids(self) -> Mapping[str, int]:
        return self._get_spectrum_ids("L")

    @cached_property
    def _camera_spectrum_ids(self) -> Mapping[str, int]:
        return self._get_spectrum_ids("C")

    def _get_spectrum_ids(self, key: str) -> dict[str, int]:
        query = f'{{ spectra(category: "{key}") {{ id owner {{ name }} }} }}'
        resp = self._send_query(query)
        data = json.loads(resp)["data"]["spectra"]
        return {_norm_name(item["owner"]["name"]): int(item["id"]) for item in data}

    def _get_dye_by_id(self, id: str | int) -> Fluorophore:
        resp = self._send_query(DYE_QUERY, {"id": int(id)})
        return DyeResponse.model_validate_json(resp).data.dye

    def _get_protein_by_id(self, id: str) -> Protein:
        resp = self._send_query(PROTEIN_QUERY, {"id": id})
        return ProteinResponse.model_validate_json(resp).data.protein


def _norm_name(name: str) -> str:
    return name.lower().replace(" ", "-").replace("/", "-")


def get_microscope(id: str = "i6WL2W") -> Microscope:
    return FPbaseClient.instance().get_microscope(id)


def get_fluorophore(name: str) -> Fluorophore:
    return FPbaseClient.instance().get_fluorophore(name)


def get_filter(name: str) -> Filter:
    return FPbaseClient.instance().get_filter(name)


def get_camera(name: str) -> Camera:
    return FPbaseClient.instance().get_camera(name)


def get_light_source(name: str) -> LightSource:
    return FPbaseClient.instance().get_light_source(name)


def get_protein(name: str) -> Protein:
    return FPbaseClient.instance().get_protein(name)


def list_proteins() -> list[str]:
    return FPbaseClient.instance().list_proteins()


def list_dyes() -> list[str]:
    return FPbaseClient.instance().list_dyes()


def list_fluorophores() -> list[str]:
    return FPbaseClient.instance().list_fluorophores()


def list_microscopes() -> list[str]:
    return FPbaseClient.instance().list_microscopes()


def list_filters() -> list[str]:
    return FPbaseClient.instance().list_filters()


def list_cameras() -> list[str]:
    return FPbaseClient.instance().list_cameras()


def list_light_sources() -> list[str]:
    return FPbaseClient.instance().list_light_sources()


def _get_or_raise_suggestion(
    query: str, possibilities: Mapping[str, Any], type_: str
) -> Any:
    """Raise a ValueError with a suggestion if a close match is found."""
    try:
        return possibilities[query.lower()]
    except KeyError as e:
        if closest := get_close_matches(query, possibilities, n=1, cutoff=0.5):
            suggest = f" Did you mean {closest[0]!r}?"
        else:  # pragma: no cover
            suggest = ""
        raise ValueError(f"{type_} {query!r} not found.{suggest}") from e


_RESPONSE_CACHE: dict[str, dict] = {}


def graphql_query(
    query: str,
    variables: dict | None = None,
    *,
    session: requests.Session | None = None,
) -> dict[str, Any]:
    """Send a generic GraphQL query to the FPbase API.

    See docs and test out queries at https://www.fpbase.org/graphql/

    Parameters
    ----------
    query : str
        A graphql query string. For example, "{ proteins { name } }"
    variables : dict | None, optional
        If the query requires variables, pass them here, by default None
    session : requests.Session | None, optional
        Optionally pass a requests session, by default, will create a new session.

    Returns
    -------
    dict[str, Any]
        JSON response from the API deserialized into a Python dictionary.

    Examples
    --------
    # get all protein names and sequences
    >>> data = fpbase.graphql_query("{proteins { name seq } }")

    # get specific fields for a specific protein
    # note that single/double quotes are NOT interchangeable here
    >>> data = fpbase.graphql_query('{protein(id: "R9NL8") { name seq } }')

    # get optical configs for a microscope, using a variable
    >>> q = "query getScope($id: String!){ microscope(id: $id){ name opticalConfigs {name} } }"
    >>> data = fpbase.graphql_query(q, {"id": "i6WL2WdgcDMgJYtPrpZcaJ"})
    """  # noqa: E501
    url = FPBASE_URL
    if (key := _hashargs(url, query, variables)) not in _RESPONSE_CACHE:
        data_bytes = _fetch_query(query, variables, session=session, url=url)
        _RESPONSE_CACHE[key] = json.loads(data_bytes)
    return _RESPONSE_CACHE[key]


def _hashargs(*args: str | dict | None | tuple) -> str:
    hasher = hashlib.md5()
    for arg in args:
        if isinstance(arg, dict):
            arg = tuple(sorted(arg.items()))
        hasher.update(str(arg).encode("utf-8"))
    return hasher.hexdigest()


def _fetch_query(
    query: str,
    variables: dict | None = None,
    *,
    session: requests.Session | None = None,
    url: str = FPBASE_URL,
) -> bytes:
    payload = {"query": query, "variables": variables or {}}
    data = json.dumps(payload).encode("utf-8")
    post = requests.post if session is None else session.post
    response = post(url, data=data, headers=_HEADERS)
    response.raise_for_status()
    return response.content
