"""Main fetching logic."""

from __future__ import annotations

import hashlib
import json
import threading
from difflib import get_close_matches
from functools import cached_property
from typing import TYPE_CHECKING, Never

import requests

from ._graphql import DYE_QUERY, MICROSCOPE_QUERY, PROTEIN_QUERY, SPECTRUM_QUERY
from .models import (
    Camera,
    DyeResponse,
    Filter,
    Fluorophore,
    Light,
    Microscope,
    MicroscopeResponse,
    Protein,
    ProteinResponse,
    Spectrum,
    SpectrumResponse,
)

if TYPE_CHECKING:
    from collections.abc import Iterable, Mapping


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

    def __init__(self, base_url: str = "https://www.fpbase.org/graphql/"):
        self.base_url = base_url
        self.session = requests.Session()
        self.session.headers.update(
            {"Content-Type": "application/json", "User-Agent": "fpbase-py"}
        )
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
        _ids = self._fluorophore_ids
        try:
            fluor_info = _ids[name.lower()]
        except KeyError as e:
            _raise_with_suggestion(name, _ids, e, "Fluorophore")

        if fluor_info["type"] == "d":
            return self._get_dye_by_id(fluor_info["id"])
        elif fluor_info["type"] == "p":
            return self._get_protein_by_id(fluor_info["id"])
        raise ValueError(f"Invalid fluorophore type {fluor_info['type']!r}")

    def get_protein(self, name: str) -> Protein:
        """Fetch protein by name.

        Examples
        --------
        >>> get_protein("EGFP")
        """
        _ids = self._fluorophore_ids
        try:
            fluor_info = _ids[name.lower()]
        except KeyError as e:
            _raise_with_suggestion(name, _ids, e, "Protein")
        if fluor_info["type"] != "p":
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
        return [item["name"] for item in json.loads(resp)["data"]["microscopes"]]

    def list_filters(self) -> list[str]:
        """List all available filters."""
        return sorted(self._filter_spectrum_ids.keys())

    def list_cameras(self) -> list[str]:
        """List all available cameras."""
        return sorted(self._camera_spectrum_ids.keys())

    def list_lights(self) -> list[str]:
        """List all available lights."""
        return sorted(self._light_spectrum_ids.keys())

    def get_filter(self, name: str) -> Filter:
        """Fetch filter by name."""
        spectrum = self._get_spectrum(name, "Filter")
        if spectrum.owner_filter is None:
            raise ValueError(f"Filter {name!r} not found.")
        return spectrum.owner_filter

    def get_camera(self, name: str) -> Camera:
        """Fetch camera spectrum by name."""
        spectrum = self._get_spectrum(name, "Camera")
        if spectrum.owner_camera is None:
            raise ValueError(f"Camera {name!r} not found.")
        return spectrum.owner_camera

    def get_light(self, name: str) -> Light:
        """Fetch light spectrum by name."""
        spectrum = self._get_spectrum(name, "Light")
        if spectrum.owner_light is None:
            raise ValueError(f"Light {name!r} not found.")
        return spectrum.owner_light

    def _get_spectrum(self, name: str, type_: str) -> Spectrum:
        possibilities: Mapping[str, int] = {
            "Filter": self._filter_spectrum_ids,
            "Light": self._light_spectrum_ids,
            "Camera": self._camera_spectrum_ids,
        }[type_]
        normed = _norm_name(name)
        try:
            filter_id = possibilities[normed]
        except KeyError as e:
            _raise_with_suggestion(name, possibilities, e, type_)

        resp = self._send_query(SPECTRUM_QUERY, {"id": int(filter_id)})
        return SpectrumResponse.model_validate_json(resp).data.spectrum

    # -----------------------------------------------------------

    def _send_query(self, query: str, variables: dict | None = None) -> bytes:
        payload = {"query": query, "variables": variables or {}}
        payload_str = json.dumps(payload, sort_keys=True)  # Convert to JSON string
        # Create a hash
        hashkey = hashlib.md5(payload_str.encode("utf-8")).hexdigest()
        if hashkey not in self._cache:
            data = json.dumps(payload).encode("utf-8")
            response = self.session.post(self.base_url, data=data)
            response.raise_for_status()
            self._cache[hashkey] = response.content
        return self._cache[hashkey]

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


def get_light(name: str) -> Light:
    return FPbaseClient.instance().get_light(name)


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


def list_lights() -> list[str]:
    return FPbaseClient.instance().list_lights()


def _raise_with_suggestion(
    query: str, possibilities: Iterable[str], e: Exception, type_: str
) -> Never:
    """Raise a ValueError with a suggestion if a close match is found."""
    if closest := get_close_matches(query, possibilities, n=1, cutoff=0.5):
        suggest = f" Did you mean {closest[0]!r}?"
    else:
        suggest = ""
    raise ValueError(f"{type_} {query!r} not found.{suggest}") from e
