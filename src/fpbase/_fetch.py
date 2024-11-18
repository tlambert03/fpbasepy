"""Main fetching logic."""

from __future__ import annotations

import hashlib
import json
import threading
from difflib import get_close_matches
from functools import cached_property
from typing import TYPE_CHECKING

import requests

from ._graphql import DYE_QUERY, FILTER_QUERY, MICROSCOPE_QUERY, PROTEIN_QUERY
from .models import (
    DyeResponse,
    FilterSpectrum,
    FilterSpectrumResponse,
    Fluorophore,
    Microscope,
    MicroscopeResponse,
    ProteinResponse,
)

if TYPE_CHECKING:
    from collections.abc import Mapping


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
        >>> get_fluorophore("mturquoise2")
        """
        _ids = self._fluorophore_ids
        if name in _ids:  # direct hit
            fluor_info = _ids[name]
        else:
            try:
                fluor_info = _ids[name.lower()]
            except KeyError as e:
                if closest := get_close_matches(name, _ids, n=1, cutoff=0.5):
                    suggest = f" Did you mean {closest[0]!r}?"
                else:
                    suggest = ""
                raise ValueError(f"Fluorophore {name!r} not found.{suggest}") from e

        if fluor_info["type"] == "d":
            return self._get_dye_by_id(fluor_info["id"])
        elif fluor_info["type"] == "p":
            return self._get_protein_by_id(fluor_info["id"])
        raise ValueError(f"Invalid fluorophore type {fluor_info['type']!r}")

    def get_filter(self, name: str) -> FilterSpectrum:
        """Fetch filter spectrum by name."""
        normed = _norm_name(name)
        try:
            filter_id = self._filter_spectrum_ids[normed]
        except KeyError as e:
            if closest := get_close_matches(
                normed, self._filter_spectrum_ids, n=1, cutoff=0.5
            ):
                suggest = f" Did you mean {closest[0]!r}?"
            else:
                suggest = ""
            raise ValueError(f"Filter {name!r} not found.{suggest}") from e

        resp = self._send_query(FILTER_QUERY, {"id": int(filter_id)})
        return FilterSpectrumResponse.model_validate_json(resp).data.spectrum

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
                lookup[item["name"].lower()] = {"id": item["id"], "type": key[0]}
                lookup[item["slug"]] = {"id": item["id"], "type": key[0]}
                if key == "proteins":
                    lookup[item["id"]] = {"id": item["id"], "type": key[0]}
        return lookup

    @cached_property
    def _filter_spectrum_ids(self) -> Mapping[str, int]:
        resp = self._send_query('{ spectra(category:"F") { id owner { name } } }')
        data: dict = json.loads(resp)["data"]["spectra"]
        return {_norm_name(item["owner"]["name"]): int(item["id"]) for item in data}

    def _get_dye_by_id(self, id: str | int) -> Fluorophore:
        resp = self._send_query(DYE_QUERY, {"id": int(id)})
        return DyeResponse.model_validate_json(resp).data.dye

    def _get_protein_by_id(self, id: str) -> Fluorophore:
        resp = self._send_query(PROTEIN_QUERY, {"id": id})
        return ProteinResponse.model_validate_json(resp).data.protein


def _norm_name(name: str) -> str:
    return name.lower().replace(" ", "-").replace("/", "-")


def get_microscope(id: str = "i6WL2W") -> Microscope:
    return FPbaseClient.instance().get_microscope(id)


def get_fluorophore(name: str) -> Fluorophore:
    return FPbaseClient.instance().get_fluorophore(name)


def get_filter(name: str) -> FilterSpectrum:
    return FPbaseClient.instance().get_filter(name)
