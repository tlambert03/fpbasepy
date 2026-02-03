"""Main fetching logic."""

from __future__ import annotations

import hashlib
import json
import re
import threading
from difflib import get_close_matches
from functools import cached_property
from typing import TYPE_CHECKING, Any, Final

import requests

from ._graphql import (
    DYE_QUERY,
    MICROSCOPE_QUERY,
    PROTEIN_QUERY,
    SPECTRUM_QUERY,
    build_multiple_query,
)
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
    from collections.abc import Mapping, Sequence

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

    def get_multiple_proteins(self, names: Sequence[str]) -> dict[str, Protein | None]:
        """Fetch multiple proteins by name in a single request.

        Parameters
        ----------
        names : Sequence[str]
            Protein names, slugs, or IDs to fetch

        Returns
        -------
        dict[str, Protein | None]
            Mapping of input name -> Protein object (or None if not found)

        Examples
        --------
        >>> proteins = client.get_multiple_proteins(["EGFP", "mCherry", "mTurquoise2"])
        >>> proteins["EGFP"].default_state.ex_max
        488.0
        """
        if not names:
            return {}

        # Resolve names to IDs using existing _fluorophore_ids cache
        items = {}
        name_to_alias = {}
        for name in names:
            try:
                fluor_info = _get_or_raise_suggestion(
                    name, self._fluorophore_ids, "Protein"
                )
                if fluor_info["type"] == "p":
                    # Use name as alias, create safe alias if needed
                    alias = _make_safe_alias(name)
                    items[alias] = {"id": fluor_info["id"]}
                    name_to_alias[name] = alias
            except ValueError:
                # Name not found, will return None for this entry
                pass

        if not items:
            return dict.fromkeys(names)

        # Extract fields from PROTEIN_QUERY template
        fields = _extract_fields_from_query(PROTEIN_QUERY)

        # Build and execute multi-item query
        query = build_multiple_query("protein", items, fields)
        resp = self._send_query(query)
        data = json.loads(resp)["data"]

        # Map aliases back to original names and parse responses
        results: dict[str, Protein | None] = {}

        for name in names:
            alias = name_to_alias.get(name)
            if alias and alias in data:
                protein_data = data[alias]
                if protein_data is None:
                    results[name] = None
                else:
                    # Wrap in response structure for validation
                    wrapped = {"data": {"protein": protein_data}}
                    results[name] = ProteinResponse.model_validate(wrapped).data.protein
            else:
                # Name wasn't found in cache
                results[name] = None

        return results

    def get_multiple_dyes(self, names: Sequence[str]) -> dict[str, Fluorophore | None]:
        """Fetch multiple dyes by name in a single request.

        Parameters
        ----------
        names : Sequence[str]
            Dye names or slugs to fetch

        Returns
        -------
        dict[str, Fluorophore | None]
            Mapping of input name -> Fluorophore object (or None if not found)

        Examples
        --------
        >>> dyes = client.get_multiple_dyes(["DAPI", "Hoechst 33342"])
        """
        if not names:
            return {}

        # Resolve names to IDs using existing _fluorophore_ids cache
        items = {}
        name_to_alias = {}
        for name in names:
            try:
                fluor_info = _get_or_raise_suggestion(
                    name, self._fluorophore_ids, "Dye"
                )
                if fluor_info["type"] == "d":
                    # Use name as alias, create safe alias if needed
                    alias = _make_safe_alias(name)
                    items[alias] = {"id": int(fluor_info["id"])}
                    name_to_alias[name] = alias
            except ValueError:
                # Name not found, will return None for this entry
                pass

        if not items:
            return dict.fromkeys(names)

        # Extract fields from DYE_QUERY template
        fields = _extract_fields_from_query(DYE_QUERY)

        # Build and execute multi-item query
        query = build_multiple_query("dye", items, fields)
        resp = self._send_query(query)
        data = json.loads(resp)["data"]

        # Map aliases back to original names and parse responses
        results: dict[str, Fluorophore | None] = {}

        for name in names:
            alias = name_to_alias.get(name)
            if alias and alias in data:
                dye_data = data[alias]
                if dye_data is None:
                    results[name] = None
                else:
                    # Wrap in response structure for validation
                    wrapped = {"data": {"dye": dye_data}}
                    results[name] = DyeResponse.model_validate(wrapped).data.dye
            else:
                # Name wasn't found in cache
                results[name] = None

        return results

    def get_multiple_microscopes(
        self, ids: Sequence[str]
    ) -> dict[str, Microscope | None]:
        """Fetch multiple microscopes by ID in a single request.

        Parameters
        ----------
        ids : Sequence[str]
            Microscope IDs to fetch

        Returns
        -------
        dict[str, Microscope | None]
            Mapping of ID -> Microscope object (or None if not found)

        Examples
        --------
        >>> scopes = client.get_multiple_microscopes(["i6WL2W"])
        """
        if not ids:
            return {}

        # Build items dict with safe aliases
        items = {}
        id_to_alias = {}
        for microscope_id in ids:
            alias = _make_safe_alias(microscope_id)
            items[alias] = {"id": microscope_id}
            id_to_alias[microscope_id] = alias

        # Extract fields from MICROSCOPE_QUERY template
        fields = _extract_fields_from_query(MICROSCOPE_QUERY)

        # Build and execute multi-item query
        query = build_multiple_query("microscope", items, fields)
        resp = self._send_query(query)
        data = json.loads(resp)["data"]

        # Map aliases back to original IDs and parse responses
        results: dict[str, Microscope | None] = {}

        for microscope_id in ids:
            alias = id_to_alias[microscope_id]
            if alias in data:
                microscope_data = data[alias]
                if microscope_data is None:
                    results[microscope_id] = None
                else:
                    # Wrap in response structure for validation
                    wrapped = {"data": {"microscope": microscope_data}}
                    results[microscope_id] = MicroscopeResponse.model_validate(
                        wrapped
                    ).data.microscope
            else:
                results[microscope_id] = None

        return results


def _norm_name(name: str) -> str:
    return name.lower().replace(" ", "-").replace("/", "-")


def _make_safe_alias(name: str) -> str:
    """Convert a name to a safe GraphQL alias."""
    alias = re.sub(r"[^a-zA-Z0-9]", "_", name)
    if not alias or not alias[0].isalpha():
        alias = f"item_{alias}"
    return alias.lower()


def _extract_fields_from_query(query_template: str) -> str:
    """Extract the fields portion from a single-item query template.

    Parameters
    ----------
    query_template : str
        A GraphQL query template like PROTEIN_QUERY

    Returns
    -------
    str
        The fields portion (everything between the innermost { and })
    """
    # Find the last occurrence of the query name (protein, dye, microscope)
    # and extract everything between its opening { and closing }
    # This works for nested queries by finding the innermost query
    lines = query_template.strip().split("\n")
    # Skip the query declaration line and get the body
    start_idx = None
    brace_count = 0
    fields_lines = []

    for i, line in enumerate(lines):
        if "{" in line and start_idx is None:
            # Found first opening brace (after query declaration)
            start_idx = i + 1
            brace_count = line.count("{") - line.count("}")
            continue

        if start_idx is not None:
            brace_count += line.count("{") - line.count("}")
            if brace_count > 0:
                fields_lines.append(line)
            else:
                # Reached the closing brace
                break

    # Clean up the fields: remove leading/trailing whitespace and closing braces
    fields = "\n".join(fields_lines).strip()
    # Remove the outer query wrapper (protein(id: $id) {...})
    # We want just the inner fields
    if "(" in fields.split("{")[0]:
        # Remove everything before the first {
        fields = "{".join(fields.split("{")[1:])
        # Remove the last }
        fields = "}".join(fields.rsplit("}", 1)[:-1])

    return fields.strip()


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


def get_multiple_proteins(names: Sequence[str]) -> dict[str, Protein | None]:
    return FPbaseClient.instance().get_multiple_proteins(names)


def get_multiple_dyes(names: Sequence[str]) -> dict[str, Fluorophore | None]:
    return FPbaseClient.instance().get_multiple_dyes(names)


def get_multiple_microscopes(ids: Sequence[str]) -> dict[str, Microscope | None]:
    return FPbaseClient.instance().get_multiple_microscopes(ids)


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
