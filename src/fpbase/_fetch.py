"""Main fetching logic."""

from __future__ import annotations

import hashlib
import json
import os
import sys
import tempfile
import threading
import time
import warnings
from difflib import get_close_matches
from functools import cached_property
from pathlib import Path
from typing import TYPE_CHECKING, Any, Final
from urllib.parse import urlsplit

import requests
from requests.adapters import HTTPAdapter
from requests.auth import AuthBase
from urllib3.util.retry import Retry

from . import __version__
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
ISSUES_URL: Final = "https://github.com/tlambert03/fpbasepy/issues"
# the version lets the server tell client versions apart (e.g. which support API keys)
USER_AGENT: Final = f"fpbase-py/{__version__} (+https://github.com/tlambert03/fpbasepy)"
_HEADERS = {"Content-Type": "application/json", "User-Agent": USER_AGENT}
# header the server can use to send a message to users (e.g. upcoming API changes)
NOTICE_HEADER: Final = "X-FPbase-Notice"
_LOCAL_HOSTS = {"localhost", "127.0.0.1", "::1"}
# server errors are retried only once: a 503/504 is usually a request that ran past
# the server's 30s timeout, and resending it just ties up another server worker
_RETRY_ONCE_STATUSES: Final = frozenset({502, 503, 504})


class _Retry(Retry):
    """Retry throttling (429) patiently, but a server error only once, after a pause."""

    def is_retry(
        self, method: str, status_code: int, has_retry_after: bool = False
    ) -> bool:
        if status_code in _RETRY_ONCE_STATUSES and any(
            h.status in _RETRY_ONCE_STATUSES for h in self.history
        ):
            return False
        return super().is_retry(method, status_code, has_retry_after)

    def get_backoff_time(self) -> float:
        backoff = super().get_backoff_time()
        if self.history and self.history[-1].status in _RETRY_ONCE_STATUSES:
            # give a restarting or overloaded server a moment
            return max(backoff, 5 * self.backoff_factor)
        return backoff


# wait and retry when throttled (429) or the server is briefly unavailable,
# honoring the server's Retry-After header
_RETRY = _Retry(
    total=5,
    backoff_factor=1,
    status_forcelist=(429, 502, 503, 504),
    allowed_methods=None,  # retry POST too: graphql queries are read-only
    respect_retry_after_header=True,
    raise_on_status=False,  # return the last response, so raise_for_status() raises
)
# how long lookup tables (e.g. all protein names) are cached on disk
DISK_CACHE_TTL: Final = 24 * 60 * 60  # seconds
# microscopes are edited by their owners, so they're cached for less time
MICROSCOPE_CACHE_TTL: Final = 60 * 60  # seconds


class FPbaseWarning(UserWarning):
    """A message from the FPbase server."""


class _ApiKeyAuth(AuthBase):
    """Send the API key as a bearer token, but only over HTTPS (or to localhost)."""

    def __init__(self, api_key: str) -> None:
        self.api_key = api_key

    def __call__(self, r: requests.PreparedRequest) -> requests.PreparedRequest:
        url = urlsplit(r.url or "")
        if url.scheme == "https" or url.hostname in _LOCAL_HOSTS:
            r.headers["Authorization"] = f"Bearer {self.api_key}"
        return r


_shown_notices: set[str] = set()


def _warn_notice(response: requests.Response, *args: Any, **kwargs: Any) -> None:
    if (notice := response.headers.get(NOTICE_HEADER)) and notice not in _shown_notices:
        _shown_notices.add(notice)
        warnings.warn(notice, FPbaseWarning, stacklevel=2)


def _new_session(api_key: str | None = None) -> requests.Session:
    session = requests.Session()
    session.headers.update(_HEADERS)
    adapter = HTTPAdapter(max_retries=_RETRY)
    session.mount("https://", adapter)
    session.mount("http://", adapter)
    if api_key := api_key or os.environ.get("FPBASE_API_KEY"):
        session.auth = _ApiKeyAuth(api_key)
    session.hooks["response"].append(_warn_notice)
    return session


def _raise_for_status(response: requests.Response) -> None:
    """Like `response.raise_for_status()`, with an explanation of common errors."""
    if response.ok:
        return
    if response.headers.get("cf-mitigated") == "challenge":
        msg = (
            "FPbase's bot protection blocked this request, which usually happens "
            "on cloud servers (e.g. AWS, Azure). If this persists, please open an "
            f"issue at {ISSUES_URL}"
        )
    else:
        msg = f"FPbase returned HTTP {response.status_code}"
        if response.status_code == 429:
            msg += " (rate limit exceeded, even after waiting and retrying)"
        if detail := _error_detail(response):
            msg += f": {detail}"
    raise requests.HTTPError(msg, response=response)


def _error_detail(response: requests.Response) -> str | None:
    """Return the error message in a JSON (DRF or GraphQL) error response."""
    try:
        data = response.json()
    except ValueError:
        return None
    if not isinstance(data, dict):
        return None
    if isinstance(detail := data.get("detail"), str):
        return detail
    if (errors := data.get("errors")) and isinstance(errors, list):
        first = errors[0]
        return str(first.get("message", first) if isinstance(first, dict) else first)
    return None


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

    def __init__(self, base_url: str = FPBASE_URL, api_key: str | None = None):
        """Create a client.

        Parameters
        ----------
        base_url : str
            URL of the FPbase GraphQL endpoint.
        api_key : str | None
            FPbase API key. Defaults to the `FPBASE_API_KEY` environment variable.
        """
        self.base_url = base_url
        self.session = _new_session(api_key)
        self._cache: dict[str, bytes] = {}

    def get_microscope(self, id: str = "i6WL2W") -> Microscope:
        """Get microscope by ID.

        Examples
        --------
        >>> get_microscope("i6WL2W")
        """
        # on disk too: a script run repeatedly shouldn't re-download a large microscope
        resp = self._send_query(
            MICROSCOPE_QUERY, {"id": id}, persist=True, ttl=MICROSCOPE_CACHE_TTL
        )
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

    def _send_query(
        self,
        query: str,
        variables: dict | None = None,
        *,
        persist: bool = False,
        ttl: float = DISK_CACHE_TTL,
    ) -> bytes:
        """Send query, caching in memory (and on disk for `ttl` s, if `persist`)."""
        if (key := _hashargs(self.base_url, query, variables)) not in self._cache:
            content = _read_disk_cache(key, ttl) if persist else None
            if content is None:
                payload = {"query": query, "variables": variables or {}}
                data = json.dumps(payload).encode("utf-8")
                response = self.session.post(self.base_url, data=data)
                _raise_for_status(response)
                content = response.content
                if persist and "errors" not in json.loads(content):
                    _write_disk_cache(key, content)
            self._cache[key] = content
        return self._cache[key]

    @cached_property
    def _fluorophore_ids(self) -> dict[str, dict[str, str]]:
        """Return a lookup table of fluorophore {name: {id: ..., type: ...}}."""
        resp = self._send_query(
            "{ dyes { id name slug } proteins { id name slug } }", persist=True
        )
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
        resp = self._send_query(query, persist=True)
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
        Optionally pass a requests session. By default, a shared session is used
        that retries throttled (429) and briefly-unavailable requests.

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


def _hashargs(*args: str | dict | tuple | None) -> str:
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
    session = session or FPbaseClient.instance().session
    response = session.post(url, data=data, headers=_HEADERS)
    _raise_for_status(response)
    return response.content


def _cache_dir() -> Path:
    if env := os.environ.get("FPBASE_CACHE_DIR"):
        return Path(env)
    if sys.platform == "win32":
        base = os.environ.get("LOCALAPPDATA", Path.home() / "AppData" / "Local")
    elif sys.platform == "darwin":
        base = Path.home() / "Library" / "Caches"
    else:
        base = os.environ.get("XDG_CACHE_HOME", Path.home() / ".cache")
    return Path(base) / "fpbase"


def _read_disk_cache(key: str, ttl: float = DISK_CACHE_TTL) -> bytes | None:
    path = _cache_dir() / f"{key}.json"
    try:
        if time.time() - path.stat().st_mtime < ttl:
            return path.read_bytes()
    except OSError:
        pass
    return None


def _write_disk_cache(key: str, content: bytes) -> None:
    # write atomically: many processes may start at once (e.g. jobs on a cluster)
    try:
        cache_dir = _cache_dir()
        cache_dir.mkdir(parents=True, exist_ok=True)
        with tempfile.NamedTemporaryFile(dir=cache_dir, delete=False) as f:
            f.write(content)
        os.replace(f.name, cache_dir / f"{key}.json")
    except OSError:
        pass  # the disk cache is best-effort (e.g. read-only home directory)
