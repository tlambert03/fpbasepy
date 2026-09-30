"""Tests of client behavior (retries, caching, API keys...), against a local server."""

from __future__ import annotations

import json
import threading
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
from typing import TYPE_CHECKING, Any

import pytest
import requests
from urllib3.util.retry import RequestHistory

import fpbase
from fpbase import _fetch
from fpbase._fetch import FPbaseClient

if TYPE_CHECKING:
    from collections.abc import Callable, Iterator
    from pathlib import Path

RESPONSE = {"data": {"dyes": [], "proteins": []}}


class FakeServer:
    def __init__(self) -> None:
        self.statuses: list[int] = []  # statuses to send before succeeding
        self.response_headers: dict[str, str] = {}
        self.error_body: dict | None = None
        self.request_headers: list[dict[str, str]] = []
        server = self

        class Handler(BaseHTTPRequestHandler):
            def do_POST(self) -> None:
                server.request_headers.append(dict(self.headers))
                self.rfile.read(int(self.headers["Content-Length"]))
                status = server.statuses.pop(0) if server.statuses else 200
                self.send_response(status)
                if status == 429:
                    self.send_header("Retry-After", "0")
                for key, value in server.response_headers.items():
                    self.send_header(key, value)
                self.end_headers()
                if status == 200:
                    self.wfile.write(json.dumps(RESPONSE).encode())
                elif server.error_body is not None:
                    self.wfile.write(json.dumps(server.error_body).encode())

            def log_message(self, *args: object) -> None:
                pass

        self.httpd = ThreadingHTTPServer(("127.0.0.1", 0), Handler)
        self.url = f"http://127.0.0.1:{self.httpd.server_address[1]}/graphql/"


@pytest.fixture
def server() -> Iterator[FakeServer]:
    srv = FakeServer()
    thread = threading.Thread(target=srv.httpd.serve_forever, daemon=True)
    thread.start()
    yield srv
    srv.httpd.shutdown()
    srv.httpd.server_close()


@pytest.fixture
def make_client(server: FakeServer) -> Iterator[Callable[..., FPbaseClient]]:
    clients: list[FPbaseClient] = []

    def _make(**kwargs: Any) -> FPbaseClient:
        clients.append(FPbaseClient(base_url=server.url, **kwargs))
        return clients[-1]

    yield _make
    for client in clients:
        client.session.close()


@pytest.fixture(autouse=True)
def cache_dir(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> Path:
    monkeypatch.setenv("FPBASE_CACHE_DIR", str(tmp_path))
    monkeypatch.delenv("FPBASE_API_KEY", raising=False)
    monkeypatch.setattr(_fetch, "_shown_notices", set())
    # don't wait between retries
    monkeypatch.setattr(_fetch, "_RETRY", _fetch._RETRY.new(backoff_factor=0))
    return tmp_path


def test_retries_throttled_requests(
    server: FakeServer, make_client: Callable[..., FPbaseClient]
) -> None:
    server.statuses = [429, 429]
    client = make_client()
    assert json.loads(client._send_query("{ dyes { id } }")) == RESPONSE
    assert len(server.request_headers) == 3


def test_raises_when_retries_exhausted(
    server: FakeServer, make_client: Callable[..., FPbaseClient]
) -> None:
    server.statuses = [429] * 10
    client = make_client()
    with pytest.raises(requests.HTTPError, match="429"):
        client._send_query("{ dyes { id } }")
    assert len(server.request_headers) == 6  # the first try, and 5 retries


@pytest.mark.parametrize("status", [502, 503, 504])
def test_server_errors_retried_once(
    server: FakeServer, make_client: Callable[..., FPbaseClient], status: int
) -> None:
    server.statuses = [status]
    client = make_client()
    assert json.loads(client._send_query("{ dyes { id } }")) == RESPONSE
    assert len(server.request_headers) == 2


def test_server_errors_not_retried_twice(
    server: FakeServer, make_client: Callable[..., FPbaseClient]
) -> None:
    # e.g. a query that runs past the server's timeout: resending it won't help
    server.statuses = [503, 503, 503]
    client = make_client()
    with pytest.raises(requests.HTTPError, match="503"):
        client._send_query("{ dyes { id } }")
    assert len(server.request_headers) == 2


def test_throttling_then_server_error(
    server: FakeServer, make_client: Callable[..., FPbaseClient]
) -> None:
    server.statuses = [429, 429, 503]
    client = make_client()
    assert json.loads(client._send_query("{ dyes { id } }")) == RESPONSE
    assert len(server.request_headers) == 4


def test_server_error_retry_waits() -> None:
    """The one retry after a server error pauses; throttling uses Retry-After."""
    retry = _fetch._Retry(backoff_factor=1)
    after = {
        status: retry.new(history=(RequestHistory("POST", "/", None, status, None),))
        for status in (429, 503)
    }
    assert after[503].get_backoff_time() == 5
    assert after[429].get_backoff_time() == 0


def test_disk_cache_shared_between_clients(
    server: FakeServer, make_client: Callable[..., FPbaseClient], cache_dir: Path
) -> None:
    query = "{ dyes { id name slug } proteins { id name slug } }"
    make_client()._send_query(query, persist=True)
    assert len(list(cache_dir.glob("*.json"))) == 1
    # a new client (e.g. a new process) reads from disk instead of the server
    assert make_client()._fluorophore_ids == {}
    assert len(server.request_headers) == 1


def test_non_persisted_queries_not_cached_on_disk(
    make_client: Callable[..., FPbaseClient], cache_dir: Path
) -> None:
    make_client()._send_query("{ dyes { id } }")
    assert not list(cache_dir.glob("*.json"))


def test_user_agent_includes_version(
    server: FakeServer, make_client: Callable[..., FPbaseClient]
) -> None:
    make_client()._send_query("{ dyes { id } }")
    assert server.request_headers[0]["User-Agent"].startswith(
        f"fpbase-py/{fpbase.__version__} "
    )


def test_no_api_key_by_default(
    server: FakeServer, make_client: Callable[..., FPbaseClient]
) -> None:
    make_client()._send_query("{ dyes { id } }")
    assert "Authorization" not in server.request_headers[0]


@pytest.mark.parametrize("from_env", [True, False])
def test_api_key(
    server: FakeServer,
    make_client: Callable[..., FPbaseClient],
    monkeypatch: pytest.MonkeyPatch,
    from_env: bool,
) -> None:
    if from_env:
        monkeypatch.setenv("FPBASE_API_KEY", "secret")
        client = make_client()
    else:
        client = make_client(api_key="secret")
    client._send_query("{ dyes { id } }")
    assert server.request_headers[0]["Authorization"] == "Bearer secret"


@pytest.mark.parametrize(
    ("url", "sent"),
    [
        ("https://www.fpbase.org/graphql/", True),
        ("http://localhost:8000/graphql/", True),
        ("http://www.fpbase.org/graphql/", False),
    ],
)
def test_api_key_only_sent_securely(url: str, sent: bool) -> None:
    request = requests.Request("POST", url).prepare()
    _fetch._ApiKeyAuth("secret")(request)
    assert ("Authorization" in request.headers) is sent


def test_server_notice_warns_once(
    server: FakeServer, make_client: Callable[..., FPbaseClient]
) -> None:
    server.response_headers = {_fetch.NOTICE_HEADER: "Get an API key!"}
    client = make_client()
    with pytest.warns(fpbase.FPbaseWarning, match="Get an API key!"):
        client._send_query("{ dyes { id } }")
    client._send_query("{ proteins { id } }")  # not warned again (warnings are errors)


def test_bot_challenge_error(
    server: FakeServer, make_client: Callable[..., FPbaseClient]
) -> None:
    server.statuses = [403]
    server.response_headers = {"cf-mitigated": "challenge"}
    with pytest.raises(requests.HTTPError, match="bot protection"):
        make_client()._send_query("{ dyes { id } }")


@pytest.mark.parametrize(
    "body",
    [{"detail": "Invalid API key."}, {"errors": [{"message": "Invalid API key."}]}],
)
def test_error_includes_server_message(
    server: FakeServer, make_client: Callable[..., FPbaseClient], body: dict
) -> None:
    server.statuses = [401]
    server.error_body = body
    with pytest.raises(requests.HTTPError, match="HTTP 401: Invalid API key"):
        make_client()._send_query("{ dyes { id } }")
