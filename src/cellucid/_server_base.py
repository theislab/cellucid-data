"""
Shared base classes for Cellucid HTTP servers.

This module contains the common functionality used by both:
- CellucidServer (serves pre-exported data)
- AnnDataServer (serves AnnData directly)

It's an internal module (prefixed with _) and not part of the public API.
"""

from __future__ import annotations

import json
import logging
import socket
import threading
from functools import partial
from http import HTTPStatus
from http.server import HTTPServer, SimpleHTTPRequestHandler
from pathlib import Path
from typing import Any, Callable

logger = logging.getLogger("cellucid.server")

DEFAULT_PORT = 8765
DEFAULT_HOST = "127.0.0.1"
CELLUCID_WEB_URL = "https://www.cellucid.com"


# =============================================================================
# CLI Output Helpers
# =============================================================================
# These functions provide consistent, professional CLI output across all servers.

def print_step(step: int, total: int, title: str):
    """Print a step header like [1/4] Title..."""
    print(f"\n[{step}/{total}] {title}...")


def print_detail(label: str, value: str):
    """Print an indented detail line."""
    print(f"      {label}: {value}")


def print_success(message: str = "Done"):
    """Print a success indicator with checkmark."""
    print(f"      ✓ {message}")


def print_server_banner(url: str, viewer_url: str):
    """Print the final server ready banner."""
    print(f"\n{'═' * 60}")
    print("  CELLUCID SERVER RUNNING")
    print(f"{'═' * 60}")
    print()
    print(f"  Local URL:    {url}")
    print(f"  Viewer URL:   {viewer_url}")
    print()
    print("  Press Ctrl+C to stop")
    print()
    print(f"{'═' * 60}")


# Global registry for event callbacks (viewer_id -> callback function)
# This allows HTTP event handlers to route events to the correct viewer
_event_callbacks: dict[str, Callable[[dict], None]] = {}


def register_event_callback(viewer_id: str, callback: Callable[[dict], None]):
    """Register a callback to receive events for a viewer."""
    _event_callbacks[viewer_id] = callback


def unregister_event_callback(viewer_id: str):
    """Unregister a viewer's event callback."""
    _event_callbacks.pop(viewer_id, None)


def route_event(viewer_id: str, event: dict) -> bool:
    """
    Route an event to the appropriate viewer.

    Returns True if the event was delivered, False otherwise.
    """
    callback = _event_callbacks.get(viewer_id)
    if callback:
        try:
            callback(event)
            return True
        except Exception as e:
            logger.error(f"Error in event callback for viewer {viewer_id}: {e}")
    return False


def find_free_port(host: str = DEFAULT_HOST, start_port: int = DEFAULT_PORT, max_attempts: int = 100) -> int:
    """
    Find a free port starting from start_port.

    Parameters
    ----------
    host : str
        Host to bind to.
    start_port : int
        Port to start searching from.
    max_attempts : int
        Maximum number of ports to try.

    Returns
    -------
    int
        A free port number.

    Raises
    ------
    RuntimeError
        If no free port is found within max_attempts.
    """
    for port in range(start_port, start_port + max_attempts):
        try:
            with socket.socket(socket.AF_INET, socket.SOCK_STREAM) as s:
                s.bind((host, port))
                return port
        except OSError:
            continue
    raise RuntimeError(f"Could not find a free port starting from {start_port}")


def ensure_port_available(host: str, port: int, quiet: bool = False) -> int:
    """
    Ensure the specified port is available, finding a new one if needed.

    Parameters
    ----------
    host : str
        Host to bind to.
    port : int
        Preferred port.
    quiet : bool
        If True, suppress log messages.

    Returns
    -------
    int
        The available port (may be different from input if port was in use).
    """
    try:
        with socket.socket(socket.AF_INET, socket.SOCK_STREAM) as s:
            s.bind((host, port))
            return port
    except OSError:
        new_port = find_free_port(host, port + 1)
        if not quiet:
            logger.info(f"Port {port} in use, using {new_port}")
        return new_port


class CORSMixin:
    """
    Mixin class providing CORS headers for HTTP handlers.

    Use in conjunction with SimpleHTTPRequestHandler or similar.
    """

    # Override in subclass: whether to allow caching
    allow_caching: bool = True

    def _get_allowed_origin(self) -> str | None:
        origin = self.headers.get("Origin")
        if not origin:
            return None

        from urllib.parse import urlparse

        try:
            parsed = urlparse(origin)
        except Exception:
            return None

        if parsed.scheme not in ("http", "https"):
            return None

        host = parsed.hostname
        if not host:
            return None

        if host in ("localhost", "127.0.0.1"):
            return origin

        # Allow the canonical hosted viewer origin.
        if origin == CELLUCID_WEB_URL:
            return origin

        # Allow https subdomains under cellucid.com (optional, but safe).
        if parsed.scheme == "https" and (host == "cellucid.com" or host.endswith(".cellucid.com")):
            return origin

        return None

    def add_cors_headers(self):
        """Add CORS headers to the response."""
        allowed_origin = self._get_allowed_origin()
        if allowed_origin:
            self.send_header("Access-Control-Allow-Origin", allowed_origin)
            self.send_header("Vary", "Origin")
        self.send_header("Access-Control-Allow-Methods", "GET, HEAD, POST, OPTIONS")
        self.send_header("Access-Control-Allow-Headers", "Content-Type, Range")
        self.send_header("Access-Control-Expose-Headers", "Content-Length, Content-Range")
        if self.allow_caching:
            self.send_header("Cache-Control", "public, max-age=3600")
        else:
            self.send_header("Cache-Control", "no-cache")

    def do_OPTIONS(self):
        """Handle CORS preflight requests."""
        self.send_response(HTTPStatus.NO_CONTENT)
        # Note: CORS headers are added by end_headers() override in subclasses
        self.end_headers()

    def send_json(self, data: dict, head_only: bool = False):
        """Send a JSON response with CORS headers."""
        body = json.dumps(data).encode("utf-8")
        self.send_response(HTTPStatus.OK)
        self.send_header("Content-Type", "application/json")
        self.send_header("Content-Length", len(body))
        # Note: CORS headers are added by end_headers() override in subclasses
        self.end_headers()
        if not head_only:
            self.wfile.write(body)

    def send_binary(
        self,
        data: bytes,
        content_type: str = "application/octet-stream",
        head_only: bool = False,
        compressed: bool = False,
    ):
        """Send binary data with CORS headers."""
        self.send_response(HTTPStatus.OK)
        self.send_header("Content-Type", content_type)
        self.send_header("Content-Length", len(data))
        if compressed:
            self.send_header("Content-Encoding", "gzip")
        # Note: CORS headers are added by end_headers() override in subclasses
        self.end_headers()
        if not head_only:
            self.wfile.write(data)

    def send_error_response(self, code: int, message: str):
        """Send an error response with CORS headers."""
        self.send_response(code)
        self.send_header("Content-Type", "text/plain")
        body = message.encode("utf-8")
        self.send_header("Content-Length", len(body))
        # Note: CORS headers are added by end_headers() override in subclasses
        self.end_headers()
        if getattr(self, "command", "") != "HEAD":
            self.wfile.write(body)

    def handle_event_post(self) -> bool:
        """
        Handle POST request to /_cellucid/events endpoint.

        This enables Frontend → Python communication for hooks in ALL environments.
        The frontend POSTs events here, and they get routed to the appropriate viewer.

        Returns True if the request was handled, False if path doesn't match.
        """
        from urllib.parse import urlparse, unquote

        parsed = urlparse(self.path)
        path = unquote(parsed.path)

        if path != "/_cellucid/events":
            return False

        try:
            # Read the POST body
            content_length = int(self.headers.get("Content-Length", 0))
            if content_length > 1024 * 1024:  # 1MB limit
                self.send_error_response(413, "Request too large")
                return True

            body = self.rfile.read(content_length)
            event = json.loads(body.decode("utf-8"))

            # Route to the appropriate viewer
            viewer_id = event.get("viewerId")
            if not viewer_id:
                self.send_json({"status": "error", "message": "Missing viewerId"})
                return True

            delivered = route_event(viewer_id, event)
            self.send_json({
                "status": "ok" if delivered else "not_found",
                "delivered": delivered
            })

        except json.JSONDecodeError as e:
            self.send_json({"status": "error", "message": f"Invalid JSON: {e}"})
        except Exception:
            logger.exception("Error handling event POST")
            self.send_error_response(500, "Internal server error")

        return True
