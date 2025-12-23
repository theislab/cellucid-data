"""
AnnData Server for Cellucid

HTTP server that serves AnnData data dynamically in the Cellucid format.
This allows direct visualization of AnnData files without pre-export.

Supports:
- h5ad files: Single HDF5 files (backed mode for lazy loading)
- zarr stores: Directory-based stores (inherently lazy-loaded)
- In-memory AnnData objects

Usage:
    from cellucid.anndata_server import serve_anndata

    # Serve an h5ad file
    serve_anndata("/path/to/data.h5ad")

    # Serve a zarr store (must be a directory)
    serve_anndata("/path/to/data.zarr")

    # Or an in-memory AnnData
    serve_anndata(adata)
"""

from __future__ import annotations

import logging
import re
import socket
import threading
import webbrowser
from functools import partial
from http import HTTPStatus
from http.server import HTTPServer, SimpleHTTPRequestHandler
from pathlib import Path
from typing import TYPE_CHECKING, Optional, Union
from urllib.parse import parse_qs, unquote, urlparse

from .anndata_adapter import AnnDataAdapter, _safe_filename_component
from ._server_base import (
    CORSMixin,
    DEFAULT_PORT,
    DEFAULT_HOST,
    CELLUCID_WEB_URL,
    ensure_port_available,
    print_step,
    print_detail,
    print_success,
    print_server_banner,
)

if TYPE_CHECKING:
    import anndata

logger = logging.getLogger("cellucid.anndata_server")


class AnnDataRequestHandler(CORSMixin, SimpleHTTPRequestHandler):
    """
    HTTP handler for serving AnnData data in Cellucid format.

    Routes:
        /dataset_identity.json - Dataset metadata
        /obs_manifest.json - Cell metadata manifest
        /var_manifest.json - Gene expression manifest
        /connectivity_manifest.json - Connectivity manifest
        /points_{dim}d.bin - Embedding coordinates
        /vectors/{fieldId}_{dim}d.bin - Vector field displacement vectors
        /obs/{field}.values.f32.bin - Continuous obs field
        /obs/{field}.codes.{ext}.bin - Categorical obs field codes
        /obs/{field}.outliers.f32.bin - Categorical outlier quantiles
        /var/{gene}.values.f32.bin - Gene expression values
        /connectivity/edges.src.bin - Edge sources
        /connectivity/edges.dst.bin - Edge destinations
        /_cellucid/health - Health check
        /_cellucid/info - Server info
        /_cellucid/events - POST endpoint for frontend events (hooks)
    """

    allow_caching = False  # Data is dynamic, don't cache

    def __init__(self, *args, adapter: AnnDataAdapter, server_info: dict, dataset_id: str = "", **kwargs):
        self.adapter = adapter
        self.server_info = server_info
        self.dataset_id = dataset_id  # Cached for path prefix stripping
        # Don't call super().__init__ with directory since we're serving virtual files
        super(SimpleHTTPRequestHandler, self).__init__(*args, **kwargs)

    def end_headers(self):
        """Add CORS headers to every response."""
        self.add_cors_headers()
        super().end_headers()

    def do_POST(self):
        """Handle POST requests (events from frontend)."""
        if self.handle_event_post():
            return
        # No other POST endpoints - return 404
        self.send_error_response(404, f"POST not supported for path: {self.path}")

    def do_HEAD(self):
        """Handle HEAD requests (for urlExists checks)."""
        self.do_GET(head_only=True)

    def do_GET(self, head_only: bool = False):
        """Handle GET requests."""
        parsed = urlparse(self.path)
        path = unquote(parsed.path).lstrip("/")

        # Strip dataset ID prefix if present (frontend may include it)
        # e.g., "my_dataset/points_3d.bin" -> "points_3d.bin"
        if self.dataset_id and path.startswith(f"{self.dataset_id}/"):
            path = path[len(self.dataset_id) + 1:]

        # Check for gzip support
        accept_encoding = self.headers.get("Accept-Encoding", "")
        supports_gzip = "gzip" in accept_encoding

        try:
            # Root path - redirect to viewer
            if path == "" or path == "index.html":
                viewer_url = f"{CELLUCID_WEB_URL}?remote=http://{self.server_info['host']}:{self.server_info['port']}&anndata=true"
                self.send_response(302)
                self.send_header("Location", viewer_url)
                self.end_headers()
                return

            if path == "_cellucid/health":
                self.send_json({
                    "status": "ok",
                    "type": "anndata",
                    "version": self.server_info.get("version", "unknown"),
                    "format": self.server_info.get("format", "unknown"),
                    "is_backed": self.server_info.get("is_backed", False),
                    "n_cells": self.adapter.n_cells,
                    "n_genes": self.adapter.n_genes,
                }, head_only)
            elif path == "_cellucid/info":
                self.send_json(self.server_info, head_only)
            elif path == "_cellucid/datasets":
                # Return a single dataset entry for the loaded AnnData
                identity = self.adapter.get_dataset_identity()
                datasets = [{
                    "id": identity.get("id", "anndata"),
                    "path": "/",
                    "name": identity.get("name", "AnnData"),
                }]
                self.send_json({"datasets": datasets}, head_only)
            elif path == "dataset_identity.json":
                self.send_json(self.adapter.get_dataset_identity(), head_only)
            elif path == "obs_manifest.json":
                self.send_json(self.adapter.get_obs_manifest(), head_only)
            elif path == "var_manifest.json":
                self.send_json(self.adapter.get_var_manifest(), head_only)
            elif path == "connectivity_manifest.json":
                manifest = self.adapter.get_connectivity_manifest()
                if manifest:
                    self.send_json(manifest, head_only)
                else:
                    self.send_error_response(404, "No connectivity data")
            elif path.startswith("points_") and path.endswith(".bin"):
                self._handle_points(path, head_only, supports_gzip)
            elif path.startswith("obs/"):
                self._handle_obs(path, head_only, supports_gzip)
            elif path.startswith("var/"):
                self._handle_var(path, head_only, supports_gzip)
            elif path.startswith("connectivity/"):
                self._handle_connectivity(path, head_only, supports_gzip)
            elif path.startswith("vectors/"):
                self._handle_vector_fields(path, head_only, supports_gzip)
            else:
                self.send_error_response(404, f"Not found: {path}")

        except Exception:
            logger.exception("Error handling %s", path)
            self.send_error_response(500, "Internal server error")

    def _handle_points(self, path: str, head_only: bool, supports_gzip: bool):
        """Handle points_Xd.bin requests."""
        # Parse dimension from filename: points_1d.bin, points_2d.bin, points_3d.bin
        # Also handle .gz suffix for explicit compressed requests
        clean_path = path[:-3] if path.endswith(".gz") else path
        match = re.match(r"points_(\d)d\.bin", clean_path)
        if not match:
            self.send_error_response(404, f"Invalid points path: {path}")
            return

        dim = int(match.group(1))
        try:
            # Compress if client supports gzip (transparent compression for better network perf)
            compress = supports_gzip
            data = self.adapter.get_points_binary(dim, compress=compress)
            self.send_binary(data, head_only=head_only, compressed=compress)
        except ValueError:
            self.send_error_response(404, "Points data not available")

    def _handle_vector_fields(self, path: str, head_only: bool, supports_gzip: bool):
        """
        Handle vector field requests under /vectors/.

        Supported paths:
        - vectors/{fieldId}_{dim}d.bin
        - vectors/{fieldId}_{dim}d.bin.gz  (raw gzip bytes, no Content-Encoding)
        """
        raw_gz = path.endswith(".gz")
        clean_path = path[:-3] if raw_gz else path
        match = re.match(r"vectors/(.+)_([123])d\.bin$", clean_path)
        if not match:
            self.send_error_response(404, f"Invalid vectors path: {path}")
            return

        field_id = match.group(1)
        dim = int(match.group(2))

        try:
            if raw_gz:
                data = self.adapter.get_vector_field_binary(field_id, dim, compress=True)
                self.send_binary(data, head_only=head_only, compressed=False)
            else:
                compress = supports_gzip
                data = self.adapter.get_vector_field_binary(field_id, dim, compress=compress)
                self.send_binary(data, head_only=head_only, compressed=compress)
        except ValueError as e:
            self.send_error_response(404, str(e))

    def _handle_obs(self, path: str, head_only: bool, supports_gzip: bool):
        """
        Handle obs field requests.

        Supported paths:
        - obs/{field}.values.f32
        - obs/{field}.codes.u8 or obs/{field}.codes.u16
        - obs/{field}.outliers.f32
        """
        # Remove 'obs/' prefix
        filename = path[4:]

        # Strip suffixes for parsing
        if filename.endswith(".gz"):
            filename = filename[:-3]

        # Compress if client supports gzip (transparent compression)
        compress = supports_gzip

        # Parse filename: {field}.{type}.{dtype}
        # Examples: cell_type.codes.u8, n_counts.values.f32, cluster.outliers.f32
        parts = filename.rsplit(".", 2)  # Split from right: [field, type, dtype]
        if len(parts) < 3:
            self.send_error_response(404, f"Invalid obs path format: {path}")
            return

        field_name = parts[0]
        data_type = parts[1]  # 'values', 'codes', or 'outliers'

        # Find the actual field key (may differ from safe filename)
        obs_keys = self.adapter.get_obs_keys()
        actual_key = None
        for key in obs_keys:
            if _safe_filename_component(key) == field_name:
                actual_key = key
                break

        if actual_key is None:
            self.send_error_response(404, f"Obs field not found: {field_name}")
            return

        try:
            if data_type == "values":
                # Continuous field values
                data = self.adapter.get_obs_continuous_values(actual_key, compress=compress)
                self.send_binary(data, head_only=head_only, compressed=compress)
            elif data_type == "codes":
                # Categorical codes
                data, categories, missing = self.adapter.get_obs_categorical_codes(actual_key, compress=compress)
                self.send_binary(data, head_only=head_only, compressed=compress)
            elif data_type == "outliers":
                # Outlier quantiles for categorical field
                data = self.adapter.get_obs_outlier_quantiles(actual_key, compress=compress)
                self.send_binary(data, head_only=head_only, compressed=compress)
            else:
                self.send_error_response(404, f"Unknown obs data type: {data_type} (expected values, codes, or outliers)")
        except Exception:
            logger.exception("Error handling obs request: %s", path)
            self.send_error_response(500, "Internal server error")

    def _handle_var(self, path: str, head_only: bool, supports_gzip: bool):
        """
        Handle var (gene expression) requests.

        Supported paths:
        - var/{gene}.values.f32
        """
        # Remove 'var/' prefix
        filename = path[4:]

        # Strip suffixes for parsing
        if filename.endswith(".gz"):
            filename = filename[:-3]

        # Compress if client supports gzip (transparent compression)
        compress = supports_gzip

        # Check format: {gene}.values.f32
        if not filename.endswith(".values.f32"):
            self.send_error_response(404, f"Invalid var path format: {path} (expected {{gene}}.values.f32)")
            return

        gene_safe = filename[:-11]  # Remove '.values.f32'

        # Find actual gene ID (may differ from safe filename due to special characters)
        gene_ids = self.adapter.get_gene_ids()
        actual_gene = None

        # First try safe filename match
        for gid in gene_ids:
            if _safe_filename_component(gid) == gene_safe:
                actual_gene = gid
                break

        # Fallback: try exact match (in case gene ID is already safe)
        if actual_gene is None and gene_safe in gene_ids:
            actual_gene = gene_safe

        if actual_gene is None:
            self.send_error_response(404, f"Gene not found: {gene_safe}")
            return

        try:
            data = self.adapter.get_gene_expression(actual_gene, compress=compress)
            self.send_binary(data, head_only=head_only, compressed=compress)
        except KeyError:
            self.send_error_response(404, "Gene not found")
        except Exception:
            logger.exception("Error handling var request: %s", path)
            self.send_error_response(500, "Internal server error")

    def _handle_connectivity(self, path: str, head_only: bool, supports_gzip: bool):
        """
        Handle connectivity edge requests.

        Supported paths:
        - connectivity/edges.src.bin - Edge source indices
        - connectivity/edges.dst.bin - Edge destination indices
        """
        filename = path.split("/")[-1]

        # Strip suffixes for parsing
        if filename.endswith(".gz"):
            filename = filename[:-3]

        # Compress if client supports gzip (transparent compression)
        compress = supports_gzip

        try:
            sources_data, dests_data, n_edges, max_neighbors = self.adapter.get_connectivity_edges(compress=compress)

            if filename == "edges.src.bin":
                self.send_binary(sources_data, head_only=head_only, compressed=compress)
            elif filename == "edges.dst.bin":
                self.send_binary(dests_data, head_only=head_only, compressed=compress)
            else:
                self.send_error_response(404, f"Unknown connectivity file: {filename} (expected edges.src.bin or edges.dst.bin)")
        except ValueError:
            self.send_error_response(404, "Connectivity data not available")
        except Exception:
            logger.exception("Error handling connectivity request: %s", path)
            self.send_error_response(500, "Internal server error")

    def log_message(self, format: str, *args):
        """Override to use Python logging."""
        logger.debug("%s - %s", self.address_string(), format % args)


class AnnDataServer:
    """
    Server for serving AnnData data in Cellucid format.

    Example:
        server = AnnDataServer(adata)
        server.start()  # Blocking

        # Or non-blocking:
        server.start_background()
        # ... do other things ...
        server.stop()
    """

    def __init__(
        self,
        data: Union[str, Path, "anndata.AnnData"],
        port: int = DEFAULT_PORT,
        host: str = DEFAULT_HOST,
        open_browser: bool = False,
        quiet: bool = False,
        **adapter_kwargs,
    ):
        """
        Initialize the server.

        Parameters
        ----------
        data : str, Path, or AnnData
            Path to h5ad file or AnnData object.
        port : int
            Port to serve on.
        host : str
            Host to bind to.
        open_browser : bool
            Whether to open browser on start.
        quiet : bool
            Suppress info messages.
        **adapter_kwargs
            Additional arguments passed to AnnDataAdapter.
        """
        self.port = port
        self.host = host
        self.open_browser = open_browser
        self.quiet = quiet

        # Step 1: Detect format
        if isinstance(data, (str, Path)):
            data_path = Path(data)
            is_zarr = str(data).endswith('.zarr') or data_path.is_dir()
            format_name = "zarr" if is_zarr else "h5ad"
            self.data_source = str(data)
            self.data_format = format_name

            if not quiet:
                print_step(1, 4, "Detecting format")
                print_detail("Path", str(data))
                print_detail("Format", format_name)
                print_success("Format detected")
        else:
            # In-memory AnnData
            self.data_source = "in-memory AnnData"
            self.data_format = "in-memory"
            if not quiet:
                print_step(1, 4, "Detecting format")
                print_detail("Source", "in-memory AnnData")
                print_success("Format detected")

        # Step 2: Load file
        if not quiet:
            print_step(2, 4, "Loading AnnData")
            backed = adapter_kwargs.get("backed", True)
            mode = "backed (lazy loading)" if backed else "in-memory"
            print_detail("Mode", mode)

        if isinstance(data, (str, Path)):
            self.adapter = AnnDataAdapter.from_file(data, **adapter_kwargs)
        else:
            self.adapter = AnnDataAdapter(data, **adapter_kwargs)

        if not quiet:
            print_success("File opened")

        # Step 3: Analyze dataset
        if not quiet:
            print_step(3, 4, "Analyzing dataset")
            self._print_dataset_info()
            print_success("Analysis complete")

        self._server: Optional[HTTPServer] = None
        self._thread: Optional[threading.Thread] = None
        self._running = False

        # Version
        try:
            from . import __version__
        except ImportError:
            __version__ = "0.0.0"

        self.server_info = {
            "version": __version__,
            "type": "anndata",
            "format": self.data_format,
            "host": self.host,
            "port": self.port,
            "n_cells": self.adapter.n_cells,
            "n_genes": self.adapter.n_genes,
            "is_backed": self.adapter.is_backed,
        }

    def _print_dataset_info(self):
        """Print information about the loaded dataset."""
        print_detail("Cells", f"{self.adapter.n_cells:,}")
        print_detail("Genes", f"{self.adapter.n_genes:,}")

        # Get embedding info
        dims = self.adapter._available_dimensions
        dim_str = ", ".join(f"{d}D" for d in sorted(dims)) if dims else "none"
        print_detail("Embeddings", dim_str)

        # Count obs field types
        obs_keys = self.adapter.get_obs_keys()
        n_categorical = sum(1 for k in obs_keys if self.adapter.get_obs_field_kind(k) == "category")
        n_continuous = sum(1 for k in obs_keys if self.adapter.get_obs_field_kind(k) == "continuous")
        print_detail("Obs fields", f"{n_categorical} categorical, {n_continuous} continuous")

        # Check connectivity
        has_conn = self.adapter.has_connectivity()
        if has_conn:
            try:
                _, _, n_edges, _ = self.adapter.get_connectivity_edges(compress=False)
                conn_str = f"yes ({n_edges:,} edges)"
            except Exception:
                conn_str = "yes"
        else:
            conn_str = "no"
        print_detail("Connectivity", conn_str)

    @property
    def url(self) -> str:
        """Get the server URL."""
        return f"http://{self.host}:{self.port}"

    @property
    def viewer_url(self) -> str:
        """Get the full URL to open the viewer."""
        return f"{CELLUCID_WEB_URL}?remote={self.url}&anndata=true"

    def start(self, blocking: bool = True):
        """Start the server."""
        if self._running:
            logger.warning("Server is already running")
            return

        # Step 4: Start server
        if not self.quiet:
            print_step(4, 4, "Starting server")

        # Ensure port is available (finds new one if needed)
        self.port = ensure_port_available(self.host, self.port, self.quiet)
        self.server_info["port"] = self.port

        # Create handler with adapter
        # Cache dataset_id for efficient path prefix stripping
        dataset_id = self.adapter.get_dataset_identity().get("id", "")
        handler = partial(
            AnnDataRequestHandler,
            adapter=self.adapter,
            server_info=self.server_info,
            dataset_id=dataset_id,
        )

        self._server = HTTPServer((self.host, self.port), handler)
        self._running = True

        if not self.quiet:
            print_success("Server ready")
            print_server_banner(self.url, self.viewer_url)

        if self.open_browser:
            webbrowser.open(self.viewer_url)

        if blocking:
            try:
                self._server.serve_forever()
            except KeyboardInterrupt:
                if not self.quiet:
                    print("\nShutting down server...")
                # Don't call shutdown() here - it would deadlock since we're in the same thread
                # Just close the server and cleanup
                self._running = False
                if self._server:
                    self._server.server_close()
                    self._server = None
                if hasattr(self, 'adapter') and self.adapter:
                    self.adapter.close()
                if not self.quiet:
                    print("AnnData server stopped")
        else:
            self._thread = threading.Thread(target=self._server.serve_forever, daemon=True)
            self._thread.start()

    def start_background(self):
        """Start the server in a background thread."""
        self.start(blocking=False)

    def stop(self):
        """Stop the server and cleanup resources."""
        self._running = False

        if self._server:
            # shutdown() tells serve_forever() to stop, but we also need
            # server_close() to release the socket immediately
            self._server.shutdown()
            self._server.server_close()
            self._server = None

        # IMPORTANT: Close the adapter to release memory and file handles
        if hasattr(self, 'adapter') and self.adapter:
            self.adapter.close()

        if not self.quiet:
            print("AnnData server stopped")

    def is_running(self) -> bool:
        """Check if the server is running."""
        return self._running

    def wait(self):
        """Wait for the server to stop."""
        if self._thread:
            try:
                while self._running:
                    self._thread.join(timeout=1)
            except KeyboardInterrupt:
                self.stop()


def serve_anndata(
    data: Union[str, Path, "anndata.AnnData"],
    port: int = DEFAULT_PORT,
    host: str = DEFAULT_HOST,
    open_browser: bool = True,
    quiet: bool = False,
    **kwargs,
) -> AnnDataServer:
    """
    Serve an AnnData object or h5ad file directly.

    This is a convenience function for quickly visualizing AnnData.
    For production use, consider using prepare instead.

    Parameters
    ----------
    data : str, Path, or AnnData
        Path to h5ad file or AnnData object.
    port : int
        Port to serve on (default: 8765).
    host : str
        Host to bind to (default: 127.0.0.1).
    open_browser : bool
        Whether to open browser (default: True).
    quiet : bool
        Suppress info messages.
    **kwargs
        Additional arguments passed to AnnDataAdapter.

    Returns
    -------
    AnnDataServer
        The running server instance.

    Example
    -------
    >>> from cellucid import serve_anndata
    >>> serve_anndata("/path/to/data.h5ad")

    >>> # Or with in-memory AnnData
    >>> import anndata as ad
    >>> adata = ad.read_h5ad("data.h5ad")
    >>> serve_anndata(adata)
    """
    server = AnnDataServer(
        data=data,
        port=port,
        host=host,
        open_browser=open_browser,
        quiet=quiet,
        **kwargs,
    )
    server.start()
    return server
