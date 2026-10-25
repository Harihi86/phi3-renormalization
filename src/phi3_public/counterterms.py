"""Portable entry points for the imported exact counterterm engine.

Generation always uses fresh directories. Verification consumes the two saved
paths, never their stored pass/fail result. No command connects to a remote host.
"""

import argparse
import gzip
import json
import socket
import sys
import time
from functools import partial
from http.server import SimpleHTTPRequestHandler, ThreadingHTTPServer
from pathlib import Path


def read_json(path):
    """Read a plain or losslessly compressed record."""
    path = Path(path)
    if path.suffix == ".gz":
        with gzip.open(path, "rt", encoding="utf-8") as stream:
            return json.load(stream)
    if path.is_file():
        return json.loads(path.read_text(encoding="utf-8"))
    with gzip.open(str(path) + ".gz", "rt", encoding="utf-8") as stream:
        return json.load(stream)


def save(path, value):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")


def build(output):
    """Generate the complete two-loop three-point closure and bilingual report."""
    from phi3lab.ct_closure import build as build_closure

    output = Path(output)
    output.mkdir(parents=True, exist_ok=False)
    return build_closure(output)["verification"]


def family(output, loops, legs):
    """Generate both complete paths and their own lower-loop source libraries."""
    from phi3lab.ct_series_a import build_family
    from phi3lab.ct_series_b import ActionLibrary
    from phi3lab.ct_series_report import write_reports
    from phi3lab.ct_series_verify import compare, validate_source_library
    from phi3lab.ct_typed import enumerate_functional

    if type(loops) is not int or loops not in range(1, 5) or legs not in (2, 3):
        raise ValueError("Supported exact families have 1..4 loops and 2 or 3 external legs")
    output = Path(output)
    output.mkdir(parents=True, exist_ok=False)
    start = time.monotonic()
    state = {"state": "running", "loops": loops, "legs": legs, "large_integrals_run": False}
    save(output / "status.json", state)
    try:
        a = build_family(loops, legs)
        save(output / "path-a.json", a)
        library = ActionLibrary()
        b = library.family(loops, legs)
        save(output / "path-b.json", b)
        for (lower_loops, lower_legs), source in library.families.items():
            if (lower_loops, lower_legs) != (loops, legs):
                save(output / "library" / f"L{lower_loops}E{lower_legs}.json", source)
        for source in library.families.values():
            validate_source_library(source, library.families)
        result = compare(a, b)
        save(output / "verification.json", result)
        audit = enumerate_functional(("phi3",) * (2 * loops + legs - 2), tuple(range(1, legs + 1)))
        save(output / "functional-audit.json", audit)
        write_reports(a, b, result, output)
        state.update(state="terminal", exit_code=0 if result["all_equal"] else 1)
        return result
    except BaseException as exc:
        state.update(state="terminal", exit_code=1, error=repr(exc))
        raise
    finally:
        state["seconds"] = time.monotonic() - start
        save(output / "status.json", state)


def verify(directory):
    """Recompute every exact cell and validate B's saved lower-loop sources."""
    directory = Path(directory)
    a = read_json(directory / "path-a.json")
    b = read_json(directory / "path-b.json")
    if "overall_columns" in a:
        from phi3lab.ct_series_verify import compare, validate_source_library

        sources = {(b["loops"], b["legs"]): b}
        if (directory / "library").is_dir():
            for path in sorted((directory / "library").iterdir()):
                if path.name.endswith((".json", ".json.gz")):
                    source = read_json(path.with_suffix("") if path.suffix == ".gz" else path)
                    sources[source["loops"], source["legs"]] = source
        for source in sources.values():
            validate_source_library(source, sources)
        return compare(a, b)
    from phi3lab.ct_verify import compare_paths

    return compare_paths(a, b)


def one_loop(output):
    """Finite one-loop examples; quadrature refinement is not a rigorous bound."""
    from phi3lab.finite import bubble_bph, triangle_bph

    output = Path(output)
    output.mkdir(parents=True, exist_ok=False)
    rows = []
    for name, evaluate in (
        ("bubble", lambda order: bubble_bph(1, 1, order=order)),
        ("triangle", lambda order: triangle_bph((1, 1, 1), 1, order=order)),
    ):
        coarse, fine = evaluate(32), evaluate(64)
        rows.append(
            {
                "graph": name,
                "order32": coarse,
                "order64": fine,
                "refinement_difference": abs(fine - coarse),
            }
        )
    data = {
        "m_squared": 1,
        "external_invariants": 1,
        "loop_measure": "d^6q/pi^3",
        "coupling_and_symmetry_weight_included": False,
        "rows": rows,
        "rigorous_quadrature_error_bound": False,
        "high_loop_acceptance_added": 0,
    }
    save(output / "one-loop.json", data)
    return data


class ReportHandler(SimpleHTTPRequestHandler):
    """Local static viewer with on-demand exact saved-path verification."""

    def do_GET(self):
        if self.path == "/api/verify":
            try:
                data = verify(self.directory)
                code = 200 if data["all_equal"] else 409
            except (ValueError, KeyError, ArithmeticError, OSError) as exc:
                code, data = 409, {"all_equal": False, "error": str(exc)}
            payload = json.dumps(data).encode("utf-8")
            self.send_response(code)
            self.send_header("Content-Type", "application/json; charset=utf-8")
            self.send_header("Content-Length", str(len(payload)))
            self.end_headers()
            self.wfile.write(payload)
            return
        return super().do_GET()


class ReportServer(ThreadingHTTPServer):
    """Use exclusive binding on Windows, matching the retained workbench."""

    allow_reuse_address = False

    def server_bind(self):
        if hasattr(socket, "SO_EXCLUSIVEADDRUSE"):
            self.socket.setsockopt(socket.SOL_SOCKET, socket.SO_EXCLUSIVEADDRUSE, 1)
        super().server_bind()


def serve(directory, language, port):
    directory = Path(directory).resolve()
    page = directory / f"index-{'cn' if language == 'zh' else 'en'}.html"
    if not page.is_file():
        raise FileNotFoundError("Build a report before starting the viewer")
    with ReportServer(
        ("127.0.0.1", port), partial(ReportHandler, directory=str(directory))
    ) as server:
        print(f"http://127.0.0.1:{server.server_port}/{page.name}", flush=True)
        try:
            server.serve_forever()
        except KeyboardInterrupt:
            pass


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest="command", required=True)
    for name, help_text in (
        ("build", "Complete two-loop three-point table, SVG and LaTeX"),
        ("family", "Complete exact 1..4-loop two-/three-point family"),
        ("one-loop", "Finite one-loop numerical examples"),
    ):
        sub = commands.add_parser(name, help=help_text)
        sub.add_argument("--output", required=True, type=Path)
        if name == "family":
            sub.add_argument("--loops", type=int, choices=range(1, 5), required=True)
            sub.add_argument("--legs", type=int, choices=(2, 3), required=True)
    sub = commands.add_parser("verify", help="Recompute exact comparison from both saved paths")
    sub.add_argument("--input", required=True, type=Path)
    sub = commands.add_parser("serve", help="Serve an existing report on localhost")
    sub.add_argument("--input", required=True, type=Path)
    sub.add_argument("--language", choices=("zh", "en"), default="zh")
    sub.add_argument("--port", type=int, default=0)
    args = parser.parse_args(argv)
    try:
        if args.command == "serve":
            serve(args.input, args.language, args.port)
            return 0
        if args.command == "build":
            result = build(args.output)
        elif args.command == "family":
            result = family(args.output, args.loops, args.legs)
        elif args.command == "verify":
            result = verify(args.input)
        else:
            result = one_loop(args.output)
        print(
            json.dumps(
                {k: v for k, v in result.items() if k not in ("cells", "differences")},
                ensure_ascii=False,
                indent=2,
            )
        )
        return 0 if result.get("all_equal", True) else 1
    except (ValueError, KeyError, ArithmeticError, OSError) as exc:
        print(f"{type(exc).__name__}: {exc}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
