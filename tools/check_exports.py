"""Inspect a generated exact report and the built wheel/source distribution."""

import argparse
import csv
import json
import re
import tarfile
import zipfile
from html.parser import HTMLParser
from pathlib import Path
from urllib.parse import unquote, urlsplit
from xml.etree import ElementTree

ROOT = Path(__file__).resolve().parents[1]


class Targets(HTMLParser):
    def __init__(self):
        super().__init__()
        self.ids = []
        self.open_targets = []
        self.urls = []

    def handle_starttag(self, tag, attrs):
        values = dict(attrs)
        if "id" in values:
            self.ids.append(values["id"])
        if "data-open" in values:
            self.open_targets.append(values["data-open"])
        for key in ("href", "src"):
            if values.get(key):
                self.urls.append(values[key])


def check(report, dist):
    summary = json.loads((report / "verification.json").read_text(encoding="utf-8"))
    assert summary["all_equal"] and summary["checked_cells"] == 77
    assert summary["numeric_tolerance_used"] is False
    pages = []
    for language in ("cn", "en"):
        path = report / f"index-{language}.html"
        parser = Targets()
        parser.feed(path.read_text(encoding="utf-8"))
        assert len(parser.ids) == len(set(parser.ids)), path.name
        cells = [item for item in parser.open_targets if item.startswith("cell-")]
        assert len(cells) == len(set(cells)) == 77
        assert set(cells) <= set(parser.ids)
        local_links = 0
        for url in parser.urls:
            parts = urlsplit(url)
            if parts.scheme or parts.netloc:
                continue
            if not parts.path:
                if parts.fragment:
                    assert unquote(parts.fragment) in parser.ids, url
                continue
            assert (path.parent / unquote(parts.path)).is_file(), (path.name, url)
            local_links += 1
        pages.append({"language": language, "cell_targets": len(cells), "local_links": local_links})
    figures = list((report / "figures").glob("*.svg"))
    for path in figures:
        assert ElementTree.parse(path).getroot().tag.endswith("svg"), path.name
    tex_inputs = 0
    for path in (report / "latex").rglob("*.tex"):
        for name in re.findall(r"\\input\{([^}]+)\}", path.read_text(encoding="utf-8")):
            target = path.parent / name
            if not target.suffix:
                target = target.with_suffix(".tex")
            assert target.is_file(), (path.name, name)
            tex_inputs += 1
    compared = []
    for name in ("path-a.json", "path-b.json", "verification.json"):
        generated = json.loads((report / name).read_text(encoding="utf-8"))
        committed = json.loads((ROOT / "examples/counterterms" / name).read_text(encoding="utf-8"))
        assert generated == committed, name
        compared.append(name)
    with (report / "table.csv").open(encoding="utf-8", newline="") as generated_stream:
        with (ROOT / "examples/counterterms/table.csv").open(
            encoding="utf-8", newline=""
        ) as committed_stream:
            assert list(csv.reader(generated_stream)) == list(csv.reader(committed_stream))
    wheel = dist / "phi3_renormalization-0.2.0-py3-none-any.whl"
    with zipfile.ZipFile(wheel) as archive:
        names = archive.namelist()
        core = [name for name in names if name.startswith("phi3lab/") and name.endswith(".py")]
        assert len(core) == 55, len(core)
        entry_points = archive.read(
            "phi3_renormalization-0.2.0.dist-info/entry_points.txt"
        ).decode()
        for command in ("phi3-ct", "phi3-demo", "phi3lab"):
            assert command + " = " in entry_points
    sdist = dist / "phi3_renormalization-0.2.0.tar.gz"
    with tarfile.open(sdist) as archive:
        names = {name.partition("/")[2] for name in archive.getnames()}
        manifest = json.loads((ROOT / "docs/source-import.json").read_text(encoding="utf-8"))
        for row in manifest["files"]:
            assert row["path"] in names, row["path"]
        for name in ("README_ZH.md", "tests/conftest.py", "tools/check_public_snapshot.py"):
            assert name in names, name
    return {
        "all_passed": True,
        "pages": pages,
        "svg_files_parsed": len(figures),
        "generated_tex_inputs_checked": tex_inputs,
        "saved_example_json_equal": compared,
        "saved_table_equal": True,
        "saved_table_comparison": "Exact CSV fields; platform line endings ignored.",
        "wheel_core_modules": len(core),
        "sdist_imported_files": len(manifest["files"]),
        "tex_compilation_performed": False,
        "visual_review_performed": False,
        "mac_compute_resumed": False,
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--dist", type=Path, default=ROOT / "dist")
    parser.add_argument("--evidence", type=Path)
    args = parser.parse_args()
    rendered = json.dumps(check(args.input, args.dist), indent=2) + "\n"
    if args.evidence:
        args.evidence.write_text(rendered, encoding="utf-8", newline="\n")
    print(rendered, end="")


if __name__ == "__main__":
    main()
