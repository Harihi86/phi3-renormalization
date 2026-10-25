"""Check source identity, documentation targets and public evidence boundaries."""

import ast
import hashlib
import json
import re
from pathlib import Path
from urllib.parse import unquote

ROOT = Path(__file__).resolve().parents[1]


def main():
    manifest = json.loads((ROOT / "docs/source-import.json").read_text(encoding="utf-8"))
    for row in manifest["files"]:
        path = ROOT / row["path"]
        assert path.is_file(), row["path"]
        assert hashlib.sha256(path.read_bytes()).hexdigest() == row["import_sha256"], row["path"]
        if path.suffix == ".py":
            ast.parse(path.read_text(encoding="utf-8"), filename=row["path"])
    documents = [ROOT / "README.md", ROOT / "README_ZH.md", ROOT / "ROADMAP.md"]
    documents += list((ROOT / "docs").rglob("*.md"))
    documents += list((ROOT / "research").rglob("README.md"))
    documents += list((ROOT / "examples").rglob("README.md"))
    links = 0
    for path in documents:
        text = path.read_text(encoding="utf-8")
        for url in re.findall(r"\]\(([^)]+)\)", text):
            if url.startswith(("https://", "http://", "#", "mailto:")):
                continue
            relative = unquote(url.split("#")[0])
            assert (path.parent / relative).exists(), (path.relative_to(ROOT), relative)
            links += 1
    tex_inputs = 0
    for path in (ROOT / "docs/derivations").rglob("*.tex"):
        text = path.read_text(encoding="utf-8")
        for name in re.findall(r"\\input\{([^}]+)\}", text):
            target = path.parent / name
            if not target.suffix:
                target = target.with_suffix(".tex")
            assert target.is_file(), (path.relative_to(ROOT), name)
            tex_inputs += 1
    sensitive_paths = re.compile(
        r"/Users/|C:(?:\\+|/)Users(?:\\+|/)|harihi-macos|ssh_config|BEGIN (?:OPENSSH |RSA )?PRIVATE KEY"
    )
    checked = 0
    for base in ("src", "research", "docs", "examples"):
        for path in (ROOT / base).rglob("*"):
            if not path.is_file() or "__pycache__" in path.parts:
                continue
            if path.suffix not in (
                ".py",
                ".cpp",
                ".hpp",
                ".md",
                ".json",
                ".tex",
                ".log",
                ".csv",
                ".svg",
            ):
                continue
            assert not sensitive_paths.search(path.read_text(encoding="utf-8-sig")), (
                path.relative_to(ROOT)
            )
            checked += 1
    result = {
        "source_files": len(manifest["files"]),
        "local_links": links,
        "tex_inputs": tex_inputs,
        "public_text_files": checked,
        "all_passed": True,
        "mac_compute_resumed": False,
    }
    print(json.dumps(result, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
