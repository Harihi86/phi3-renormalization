"""Public packaging, saved-artifact verification and local viewer boundaries."""

import gzip
import hashlib
import json
import shutil
import threading
from functools import partial
from pathlib import Path
from urllib.request import urlopen

import pytest

from phi3_public.counterterms import ReportHandler, ReportServer, build, family, read_json, verify

ROOT = Path(__file__).resolve().parents[1]
EXAMPLE = ROOT / "examples/counterterms"


def test_public_example_recomputes_all_77_cells():
    result = verify(EXAMPLE)
    assert result["all_equal"] and result["checked_cells"] == 77
    assert result["numeric_tolerance_used"] is False


def test_saved_pass_flag_cannot_hide_changed_path(tmp_path):
    for name in ("path-a.json", "path-b.json", "verification.json"):
        shutil.copyfile(EXAMPLE / name, tmp_path / name)
    b = read_json(tmp_path / "path-b.json")
    b["rows"][0]["terms"][0]["coefficient"] = "0"
    (tmp_path / "path-b.json").write_text(json.dumps(b), encoding="utf-8")
    assert read_json(tmp_path / "verification.json")["all_equal"]
    assert not verify(tmp_path)["all_equal"]


def test_existing_output_is_preserved(tmp_path):
    sentinel = tmp_path / "original.txt"
    sentinel.write_text("original evidence", encoding="utf-8")
    with pytest.raises(FileExistsError):
        build(tmp_path)
    assert list(tmp_path.iterdir()) == [sentinel]
    assert sentinel.read_text() == "original evidence"


def test_lossless_compressed_records_support_both_spellings(tmp_path):
    content = {"coefficient": "-1/2", "operators": ["C2m", "C2z"]}
    path = tmp_path / "record.json.gz"
    with gzip.open(path, "wt", encoding="utf-8") as stream:
        json.dump(content, stream)
    assert read_json(path) == content
    assert read_json(tmp_path / "record.json") == content


def test_family_verifies_actual_lower_loop_source_library(tmp_path):
    output = tmp_path / "family"
    result = family(output, 2, 2)
    assert result["all_equal"] and result["checked_cells"] == 14
    assert verify(output)["all_equal"]
    b = read_json(output / "path-b.json")
    reference = next(iter(b["source_library"].values()))
    source_path = output / "library" / f"L{reference['loops']}E{reference['legs']}.json"
    source = read_json(source_path)
    term = next(
        t
        for r in source["rows"]
        if r["id"] == reference["source_id"]
        for t in r["terms"]
        if t["symbol"] == reference["symbol"]
        and t.get("projection_of") == reference["projection_of"]
    )
    term["coefficient"] = "999"
    source_path.write_text(json.dumps(source), encoding="utf-8")
    with pytest.raises(ValueError, match="source projection"):
        verify(output)


def test_http_verification_uses_installed_core_and_closes_cleanly():
    server = ReportServer(("127.0.0.1", 0), partial(ReportHandler, directory=str(EXAMPLE)))
    worker = threading.Thread(target=server.serve_forever, daemon=True)
    worker.start()
    try:
        with urlopen(f"http://127.0.0.1:{server.server_port}/api/verify", timeout=20) as response:
            result = json.load(response)
        assert result["checked_cells"] == 77 and result["all_equal"]
    finally:
        server.shutdown()
        server.server_close()
        worker.join(timeout=5)
    assert not worker.is_alive()


def test_source_import_hashes_and_ledger_categories():
    manifest = read_json(ROOT / "docs/source-import.json")
    assert len({row["path"] for row in manifest["files"]}) == len(manifest["files"])
    for row in manifest["files"]:
        assert hashlib.sha256((ROOT / row["path"]).read_bytes()).hexdigest() == row["import_sha256"]
    ledger = read_json(ROOT / "docs/evidence/precision-ledger.json")
    assert ledger["counts"] == {"strict": 24, "empirical": 48, "unresolved": 90}
    assert len({row["id"] for row in ledger["rows"]}) == 162
    assert ledger["public_export"]["compute_campaign_state"] == "paused by owner"
