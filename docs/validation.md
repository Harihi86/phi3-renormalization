# Observed validation for the 0.2.0 source import

Local host: Windows, CPython 3.12. No Mac computation was resumed.

- Full final suite: **59 passed**, no failures or warnings; [actual log](evidence/tests-final.log) and [exit receipt](evidence/tests-final-status.json).
- Fresh installed exact build: **7 rows, 11 columns, 77 exactly equal cells**, 23 nonzero cells and 32 symbolic monomials, with no numerical tolerance.
- Saved path files were independently re-compared by `phi3-ct verify`; changing a coefficient is rejected even when a stored verification file still says pass.
- B was rebuilt with A's catalogue/subgraph/contraction/old insertion-column functions disabled. Signed coefficients, port and half-edge maps, host-loop boundary momenta and lower-source references are covered by mutation tests.
- Complete exact lower-family generation through selected three-loop families is exercised. Four-loop results are historical [status extracts](evidence/historical-families.json), not a new four-loop run.
- Local HTTP verification, preservation of existing output directories, compressed JSON and installed-package imports passed.
- Standard-library one-loop examples at unit mass/invariants returned bubble 0.015592809572239047 and triangle -0.11087590197974906. Order32/64 differences are refinement diagnostics, not rigorous bounds: [actual output](evidence/one-loop.log).
- Wheel and source distribution builds succeeded. Source hashes, document links, TeX inputs, complete HTML cell targets and SVG syntax are also checked; see the evidence JSON files below.
- The original `phi3-demo` command also passed its installed-package smoke run; [actual log](evidence/demo-final.log).

Reproduce from a clean checkout:

```bash
python -m venv .venv
# Activate the new environment.
python -m pip install ".[dev]"
python -m pytest -q
python -m ruff check src tests tools
python -m ruff format --check src tests tools
python tools/check_public_snapshot.py
phi3-ct build --output runs/ct-new
phi3-ct verify --input runs/ct-new
phi3-ct one-loop --output runs/one-loop-new
python -m build
python tools/check_exports.py --input runs/ct-new
```

Tests use fresh repository-local temporary directories, avoiding other tasks'
shared temporary histories. Earlier attempts are recorded in
[validation-attempts.json](evidence/validation-attempts.json): the first
nonisolated package installation lacked its build backend; a first explicit
test scratch directory lacked its parent. Both were corrected. Initial shared
temporary-directory permission warnings and failed-attempt logs are preserved
with absolute local paths redacted; they are not presented as passed runs.

[Package/export check](evidence/package-and-export-check.json) and
[public source check](evidence/public-snapshot-check.json) record bounded
validation. Research native candidates and all historical TeX pages were not
rebuilt or visually reviewed in this update. The 162-integral ledger remains
24 strict analytic / 48 empirical / 90 unresolved. GitHub workflow configuration
alone is not evidence that a remote CI run passed; inspect the actual run for
the published commit.

The initial 0.1 release had 17 passing bare-demo tests on macOS. That historical
result is superseded for this import by the fresh Windows suite above; no
current dependency-security audit or fresh Mac execution is claimed.
