# Observed validation

Date: 2026-10-23. Host: macOS arm64, CPython 3.12.0.

- `python -m pytest -q`: **17 passed**, 0 failed, 0 errors, 0 skipped.
- `python -m ruff check src tests`: PASS.
- `python -m ruff format --check src tests`: PASS.
- Executed demonstration: one-loop bubble/triangle; seven two-loop three-point topologies in declared convention.

`python -m build`: PASS for wheel and source archive. A separate fresh virtual
environment installed the wheel and ran the documented demo successfully.
Dependency audit of the validation environment after updating its installation
tool: no known vulnerabilities reported by pip-audit on 2026-10-23. GitHub
Actions configuration is not itself evidence that a remote run passed; consult
the actual Actions runs for the published commit.

No amplitude integration, fresh typed counterterm generation or full BPH theorem validation was run for this public subset.
