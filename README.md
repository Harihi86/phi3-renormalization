# φ³ Renormalization

[中文说明](README_ZH.md) · [Theory and source map](docs/theory.md) · [Worked example](docs/L2E3G4.md)

Research software for renormalization in six-dimensional scalar φ³ theory,
with reproducible graph generation, counterterm construction and exact verification.

**0.2.0** imports the current Phi3Lab computational core: generating-functional
counting, graph isomorphisms, exact rational weights, recursive BPH subtraction,
independently generated counterterm-action graphs, and cell-by-cell comparison.
Existing numerical kernels, research mathematics and editable LaTeX derivations
are included, with their validation scope stated below.

## Install and run

Python **3.11+**. Exact combinatorics and the finite one-loop examples require
only the standard library.

```bash
python -m venv .venv
# Activate .venv using your platform's normal command.
python -m pip install ".[dev]"
phi3-ct build --output runs/counterterms
phi3-ct verify --input runs/counterterms
phi3-ct serve --input runs/counterterms --language en
```

Open the printed localhost address; use `--language zh` for Chinese.
Generation requires a new output directory. These commands do not connect to
SSH, launch remote workers or run a large integration campaign.

The default build generates both paths afresh and produces the complete
**7-row × 11-column = 77-cell** two-loop three-point table, including structural
zeros. Outputs include `path-a.json`, `path-b.json`, `table.csv`,
`verification.json`, full Wick audits, Chinese/English HTML, SVG and editable
TikZ/LaTeX. [Compact reproducible example](examples/counterterms/README.md).

![Graph with a two-point subdivergence](examples/counterterms/figures/L2E3G4.svg)

## Implementation and principles

| Component | Source |
|---|---|
| Bare multigraphs and isomorphisms | [`graphs.py`](src/phi3lab/graphs.py), [`higher_graphs.py`](src/phi3lab/higher_graphs.py) |
| Generating-functional and typed Wick counting | [`ct_wick.py`](src/phi3lab/ct_wick.py), [`ct_typed.py`](src/phi3lab/ct_typed.py) |
| Path A: recursive `C = −T R̄` and Henge checks | [`ct_bph.py`](src/phi3lab/ct_bph.py), [`ct_series_a.py`](src/phi3lab/ct_series_a.py) |
| Path B: independent action expansion and port gluing | [`ct_action.py`](src/phi3lab/ct_action.py), [`ct_series_b.py`](src/phi3lab/ct_series_b.py) |
| Exact structural/coefficient/momentum validation | [`ct_verify.py`](src/phi3lab/ct_verify.py), [`ct_series_verify.py`](src/phi3lab/ct_series_verify.py) |
| Numerical implementation | [`automatic_bph.py`](src/phi3lab/automatic_bph.py) and retained `finite`, `primitive`, `nested_bubble`, `nested_triangle`, `overlap`, `two_loop` modules |
| Further mathematical sources | [326 Python research modules](research/README.md), [168 native source/header files](research/native-precision/README.md), [53 LaTeX sources](docs/derivations/README.md) |

Path B does not obtain insertion columns from A's contractions or
`renormalization.insertion_columns`. Regression tests forbid those calls while
rebuilding B. The paths share graph data structures and isomorphism utilities;
their graph generation and coefficient derivations are independent.

The engine uses the massive **six-dimensional Euclidean** action with
`+gφ³/3!`, propagator `1/(q²+m²)`, `m²>0`, and loop measure `d⁶q/π³`.
Bare vertices are `−g`; counterterm vertices are minus their action coefficients.
Two-point Taylor subtraction has **total vector degree 2**, split into mass
`P⁰` and kinetic `P²`; three-point subtraction has degree zero. Tadpoles are
retained and the finite linear counterterm condition is unspecified.

[L2E3G4](docs/L2E3G4.md) traces the bare coefficient `−1/2`, mass and kinetic
insertion coefficients `+1/2`, and overall subtraction. Its subgraph boundary
momentum is a host loop variable: only the subgraph is Taylor-scaled.
[Complete conventions](docs/conventions.md).

## Scope and numerical status

```bash
phi3-ct family --loops 2 --legs 2 --output runs/L2E2
phi3-ct verify --input runs/L2E2
phi3-ct one-loop --output runs/one-loop
python -m pytest -q
```

The family engine supports `--loops 1..4` and `--legs 2|3`. Higher-loop generation
is opt-in and costs more. Historical complete-family records have these sizes;
the four-loop campaign is not rerun by this import:

| Family | Bare rows | Columns | Exact cells |
|---|---:|---:|---:|
| L1E2 | 1 | 3 | 3 |
| L1E3 | 1 | 2 | 2 |
| L2E2 | 2 | 7 | 14 |
| L2E3 | 7 | 11 | 77 |
| L3E2 | 10 | 40 | 400 |
| L3E3 | 56 | 162 | 9,072 |
| L4E2 | 64 | 385 | 24,640 |
| L4E3 | 540 | 2,553 | 1,378,620 |

Counts use fixed external labels and amputated connected 1PI graphs; they are
not full generating-functional term counts. Missing external reference diagrams
remain pending graph-by-graph comparison.

The separate [162-integral research ledger](docs/evidence/precision-ledger.json)
records **24 strict analytic results, 48 empirical dual-kernel acceptances and
90 unresolved cases**. Empirical acceptance is not a rigorous quadrature bound;
exact combinatorial equality is not integral accuracy. Heavy research
computation remains **paused**. This source import adds no accepted integrals.
[Status and limitations](docs/status.md).

Optional numerical dependencies: `python -m pip install ".[numerics]"`.
The wheel contains the 55-module computational core. Research sources and TeX
are supplied in the checkout and source distribution. Native research bindings,
machine-specific job runners and bulk raw samples are not bundled.

## Verification and provenance

[Observed validation](docs/validation.md) separates fresh tests from historical
evidence. [Source import manifest](docs/source-import.json) records original and
public hashes. The earlier `phi3-demo` and `phi3_public` API remain available.

```bash
python -m ruff check src tests
python -m ruff format --check src tests
python -m pytest -q
python -m build
```

Imported research code retains its original formatting and is excluded from
format-only rewrites; exact regression tests exercise its behavior.
[Publication scope](docs/publication.md) · [Roadmap](ROADMAP.md).

References: Srednicki, *Quantum Field Theory*
([author's page](https://www.physics.ucsb.edu/~mark/qft.html)); Kennedy,
*A Simple Proof of the BPH Theorem* ([hep-th/9612113](https://arxiv.org/abs/hep-th/9612113));
Kennedy, Binoth and Rippon, *Automating Renormalization of Quantum Field Theories*
([arXiv:0712.1016](https://arxiv.org/abs/0712.1016)).
No new general all-loop convergence theorem is claimed.

MIT-licensed original software; external books, papers and dependency sources
are cited or installed separately, not vendored.
