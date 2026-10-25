# From the generating functional to the correspondence table

[中文](theory-zh.md) · [Conventions](conventions.md) · [L2E3G4](L2E3G4.md)

## 1. Count before restricting the graph class

The interaction exponential supplies the species factorial and a factor
`1/(r!)` per r-valent vertex. The free Gaussian supplies half-edge pairings.
`ct_wick.enumerate_wick` records their exact multiplicity, vacuum components,
source pairings and the reason a class is accepted or excluded. Connectedness,
amputation and the internal-bridge 1PI restriction are separate operations.
Fixed external labels remain fixed under graph isomorphisms. Indexed parallel
edges and self-loop half-edges contribute their own automorphisms.

The signed coefficient is the positive rational Wick weight multiplied by
vertex signs and coupling powers. This separation prevents mistaking `1/S`
for the full amplitude coefficient. `ct_typed` extends the enumeration to
counterterm vertex species and higher families.

## 2. Path A actually recurses

`ct_bph.recursion` constructs each proper edge-set divergent subgraph, its
boundary ports, its recursively subtracted interior and the contracted host.
It records `C(γ)=−T_γ R̄(γ)` and `R(G)=(1−T_G)R̄(G)` rather than inferring them
from a finished insertion table. Compatible routing distinguishes inner and
outer loop bases. The Taylor operation scales only boundary dependence inside
the chosen subgraph. A host loop variable entering that boundary is retained
in the host.

Disjoint subgraphs may form a spinney. Nested counterterms arise from the
recursive child, with the inner projector applied first. Overlapping subgraphs
are retained individually but never multiplied as one compatible spinney.
Henge edge maps and a separately expanded forest oracle check the recursion.
`ct_series_a` extends this representation through four loops.

## 3. Path B starts with the action

`ct_action.build_path_b()` takes no target catalogue. It enumerates typed
action sectors at the required coupling order, Wick-pairs each sector, applies
the connected/1PI restrictions, and independently builds lower-loop local
sources. Every source is glued into the appropriate counterterm vertex through
explicit port bijections. Canonical vertex, edge and half-edge maps identify
the expanded original graph. Multiplicity is divided by the port factorial;
it is not added as an extra unlabelled-leg symmetry factor.

`ct_series_b.ActionLibrary` recursively constructs its own lower-loop local
projections. Stored source-library references are checked against those actual
families. Its inputs are types, coupling order and action rules, not Path A's
subgraphs, contractions or coefficients. Sharing graph representation and
canonicalization does not make the generation paths independent numerical
integrators; that is a different question.

## 4. Compare the same mathematical objects

`ct_verify` and `ct_series_verify` compare rational operator polynomials in each
row and column, including structural zeros. They validate graph labels/types,
Taylor degree, source identity, half-edge maps and loop-basis transformations.
A numerical tolerance cannot turn a failed equality into a pass. Regression
tests alter signs, multiplicities, labels and mappings and require rejection;
they also disable A's generators while rebuilding B.

The final local column is obtained by applying minus the overall projector to
the bare plus lower-counterterm terms. The mass and kinetic overall columns in
a two-point family remain distinct. The explicit [bubble insertion example](L2E3G4.md)
shows how this changes signs without deleting the host loop momentum.

## 5. Numerical evaluation is a separate layer

`automatic_bph` compiles graph incidence and routing into Gaussian Taylor
programs; `finite` and the topology-specific modules retain independent finite
reference constructions. The later research sources perform exact polynomial
binding, analytic radial integrations, rational changes of variables and
interval primitive evaluation. An exact change of variables, a coordinate-face
estimate, an arithmetic interval and a whole-integral error bound are different
claims. Their current acceptance inventory is in [status.md](status.md).

The TeX [compiler project](derivations/workbench/automatic-analysis-project/main.tex)
and [later derivations](derivations/README.md) contain the longer original
arguments. No general all-loop theorem is inferred from finite test coverage.
