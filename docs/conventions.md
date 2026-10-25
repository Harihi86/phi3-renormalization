# Executable subtraction conventions

## Action, signs and graph classes

All subtraction/numerical formulas use six-dimensional Euclidean signature,
`D(q)=1/(q²+m²)`, `m²>0`, and one `d⁶q/π³` measure per loop. The action contains
`+gφ³/3!`, `+δm²φ²/2!`, `+δZ(∂φ)²/2`, `+δgφ³/3!` and a linear term `aφ`.
Expansion of `exp(−S_int)` gives vertices `−g`, `−δm²`, `−δZ P²`, `−δg` and `−a`.
The positive weight `1/S` excludes coupling, signs, propagators and measure.

The older bare demo's Minkowski description does not insert factors of `i`
into this Euclidean implementation; its graph operations are metric-independent.

`Z[J]=exp(−Σ_r λ_r/r! ∫(δ/δJ)^r) exp(JDJ/2)` is normalized by `Z[0]`.
Connected functions come from `W[J]=log(Z[J]/Z[0])`; differentiate sources before
setting `J=0`. Amputation removes external free propagators. 1PI additionally
excludes internal bridges. Functional audits retain direct source pairings and
vacuum components before these restrictions.

Edges are indexed multisets; divergent subgraphs are edge subsets, not induced
vertex subsets. Every omitted incident half-edge supplies a boundary port.
External labels are fixed pointwise, with all momenta incoming and summing to
zero. Self-loops and tadpoles are retained. No normal-ordering/zero-VEV condition
is imposed. The finite linear term is unspecified; one-point insertion hosts
have bridges and do not enter the target amputated 1PI three-point family.

For connected bare trivalent graphs, `3V=2I+E`, `L=I−V+1` and `ω=6−2E`.
Typed counterterm graphs use their own valences and coupling orders.

## Projectors and independent generation

Two-point Taylor subtraction is total vector degree ≤2 at independent boundary
momentum `P=0`, with local basis `1,P²`. It is not degree two in `P²`.
Three-point subtraction is degree zero in two independent boundary vectors.
The third vector is minus their sum.

Inside each selected subgraph, retain its inner loop variables and replace
boundary vectors by `tP`; extract powers through `t^ω`. Boundary momenta can
contain host loop variables. Host propagators are not scaled by this operation.

Path A recursively constructs `C(γ)=−T_γ R̄(γ)` and applies
`R(G)=(1−T_G)R̄(G)`. Nested operators act inside out; disjoint subgraphs may be
multiplied, overlapping subgraphs may not. Henge and forest checks audit this
actual recursion.

Path B independently enumerates typed action sectors and Wick pairings,
constructs its own lower-loop counterterm sources and glues their ports.
No target catalogue, contracted insertion column or target coefficient is
supplied by A. Graph representation and canonicalization helpers are shared.

Every comparison cell is an exact rational operator polynomial. Signs, vertex
types, fixed labels, half-edge bijections and loop-momentum changes are checked.
Structural zeros remain explicit; numerical tolerances are not used.

Local divergent integrals remain symbolic under a common UV regulator until
subtractions are combined. Table equality is neither a numerical evaluation of
those separate divergent terms nor a general all-loop convergence proof.
