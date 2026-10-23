# Declared scientific conventions

## Theory and power counting

The reference interaction is `-g φ³ / 3!` in the Minkowski Lagrangian, with kinetic
term `+(∂φ)²/2` and mass term `-m²φ²/2`. The metric is `(+,-,...,-)` in d dimensions;
in four dimensions this is `(+,-,-,-)`. The graph algorithms themselves are metric
independent. Power counting uses Euclidean UV scaling of a scalar propagator
`1/(k²+m²)`, with `m > 0` when future zero-momentum subtractions require infrared
control. The current executable does not evaluate a Minkowski or Euclidean integral.

The demonstrated dimension is `d=6`, where `[φ]=(d-2)/2`, `[g]=(6-d)/2` and
`ω=dL-2I=6-2E` for connected bare trivalent graphs. Numerator insertions or typed
counterterm vertices require different power counting and are outside this demo.

## Graphs and weights

Internal edges are an indexed multiset. Subgraphs are sets of edge indices, not
induced vertex subsets. Every omitted incident half-edge becomes a separate
boundary port. Two-point and three-point external legs are amputated, distinctly
labelled and held fixed in canonicalization and symmetry calculations. All
external momenta are incoming with sum zero; no source-J attachments are counted.

The generator attaches each external leg to a distinct vertex in the supported
1PI scope. Parallel edges are retained; self-loops are allowed. Vacuum graphs are
excluded. A self-loop can contribute to the one-vertex one-point case, while
bridges are excluded by the 1PI test.

`weight = 1 / |Aut(G)|` is a positive rational combinatorial weight, with vertex,
parallel-edge and self-loop half-edge automorphisms included. It excludes powers
of g, factors of i, propagators and integration normalization. A fixed-label
bubble has weight 1/2; a fixed-label triangle has weight 1. Unfixed-label factors
must not be substituted into the fixed-label demonstration.

## Renormalization boundary

No subtraction prescription is executed by this release. A future massive
Euclidean momentum-subtraction implementation will need to declare Taylor
expansion points, subtraction degrees, overall versus proper-subgraph terms,
local operator basis and counterterm signs explicitly. A formal statement
`R(G)=(1-t_G) Rbar(G)` is a definition under a declared prescription, not an
implemented numerical result here. A graph-level UV degree is not evidence that
a subtracted integral converges.

BPH recursion, Henge decompositions and Zimmermann forest representations must
be compared under matched conventions. Henge components need not all be UV
divergent. Claims about all-loop finiteness require scale bounds uniform in
subgraph external momenta; they cannot follow from a finite list of computations.
