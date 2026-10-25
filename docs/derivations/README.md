# Editable mathematical derivations

This directory contains 53 original LaTeX sources. Four historical workbench
projects cover two-loop analysis, generic recursive Gaussian BPH compilation,
all supported two-loop numerical topologies and extreme-scale formulations.
Thirty further project entry points cover the later precision research.

Each project has its own `main.tex`. Supporting `\input` sources are retained.
The Chinese workbench projects use `ctexart` and Fandol fonts; compile them with
XeLaTeX (or a compatible Unicode TeX environment). Other projects use their
declared standard LaTeX packages. The generated exact example has separate
Chinese/English master files under `examples/counterterms/latex`.

Historical version labels and progress statements describe the original
document date. The current acceptance inventory is in
[status.md](../status.md). In particular, a coordinate-face bound or a local
native interval check does not establish a new whole-integral certificate.
No claim is made that this import recompiled and visually reviewed every page
of all historical projects.

Useful entry points: [recursive compiler](workbench/automatic-analysis-project/main.tex),
[two-loop analysis](workbench/full-analysis-project/main.tex),
[projective forests](projective-forest-remainders/main.tex),
[paired face estimates](paired-face-orders/main.tex),
[mixed/cyclotomic primitives](mixed-and-cyclotomic-derivation/main.tex),
[hybrid normalization and recovery](hybrid-primitive-recovery/main.tex).
