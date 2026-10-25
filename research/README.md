# Mathematical research source snapshot

`precision/` contains 326 original Python mathematical modules. The selection
keeps local imports closed, including the exact Gaussian polynomial-binding
helper, while excluding machine-bound job runners, remote-control scripts and
fixed private job bindings. Its algorithms cover exact polynomial/Schur
comparisons, radial reductions, scalar primitives, rational transports,
coordinate-face estimates, pairing and independent reference constructions.

These files retain their original names and source hashes. Some are earlier
candidate versions; their presence is not an acceptance claim. The public
[ledger](../docs/evidence/precision-ledger.json) records actual research status.
In particular, coordinate-face estimates are not a joint multiscale proof,
native arithmetic intervals are not automatically quadrature error bounds,
and candidate80 native tests had not launched when computation was paused.

The source catalogue and dependency closure are in
[source-import.json](../docs/source-import.json). Install optional mathematical
dependencies with `python -m pip install ".[numerics]"` and add
`research/precision` to the import path when studying these standalone modules.
They are source references rather than the default installed command API.

`native/` contains the original generic C++ executor source; `native-precision/`
contains 168 candidate C++/header sources. Generated bindings, original samples
and host build products are not copied. See their individual notices and
[native source scope](native-precision/README.md).

The portable, routinely tested path is the installed `phi3-ct` command and the
`phi3lab` package. The source import does not launch a research campaign.
