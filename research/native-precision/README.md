# Native mathematical source archive

These 168 original C++ and header files preserve the research implementations,
including independent interval arithmetic kernels and experimental variants.
Original filenames and hashes are in `docs/source-import.json`.

This is source material, not a prebuilt backend. Some translation units require
generated graph-specific tables or frozen moment headers from research jobs,
and external FLINT/Arb or Boost headers and libraries. Those bindings, binaries,
external dependency trees and machine launch scripts are not bundled. Building
every historical variant with a single wildcard command is not supported.

Several variants are superseded or experimental. Consult the derivation and
the recorded ledger; presence here is not evidence that a numerical candidate
passed. No native job is automatically launched by package installation, test
collection or the public exact CLI.
