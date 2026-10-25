# Scientific status at the source import

The latest preserved research ledger is `precision-outcomes-index-59`:
24 strict analytic entries, 48 empirical dual-kernel entries and 90 unresolved
entries, totaling 162 unique original graph IDs. Its compact public export is
[precision-ledger.json](evidence/precision-ledger.json). Original summary and
supervisor hashes are retained there; private job payloads are not bundled.

Strict analytic and empirical categories are distinct. The empirical gates
compare independent constructions with fixed sample/refinement criteria; they
do not provide a rigorous quadrature error bound. Exact combinatorial equality
is checked separately and does not turn an empirical integral into a strict one.

At the owner's pause, eight active compute workers and their native children
were stopped. Seven diagnostic jobs were incomplete. Candidate78's diagnostic
finished, but its one-time selection was empty; production did not start.
Candidate80's proof and native build finished, while its native test and test
supervisor directories did not exist: those tests were not actually launched.
The hybrid candidate79 diagnostic remained incomplete. No partial diagnostic
result is used to select production or increase acceptance counts.

This import runs local installation, bounded exact tests and small one-loop
examples only. It does not restart Mac jobs, replay uncertain operations,
alter original thresholds or add integral acceptances.

The exact four-loop source and historical family sizes are included. Default
public tests exercise lower families and selected three-loop exact closure;
the expensive complete four-loop generation is opt-in. Missing external
reference diagrams remain pending external comparison. The finite linear
tadpole condition is unspecified. No all-loop uniform convergence theorem is
inferred from these finite computations.
