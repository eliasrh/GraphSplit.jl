# Changelog

## 0.4.0 — 2026-09-07

- Applied solved relative origin-time corrections to configured catalog calendar
  fields, including full Gregorian rollover and a no-calendar mode for point
  inputs.
- Added cumulative `catalog_preloc_dxdydzt0.txt` and
  `catalog_dd_dxdydzt0.txt` files in metres/metres/metres/seconds for every
  persistent event ID, with zeros for events never relocated.
- Made a DD-only iterated graph restore its cumulative relative-time state from
  the matching prior shift file, avoiding a misleading high first residual and
  double application of corrected calendar time.
- Added regression coverage for positive and negative date rollover, inactive
  events, invalid calendar rows, ID-reordered shift input, cumulative two-pass
  shifts, and exact restart prediction equivalence.
- Integrated catalog timing, shift-file interpretation, and the sparse-network
  iterated-graph use case into the manual, configuration reference, file-format
  guide, examples, parity notes, and validation guide.

## 0.3.0 — 2026-09-01

- Added an exact fixed-depth constraint for all events or explicit serial-ID
  subsets, with optional per-event depths from a partial reference catalog.
- Added a physical minimum-depth constraint using an active-set
  Gauss-Newton solve that re-solves coupled x/y/t0 updates and releases depth
  when the next trial points back into the admissible interior.
- Added an opt-in unconstrained Stage-1 pilot and reflected bounded restart for
  sparse-station shallow/deep mirror ambiguity; Stage 2 is built only from the
  bounded restart solution.
- Extended travel-time-table compatibility checks to cover the physical bound,
  reflected seed depths, and a configurable forbidden-side pilot margin;
  clamping is disabled during the branch-search pilot.
- Added `depth_constraint_status.csv` and bootstrap bound-active fractions
  without changing any scientific catalog columns.
- Marked one-sided bound-active linearized z covariance as undefined while
  retaining conditional x/y covariance; exact fixed depths retain a documented
  zero conditional z variance.
- Reorganized the README as a first-user workflow and updated the manual,
  configuration reference, file formats, examples, parity notes, and tests.

## 0.2.0 — 2026-09-01

- Added conditional nonlinear block bootstrap uncertainty with station-phase
  groups as the default resampling unit and whole stations as an optional test.
- Preserved every longitude, latitude, depth, and relative-origin-time sample
  in wide serial-ID-keyed tables and recorded exact block multiplicities;
  optional per-replicate catalogs remain pure.
- Added `booterrxyz.txt` percentile or standard-deviation summaries and explicit
  `NaN` handling for events that lose resampled support.
- Added scalable randomized regularized inverse-Hessian covariance in the
  separate `linerrxyz.txt` sidecar.
- Kept all primary/prelocation catalogs unchanged and documented the limits of
  both conditional uncertainty interpretations.
- Added `catalog`, `common_centroid`, and `common_manual` initialization so
  Stage 1 can perform a complete location without single-event seed
  hypocenters; common modes cannot bypass Stage 1.
- Extended lookup-table coverage checks to include both the input catalog and
  an overridden common seed.
- Clarified that hard pinning one or a few displaced events is not a reliable
  catalog-translation method because robust downweighting and graph separation
  can remove their leverage.

## 0.1.3 — 2026-08-13

- Fixed native travel-time table header reads on Julia 1.12 by replacing the
  removed `read(io, UInt8, n)` method with version-portable exact byte reads.
- Added explicit geometry and model-hash assertions to the native-table test.
- Ignored locally generated Julia manifests so one runtime's resolution is not
  accidentally committed to the multi-version package.

## 0.1.2 — 2026-08-13

- Added a plain-language reference for all 117 TOML inputs, with tuning order,
  units, zero behavior, interactions, troubleshooting, and worked examples.
- Added prelocation and relocation damping to the minimal TOML.
- Clarified that a trusted pin catalog may contain only the pinned serial IDs;
  matching is by persistent ID rather than row index or row order.
- Documented exact-count versus average-degree graph augmentation and reject
  enabled configurations with no target or two competing targets.
- Expanded input-file examples and added configuration documentation tests.

## 0.1.1 — 2026-08-13

- Fixed Julia soft-scope warnings in both travel-time/relocation launchers.
- Fixed the benchmark comparator's `neighboring_pairs` scope error.
- Added command-line and benchmark-comparator regression coverage.

## 0.1.0 — 2026-08-12

- First maintained Julia implementation.
- Added native Cartesian and radial memory-mapped travel-time tables.
- Added TOML-driven two-stage robust sparse relocation.
- Added exact dependency-free kNN event graph and matrix-free PCG solver.
- Added full/prelocation filtered catalogs, graph metadata, and solver history.
- Added exact serial-ID pinning and optional trusted same-ID pin catalog.
- Retained the coherent theta-bias scan as an experimental disabled feature.
- Added the standalone Yifan-style benchmark comparator.
