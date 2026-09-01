# Changelog

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
