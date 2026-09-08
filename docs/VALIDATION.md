# Validation status

## Automated checks in this repository

`test/runtests.jl` covers:

- local latitude/longitude coordinate round trips;
- spherical central-angle derivatives against finite differences;
- exact k-d-tree neighbors against brute force;
- constant-velocity Cartesian and radial travel-time derivatives;
- the fast-sweeping field against straight-ray constant-velocity times;
- native table header, memory-map, interpolation, and analytic gradients;
- serial-ID observation construction after catalog reordering;
- forward/adjoint consistency of the matrix-free Jacobian;
- PCG recovery for a known positive-definite system;
- launcher argument parsing and a zero-error benchmark-comparator smoke test;
- complete TOML documentation coverage, visible minimal damping, and
  unambiguous augmentation/uncertainty validation;
- common-centroid/manual initialization, its Stage-1 requirement, and lookup
  coverage of overridden manual seeds;
- Gregorian origin-time rollover, disabled/invalid calendar handling, cumulative
  x/y/z/t0 sidecars, serial-ID restart joins, and restart prediction equivalence;
- exact fixed-depth parameter removal, positive/negative vertical conventions,
  reflected-depth feasibility, and an active-set crossing/release regression;
- exact station-phase bootstrap multiplicities, sample quantiles/covariances,
  and serial-ID-keyed uncertainty sidecar formats.

Run it with `julia --project=. -e 'using Pkg; Pkg.test()'`. GitHub Actions runs
the same suite on Julia 1.10 and the current stable Julia release.

## Numerical expectations

The fast-sweeping solver is first order. In a homogeneous model, interpolation
error should decrease with grid spacing; exact agreement away from grid rays is
not expected on a coarse grid. Relocation gradients are derivatives of the
stored trilinear table, so solver values and Jacobians remain internally
consistent.

For a new field dataset, validate in this order:

1. build a table and inspect its reported grid size and coverage;
2. run a basic solution and check convergence and graph support;
3. compare Cartesian and radial solutions if aperture makes curvature relevant;
4. when testing an iterated graph, rerun from the first `catalog_dd.txt` with
   prelocation disabled and its matching `catalog_dd_dxdydzt0.txt` supplied as
   `io.restart_shift_file`;
5. test pins only when their serial IDs and reference locations are trusted;
6. if using a physical depth bound, compare an unconstrained diagnostic run,
   inspect reflected IDs and final bound-active events, and do not interpret a
   boundary pile-up as depth resolution;
7. for known-depth explosions, verify fixed IDs/depths in
   `depth_constraint_status.csv` and remember that their z uncertainty is
   conditional on the constraint;
8. treat bias-corrected results as an experiment, and retain the uncorrected
   solution as the baseline.
9. for a published uncertainty result, inspect bootstrap sample clouds and the
   valid-replicate fraction rather than relying only on covariance summaries.

## Current execution boundary

The source was parsed as Julia and the configuration files were parsed as TOML
in the build environment used to assemble this release. That environment did
not contain a Julia runtime, so the repository test suite could not be executed
there. A subsequent user-run Julia execution on macOS completed the bundled
1,000-event Cartesian example end to end, including native table construction
and both relocation stages. The first full `Pkg.test()` run on Julia 1.12.3
passed 259 checks and exposed a removed byte-vector `read` method in the native
table reader; v0.1.3 replaces it with the version-portable `read!` API.
Independent translations of the homogeneous fast-sweeping equations
converged in two outer sweeps on the test grid; the Cartesian coordinate-axis
errors were below `4e-16 s`, and the radial same-depth surface error against the
constant-velocity chord time was below `3e-8 s`. An independent evaluation of
the bundled benchmark seed also produced the smoke-test metrics recorded in the
benchmark README.

For v0.2.0, GitHub Actions completed all 315 checks successfully on Julia 1.10
and the current stable Julia release on 2026-09-01. The two-version CI workflow
remains the authoritative release gate. The complete v0.4.0 working tree passes
386 checks locally on Julia 1.10.10, including the depth-constraint and catalog
time/shift restart regressions.
