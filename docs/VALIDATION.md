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
  unambiguous augmentation target validation.

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
4. rerun from the first `catalog_dd.txt` with prelocation disabled;
5. test pins only when their serial IDs and reference locations are trusted;
6. treat bias-corrected results as an experiment, and retain the uncorrected
   solution as the baseline.

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

The CI workflow is the authoritative Julia runtime gate before tagging a
release. Do not remove this note until `Pkg.test()` has passed on an installed
Julia 1.10+ runtime.
