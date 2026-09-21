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

## Interpreting validation

The automated suite checks numerical operations, file formats, constraints and
restart behavior. It does not establish location accuracy or uncertainty
calibration for a new dataset. Retain the exact software revision, TOML and
input files for each reported run, and use the dataset checks above alongside
the solver's convergence diagnostics.

GitHub Actions reports results for each tested revision on Julia 1.10 and the
current stable release. A local `Pkg.test()` run is also recommended before a
new analysis, particularly after changing Julia versions. The comparator's
smoke test generates its own small truth dataset; it does not depend on the
external Yu true-location file.
