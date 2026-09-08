# Implementation notes and MATLAB parity

The Julia implementation keeps the mathematical core and file contract of the
provided MATLAB `graphsplit_theta_dd_v3_radial_v1.1` working version while
making the public interface smaller and safer.

## Preserved behavior

- Two-stage theta-reference prelocation followed by sparse theta DD relocation.
- DDSync theta degree and uncertainty filtering.
- Same-reference requirement for every Stage-2 theta pair.
- Cartesian and radial axisymmetric travel-time construction from a layered
  1-D P/S model, including station elevation.
- Robust Gauss–Newton/Huber IRLS with eventwise 4×4 preconditioning.
- Mutual local event graphs, degree/radius limits, optional augmentation, and
  zero-mean or exact pinned gauges.
- Primary `catalog_preloc*` and `catalog_dd*` output names.
- The coherent theta-bias scan is available but experimental and disabled.
- Regularized randomized inverse-Hessian uncertainty is retained as a separate
  `linerrxyz.txt` sidecar rather than appended to the catalog.
- A maintained station-phase block bootstrap preserves every sampled location
  in wide, serial-ID-keyed tables.
- Catalog, common-centroid, and manually specified common starting locations
  provide the HypoDD-style initialization choices without separate runners.

## Deliberate Julia changes

- One recursive TOML configuration replaces MATLAB config scripts and separate
  run examples.
- Every theta/catalog join uses the final serial event ID. This remains correct
  after catalog row reordering and filtering.
- Exact dependency-free k-d-tree search replaces any dense all-pairs distance
  path. The production solver is matrix-free PCG; direct construction is only a
  bounded small-test fallback.
- Native `.gstt` tables store provenance and geometry in their header and map P
  and S cubes from disk. Legacy MATLAB `.mat` tables are not read.
- Travel-time value and gradient use one consistent trilinear interpolant. The
  MATLAB prototype could construct derivative cubes and use MATLAB-specific
  interpolation variants; the Julia form avoids two extra large arrays and is
  deterministic across installations.
- Trusted pin catalogs must share serial IDs. MATLAB's optional heuristic
  time/location/magnitude matching is not part of the maintained core because
  silent pin mismatches are more damaging than requiring an explicit ID join.
- Common initialization is validated to require Stage 1, so Stage 2 never
  constructs a nearest-neighbor graph directly from coincident seeds.
- Physical minimum depth uses an active-set Gauss-Newton update rather than
  travel-time-table clamping. An explicit reflected Stage-1 restart can search
  the admissible counterpart of a sparse-network mirror solution before graph
  construction.
- Exact fixed-depth constraints may apply to all events or serial-ID subsets
  independently of the gauge; x, y, and relative origin time remain free.
- Constraint status and bootstrap bound-active fractions are separate
  serial-ID-keyed sidecars, preserving the primary catalog format.
- Relocated catalogs apply relative origin-time corrections to configured
  calendar fields with full rollover, matching the practical HypoDD convention.
  Cumulative local x/y/z/t0 shifts are also written in serial-ID-keyed sidecars;
  the MATLAB working version kept `t0` only in memory.
- A DD-only pass must pair a previous relocated catalog with its cumulative
  shift sidecar. This reconstructs the time state used in theta predictions and
  prevents an apparently relocated restart from silently resetting `t0`.

## Intentionally excluded experimental diagnostics

The MATLAB working tree also contains prototype group-influence, SAC, and
plotting utilities. Those are not in the authoritative core release because
they are not required for relocation and several depend on interpretation or
environment-specific choices. The uncertainty methods are now maintained Julia
features, with explicit output and interpretation contracts rather than the
MATLAB prototype's appended catalog columns.
