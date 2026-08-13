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
- Catalog time fields are preserved rather than modified by relative `t0`.
- The coherent theta-bias scan is available but experimental and disabled.

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

## Intentionally excluded experimental diagnostics

The MATLAB working tree contains prototype uncertainty, group-influence, SAC,
and plotting utilities. They are not in the authoritative core release because
they are not required for relocation, several depend on interpretation choices,
and including them would enlarge an otherwise standard-library-only API. The
CSV graph/support diagnostics and experimental bias scan cover the maintained
diagnostic surface. This is an explicit scope decision, not an accidental
partial port.
