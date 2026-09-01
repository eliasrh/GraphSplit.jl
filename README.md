# GraphSplit (Julia)

GraphSplit relocates earthquake catalogs from the synchronized `theta` and
`thetaStd` files produced by [DDSync](https://github.com/eliasrh/DDSync). It is
designed for large catalogs: the event graph is sparse, the nonlinear solve is
matrix-free, and native travel-time tables are memory-mapped instead of loaded
twice into RAM.

This repository is the maintained implementation. It uses only Julia standard
libraries and is driven by one TOML file.

## Quick start

Requirements: Julia 1.10 or later and completed DDSync `theta/` and
`thetastd/` directories.

```bash
cp config/graphsplit_template.toml graphsplit.toml
# Edit graphsplit.toml, then run:
julia --project=. run_graphsplit.jl graphsplit.toml
```

On the first run GraphSplit builds a native travel-time table if it does not
already exist. Later runs inspect the table header, detect Cartesian versus
radial geometry, validate its model hash, spacing, and coverage, and reuse or
rebuild it as needed.

To build or validate only the travel-time table:

```bash
julia --project=. build_travel_times.jl graphsplit.toml
julia --project=. build_travel_times.jl graphsplit.toml --force
```

No package installation is needed. `--project=.` only selects this repository's
Julia project.

## Required inputs

- `catalog.txt`: DDSync catalog format. Defaults assume latitude, longitude,
  and depth in columns 7, 8, and 9, with the positive integer serial event ID
  in the final column.
- `stations.txt`: `STA latitude longitude elevation_m`; the five-column
  `NET STA latitude longitude elevation_m` form is also accepted.
- `theta/theta_<STA>_<P|S>.txt`: `EventID theta refEventID`.
- `thetastd/std_theta_<STA>_<P|S>.txt`: normally
  `EventID std_theta refEventID degree`, with optional extra DDSync columns.
- `vm.txt`: either `depth_km Vp_km_s Vs_km_s` or
  `depth_km radius_km Vp_km_s Vs_km_s`. Repeated depths preserve velocity
  discontinuities.

Rows are matched by the serial event ID, not row order. Reordered or filtered
catalogs are therefore safe as long as their final IDs are unchanged.

## Outputs

The configured output directory contains:

- `catalog_preloc.txt` and `catalog_preloc_filt.txt` from the theta-reference
  star solve;
- `catalog_dd.txt` and `catalog_dd_filt.txt` from sparse graph relocation;
- `catalog_dd_graphmeta.csv`, with component, degree, and support diagnostics;
- `solver_history.csv`, with nonlinear and PCG convergence history;
- `run_summary.toml`, with the run's core counts and final residuals.

When requested with `[uncertainty]`, catalogs remain unchanged and additional
sidecars are written:

- `linerrxyz.txt`, containing regularized linearized x/y/z standard deviations
  and covariance terms for filtered serial IDs;
- `booterrxyz.txt`, containing bootstrap interval/covariance summaries;
- `bootstrap_samples_lon.txt`, `bootstrap_samples_lat.txt`,
  `bootstrap_samples_depth_km.txt`, and `bootstrap_samples_t0_s.txt`, with the
  all-data solution and every station-phase bootstrap realization;
- `bootstrap_replicates.txt`, `bootstrap_block_counts.txt`, and
  `bootstrap_metadata.toml`, recording convergence, exact block multiplicities,
  and the interpretation of the bootstrap.

The catalog date/time and auxiliary columns are preserved. GraphSplit's
internal origin-time adjustments are relative nuisance parameters and are not
written into columns 1–6, matching the MATLAB implementation.

To retain empirical location clouds without changing the catalog format:

```toml
[uncertainty]
method = "bootstrap"

[uncertainty.bootstrap]
replicates = 100
resampling_unit = "station_phase"
write_samples = true
```

The bootstrap is optional because it reruns Stage 2 once per replicate. Use
`method = "linearized"` for the cheaper formal estimate, or `method = "both"`.

## Cartesian and radial travel times

The native builder solves a first-order Godunov eikonal equation by fast
sweeping for P and S waves. Both builders retain station elevation as the third
lookup axis:

- Cartesian: `T(horizontal_range_m, event_depth_m, station_depth_m)`.
- Radial: `T(central_angle_rad, event_depth_m, station_depth_m)` for an
  axisymmetric spherical Earth.

Choose the builder used when no compatible table exists:

```toml
[travel_time]
geometry = "auto"
build_geometry = "radial" # or "cartesian"
```

`geometry = "auto"` is recommended. Geometry is read from the native `.gstt`
header. MATLAB `.mat` lookup tables are intentionally not accepted: the native
format records provenance and can be memory-mapped safely.

## Basic, common-seed, iterated, and pinned runs

The same runner covers all four cases; there are no divergent example
programs to maintain.

- Basic: use the template TOML unchanged apart from paths and grid/graph scale.
- Common seed: set `initialization.mode = "common_centroid"` to start every
  event at the input catalog centroid, or `"common_manual"` to specify one
  latitude, longitude, and depth. Stage 1 is required and separates the events
  before Stage 2 constructs the sparse graph. This permits a complete location
  without single-event catalog hypocenters.
- Iterated: run once, then use the first pass's `catalog_dd.txt` as the second
  pass `io.catalog_file`, select a new output directory, and set
  `run.prelocation = false`. A broad first graph and tighter second graph are a
  useful schedule, not a separate algorithm.
- Pinned: set `gauge.mode = "pin"`, list exact serial IDs in
  `gauge.pin_event_ids`, and choose `gauge.pin_fields`. Optionally provide a
  trusted same-ID catalog in `gauge.pin_reference_catalog`.

Pinning one or a few events after moving only those events is not a reliable
way to translate an entire catalog. The pin itself is exact, but the resulting
large residuals can be strongly Huber-downweighted and the displaced pins can
be isolated when the Stage-2 graph is built. Use common initialization when
the problem is an inadequate seed catalog; use pins only as explicit trusted
constraints within a consistently initialized, connected solution.

Copy-ready TOML fragments are in [examples/README.md](examples/README.md). Every
TOML key—including units, zero behavior, interactions, tuning advice, and
worked examples—is covered in the plain-language
[configuration reference](docs/CONFIGURATION_REFERENCE.md). The algorithmic
reasoning is in [docs/MANUAL.md](docs/MANUAL.md).

## Experimental bias diagnostic

The coherent theta-bias fit is retained as an experimental, disabled option.
It is deliberately absent from the minimal template because there is not yet
clear evidence that applying it improves locations. Every switch is exposed in
`config/graphsplit_complete.toml`. Enable the scan only to write a diagnostic
model; apply a reviewed model in a later run with
`experimental.bias.apply_model_file`.

## Yifan benchmark utility

The independent utility in `benchmark/yifan2025/` compares any relocated
catalog to `truelocs.txt` using the paper's horizontal/depth accuracy,
neighbor-pair precision within 2 km, and point-cloud Chamfer approach:

```bash
julia benchmark/yifan2025/compare_catalogs.jl \
  graphsplit_output/catalog_dd.txt \
  benchmark/yifan2025/input/truelocs.txt
```

It is intentionally outside the GraphSplit module and does not affect a run.

## Documentation and tests

- [Manual](docs/MANUAL.md)
- [Complete configuration reference](docs/CONFIGURATION_REFERENCE.md)
- [File formats](docs/FILE_FORMATS.md)
- [Implementation and MATLAB parity](docs/MATLAB_PARITY.md)
- [Validation notes](docs/VALIDATION.md)

Run the test suite with:

```bash
julia --project=. -e 'using Pkg; Pkg.test()'
```

## Citation and license

Citation metadata are in `CITATION.cff`. GraphSplit is distributed under the
GraphSplit Non-Commercial License v1.0 in `LICENSE`.
