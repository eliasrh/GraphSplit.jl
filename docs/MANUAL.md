# GraphSplit manual

## 0. Before you do anything else

GraphSplit does not read `dt.cc` directly. First run
[DDSync](https://github.com/eliasrh/DDSync) on the differential-time data and
retain its `theta/` and `thetastd/` directories. See Elías Rafn Heimisson and
Yifan Yu, “DDSync: Graph-Based Denoising of Differential Travel-Time
Observations with Applications to Pick Reconstruction and Path-Difference
Tomography,” *Seismological Research Letters* (2026),
[doi:10.1785/0220260086](https://doi.org/10.1785/0220260086).

GraphSplit operates on the synchronized station-phase arrival-time potentials
rather than returning to the much larger `dt.cc` representation. This is a key
part of its computational efficiency. Supply the same starting catalog used by
DDSync, or a reordered/filtered catalog that preserves the persistent serial
IDs in its final column. The serial ID—not the row number—is the join key.

## 1. What the two stages solve

DDSync represents a kept station–phase observation by a synchronized potential
`theta`. For event `i`, reference event `r`, and one station–phase group,
GraphSplit predicts

```text
theta_i - theta_r = (t0_i - t0_r) + (T_i - T_r).
```

Stage 1 uses these reference-centered star observations to improve the seed
catalog. Stage 2 builds a sparse local event graph and uses only graph edges
whose two theta entries share the same DDSync reference:

```text
theta_i - theta_j = (t0_i - t0_j) + (T_i - T_j).
```

Both stages use robust Gauss–Newton/Huber IRLS. The Jacobian is never assembled
as a general sparse matrix. Its forward and adjoint products are evaluated
observation by observation, and the normal equations are solved by PCG with a
per-event 4×4 block-Jacobi preconditioner. Memory is linear in events,
observations, graph edges, and theta entries—not quadratic in the event count.

## 2. Prepare a run

Keep a project directory such as:

```text
my_run/
  graphsplit.toml
  catalog.txt
  stations.txt
  vm.txt
  theta/
  thetastd/
```

Copy `config/graphsplit_template.toml` into it. Paths in the TOML are resolved
relative to that TOML, so the runner can be launched from anywhere:

```bash
julia --project=/path/to/GraphSplitJulia \
  /path/to/GraphSplitJulia/run_graphsplit.jl \
  /path/to/my_run/graphsplit.toml
```

GraphSplit recursively merges the supplied file over internal defaults. It
warns about unknown keys rather than silently ignoring misspellings. The
complete default surface is in `config/graphsplit_complete.toml`. The
plain-language [configuration reference](CONFIGURATION_REFERENCE.md) explains
every key, including units, zero behavior, parameter interactions, and worked
examples.

## 3. Starting without single-event catalog locations

The default `initialization.mode = "catalog"` uses each input hypocenter as its
own seed. When those hypocenters are unavailable or intentionally excluded,
GraphSplit can begin every event at one common location and perform the
complete location from DDSync theta data:

```toml
[initialization]
mode = "common_centroid"
```

`common_centroid` uses the mean local x, y, and depth of the input catalog.
Alternatively, specify an approximate cluster location explicitly:

```toml
[initialization]
mode = "common_manual"
latitude = 64.0200
longitude = -21.2100
depth_km = 3.0
```

These are starting guesses, not fixed locations or additional observations.
Stage 1 uses the station geometry, travel-time model, and theta-reference
differences to separate the coincident events and determine their locations.
GraphSplit therefore requires `run.prelocation = true` for both common modes.
Stage 2 is built only after Stage 1; constructing a nearest-neighbor graph
directly from coincident points would give arbitrary pairs.

The input catalog is still required for persistent event IDs and preserved
date/time/auxiliary columns. With a common mode, its latitude, longitude, and
depth columns no longer provide individual seeds. The native lookup-table
coverage check includes both the input catalog footprint and the overridden
common seed.

## 4. Travel-time model

### Native table lifecycle

With `travel_time.type = "lookup"`, GraphSplit:

1. reads the velocity model and computes its SHA-256 digest;
2. opens an existing `.gstt` table and reads its geometry and axes;
3. checks model digest, requested geometry, grid spacing, event/station depth,
   and source–station distance coverage;
4. reuses it if compatible, otherwise rebuilds it when `auto_rebuild = true`.

The table stores Float32 travel times for P and S. Coordinates and nonlinear
calculations remain Float64. The two cubes are memory-mapped, so the operating
system pages only the portions used. Set `lookup.maximum_table_gib` to a
resource limit; a larger table is refused unless `allow_large_table = true`.

### Geometry choice

Cartesian geometry is the local flat-Earth approximation and is appropriate
when the aperture is small enough for the intended precision. Radial geometry
uses exact spherical source–station central angle in the lookup and evaluates
the spherical chain rule in the relocation Jacobian. It is the safer default
for wider networks.

`horizontal_step_m` is a range increment in Cartesian mode and a surface-arc
increment, converted internally to radians, in radial mode.

### Vertical convention

The default is depth positive down. Station lookup depth is `-elevation_m`, so
stations above sea level have negative station depth. `event_vertical` can
adapt unusual catalogs, but it should normally remain `positive_depth`.

### Constant-velocity mode

`travel_time.type = "constant"` bypasses table construction and uses
`vp_ms`/`vs_ms`. This is useful for tests and controlled synthetic problems,
not as a default field model.

## 5. Observation filters

- `minimum_theta_degree`: requires the DDSync degree column when it exists.
- `maximum_sigma_s`: optional upper limit on the combined theta uncertainty;
  zero disables it.
- `minimum_observations_per_pair`: requires this many station–phase
  contributions on an event graph edge.
- `minimum_observations_per_event`: iteratively removes observations attached
  to weakly supported events; zero disables it.
- `minimum_component_size`: rejects DD observations in small geometry-graph
  components; zero disables it.

When no thetaStd file exists, GraphSplit uses the configured sigma floor and
cannot apply a degree filter for that group.

## 6. Event graph

`graph.neighbors` selects exact k-nearest candidates in `xy`, `xyz`, or
depth-scaled `xyz_scaled` coordinates. `mutual = true` keeps an edge only when
each event selects the other, which limits hubs. `maximum_degree` is then
enforced globally.

`maximum_distance_km` adds a hard radius. Zero disables it. The optional
minimum-degree pass can restore under-connected events from the candidate set.
`ensure_connected` can add available candidate edges between components, but it
cannot bridge farther than `maximum_bridge_distance_km` when that limit is set.

Spectral-style augmentation is available in the complete TOML. It uses an
approximate Fiedler vector to prioritize a small number of graph-stiffening
edges. Leave it disabled until the base graph diagnostics show a need.

## 7. Solver and gauges

The default gauge constrains only mean origin-time adjustment. That leaves the
catalog's spatial frame determined by the travel-time geometry and stations,
while removing the exact relative-time null mode.

For a weakly anchored first pass, `zero_mean = "xyzt0_scaled"` adds weak mean
spatial-update constraints after converting metres to seconds with
`reference_velocity_ms`. Avoid the unscaled form unless its units are
deliberately understood.

Pinned mode makes chosen update entries inactive exactly:

```toml
[gauge]
mode = "pin"
pin_event_ids = [101, 202]
pin_fields = "xyz"
pin_reference_catalog = "trusted_catalog.txt"
```

If no trusted catalog is supplied, those events remain at their input-catalog
locations under `initialization.mode = "catalog"`, or at the overridden common
seed under a common initialization mode. A trusted catalog must contain the
same serial IDs for every pin,
but it may contain only those pinned rows. IDs—not row numbers or row order—are
matched. For example, a three-row trusted file containing IDs 3, 50, and 67 is
sufficient for `pin_event_ids = [3, 50, 67]`.

A hard pin cannot be weakened by Huber weighting, but pinning is not a
catalog-translation operator. Moving only one or a few events to a distant
absolute location and pinning them creates large residuals on their links to
the swarm. Robust weighting can then give those links very little leverage,
and the Stage-2 graph constructed after prelocation can put the displaced pins
in another component. Damping does not repair that loss of connection.

Accordingly, do not use a displaced single-event or few-event pin merely to
shift the absolute catalog. Prefer `common_centroid` or `common_manual` when
the goal is to locate without individual seed hypocenters. If a displaced-pin
sensitivity test is nevertheless attempted, use `pin_fields = "xyz"`, inspect
the pin component in `catalog_dd_graphmeta.csv`, and set both
`prelocation.huber_k` and `relocation.huber_k` to a very large value so the
connecting residuals are not treated as outliers. `gauge.constraint_weight`
has no effect on fields removed by hard pinning, and pinning `t0` is normally
unnecessary.

`linear_solver = "direct"` exists for small tests and is limited to 4000 free
parameters. Use PCG for production. Step clipping is disabled by default;
`max_event_step_m` and `max_origin_step_s` are emergency safeguards rather than
a convergence strategy.

## 8. Recommended two-pass workflow

A second run is simply another GraphSplit invocation. It does not need special
code.

First pass:

```toml
[run]
prelocation = true

[graph]
neighbors = 60
maximum_degree = 120
maximum_distance_km = 8.0

[gauge]
zero_mean = "xyzt0_scaled"
constraint_weight = 1.0
```

Second pass:

```toml
[io]
catalog_file = "../pass1/catalog_dd.txt"
output_dir = "pass2"

[run]
prelocation = false

[graph]
neighbors = 30
maximum_degree = 60
maximum_distance_km = 5.0

[gauge]
zero_mean = "origin_time"
constraint_weight = 10.0
```

This uses the first relocation as the seed, avoids repeating Stage 1, and lets
the local graph contract around improved hypocenters. The numerical values are
a starting schedule, not universal defaults.

## 9. Uncertainty estimates

Uncertainty output is deliberately separate from the four location catalogs.
This keeps `catalog_dd.txt` compatible with DDSync/HypoDD-style workflows and
prevents dozens of diagnostic columns from becoming part of the catalog format.

### Station-phase block bootstrap

The maintained empirical estimate resamples complete station-phase theta files.
For `G` retained groups, a replicate draws `G` groups with replacement. A group
selected twice contributes twice to the robust objective; a group not selected
does not contribute. P and S at the same station remain separate by default.
This is preferable for sparse networks, where dropping a whole station may
reduce the effective geometry to only two stations. Whole-station resampling is
available as a more conservative sensitivity test.

Each replicate starts from the final all-data solution, holds the nominal event
graph and support filtering fixed, and reruns the nonlinear Stage-2 solve. The
fixed graph makes the result interpretable as sensitivity to station-phase
support instead of a mixture of support and graph-selection changes. Starting
at the all-data solution is only an optimizer warm start; it does not constrain
the replicate to remain there.

The wide sample files contain `EventID`, `all_data`, and one column per
replicate. Longitude, latitude, depth, and relative `t0` are separate files.
These samples—not an ellipse—are the primary result, because weak networks can
produce asymmetric, elongated, or multimodal clouds. `NaN` means that the event
lost all resampled support or that the replicate did not converge. Optional
per-replicate catalogs contain only events active in that realization.

`booterrxyz.txt` is a convenience summary in the local east/north/depth frame.
It reports valid-sample counts, central offsets and bounds, and all six unique
covariance terms. Bounds may be central percentiles or mean ± a selected number
of standard deviations. Always inspect sample clouds for important events.

### Regularized linearized uncertainty

The scalable alternative freezes the final robust linearization and estimates
the spatial diagonal blocks of

```text
(J' W J + damping + gauge constraints)^-1
```

with randomized matrix-free solves. `linerrxyz.txt` contains x/y/z standard
deviations and covariance terms for filtered IDs. This is much cheaper than
many nonlinear relocations and remains practical for very large catalogs, but
it is conditional on the chosen graph, weights, damping, velocity model, and
linearization. The word “regularized” is important: damping and pins can make
the formal spread small without eliminating real model error.

Neither estimate includes velocity-model uncertainty. The GraphSplit-only
bootstrap also conditions on DDSync's theta estimates. A complete end-to-end
bootstrap would start from the original differential-time observations and
rerun both DDSync and GraphSplit.

## 10. Experimental bias workflow

`experimental.bias.enabled = true` writes `theta_bias_report.csv`. For each
station–phase group with enough observations, it robustly fits a no-intercept
linear spatial residual field and reports the before/after RMS. A row is marked
applicable only when it exceeds `minimum_improvement_fraction`.

Nothing is applied during the scan run. To test a reviewed report, set
`experimental.bias.apply_model_file` in a later run. This explicit two-step
design prevents an exploratory correction from silently changing the primary
solution.

## 11. Reading the diagnostics

Check `solver_history.csv` for decreasing robust RMS and spatial steps. A PCG
warning means the outer step used the best available inner iterate, not that the
entire run necessarily failed. Repeated warnings usually indicate a weak graph,
an overly tight inner tolerance, a poor seed, or damping that is too small.

`catalog_dd_graphmeta.csv` distinguishes geometric graph degree from retained
DD pair degree and raw observation count. Events missing from
`catalog_dd_filt.txt` remain in the full catalog unchanged or weakly constrained
but did not participate in the retained DD system.

## 12. Reproducibility checklist

Archive the TOML, input catalog, station and velocity files, DDSync output,
native table header/model file, `run_summary.toml`, graph metadata, and solver
history. The serial ID is the join key across all of them.
