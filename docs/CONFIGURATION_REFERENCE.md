# GraphSplit configuration reference

This page explains every setting accepted by `graphsplit.toml`. It is written
for users who want to understand what will change their earthquake locations,
without assuming knowledge of graph theory or inverse-problem terminology.

GraphSplit starts with its built-in defaults and replaces only the values
listed in your TOML file. Paths are resolved relative to the TOML file—not the
terminal's current directory. Unknown names produce a warning, which helps
catch spelling mistakes. The copy-ready minimal file is
`config/graphsplit_template.toml`; the file containing every setting and its
default is `config/graphsplit_complete.toml`.

## 1. Which settings should I care about first?

Most users should tune settings in this order.

### Highest scientific impact

1. **Travel-time model and geometry:** `travel_time.velocity_model_file`,
   `travel_time.build_geometry`, and the lookup spacing determine how travel
   times and their spatial gradients are represented.
2. **Starting locations and prelocation:** `io.catalog_file`,
   `initialization.mode`, and `run.prelocation` determine the seed supplied to
   the final relocation. Common initialization allows location without
   individual catalog hypocenters.
3. **Spatial damping:** `prelocation.damping_lambda` and
   `relocation.damping_lambda` suppress poorly constrained motion. If the
   solution moves coherently or unrealistically away from a plausible seed,
   increase damping before changing many other settings.
4. **Which event pairs are used:** the `[graph]` settings control which nearby
   earthquakes may be compared. They influence resolution and computing cost.
5. **Which measurements are trusted:** the `[observations]` settings control
   DDSync degree, uncertainty, pair-support, and event-support filtering.
6. **Reference frame or trusted events:** `[gauge]` determines whether the
   solution is free to translate weakly, has a zero-mean update, or is tied to
   exact trusted locations.
7. **External depth information:** `[constraints.depth_bound]` prevents
   physically impossible shallow solutions; `[constraints.fixed_depth]`
   removes depth from the solve when it is known independently.

### Usually changed only after inspecting diagnostics

- `step_damping`, Huber settings, step limits, and stopping thresholds;
- graph connectivity repair and graph augmentation;
- PCG tolerance, iteration limit, and preconditioner.

### Primarily performance or bookkeeping

- travel-time table memory/threading controls;
- output switches and console verbosity.

## 2. Common symptoms and the first setting to inspect

| Symptom | First checks | Typical response |
| --- | --- | --- |
| Catalog drifts coherently or moves implausibly far | Seed quality, gauge, `solver_history.csv` step sizes | Increase the relevant `damping_lambda` by a factor of 2–10; consider `gauge.zero_mean = "xyzt0_scaled"` for a weak first pass |
| Only a small fraction of events relocate | `catalog_dd_graphmeta.csv`, DDSync degree, pair support | Relax `minimum_theta_degree` or `minimum_observations_per_pair` cautiously; enlarge the graph if pairs are missing |
| Many graph components or zero-degree events | Graph degree/radius, seed locations | Increase `graph.neighbors`, relax `graph.maximum_distance_km`, or set `mutual = false`; use connectivity repair only after checking the cause |
| No reliable single-event seed locations | `initialization.mode`, Stage-1 coverage | Use `common_centroid` or a reasonable `common_manual` seed and keep `run.prelocation = true` |
| A DD-only second pass starts with a much larger residual than the preceding solution | `io.restart_shift_file` and the catalog/shift pairing | Use the `catalog_dd.txt` and `catalog_dd_dxdydzt0.txt` produced by the same pass; the shift file restores the cumulative relative origin times |
| PCG repeatedly reaches its iteration limit | Graph support, damping, preconditioner | Increase damping, retain `block_jacobi`, or increase `inner_max_iterations`; do not begin by making `inner_tolerance` tighter |
| Locations pile up at a table boundary | Lookup depth/range and clamping | Rebuild a larger table; clamping prevents a crash but does not make boundary locations reliable |
| Sparse-station events mirror above the surface | `constraints.depth_bound`, Stage-1 depth geometry | Add an explicit physical minimum depth; enable the reflected Stage-1 restart when the deeper mirror branch is expected |
| Explosions have a known surface or bench depth | `constraints.fixed_depth` | Fix z exactly for all events or selected serial IDs while continuing to solve x, y, and t0 |
| Run is too slow or table is too large | Lookup spacings and station-depth slices | Coarsen the table carefully, especially `station_depth_step_m`; compare locations before accepting reduced resolution |

The numerical values below are defaults, not universal recommendations.

## 3. Files and directories: `[io]`

The five ordinary paths and the optional restart path may be relative to
`graphsplit.toml` or absolute.

| Setting | Default | Meaning |
| --- | --- | --- |
| `io.catalog_file` | `"catalog.txt"` | Earthquake catalog. It supplies the persistent serial event ID and preserved columns; by default it also supplies individual seed latitude, longitude, and depth. `[initialization]` can override those individual seeds in memory. |
| `io.restart_shift_file` | `""` | Cumulative `catalog_dd_dxdydzt0.txt` paired with `io.catalog_file` for an advanced DD-only second pass. It is required when `run.prelocation = false`, ignored only for a table-only invocation, and must be empty for an ordinary run with Stage 1. |
| `io.stations_file` | `"stations.txt"` | Station coordinates and elevations. Station codes must exactly match the codes in theta filenames. |
| `io.theta_dir` | `"theta"` | Directory containing DDSync `theta_<STA>_<P\|S>.txt` files. |
| `io.thetastd_dir` | `"thetastd"` | Directory containing matching DDSync `std_theta_<STA>_<P\|S>.txt` uncertainty/degree files. Missing individual thetaStd files are allowed, but their observations then use the sigma floor and cannot be degree-filtered. |
| `io.output_dir` | `"graphsplit_output"` | Directory receiving catalogs and diagnostics. It is created when necessary. |

See [FILE_FORMATS.md](FILE_FORMATS.md) for column-by-column examples.

## 4. Run control: `[run]`

| Setting | Default | Meaning |
| --- | --- | --- |
| `run.prelocation` | `true` | Run Stage 1 before the sparse pair relocation. Set to `false` only for a DD-only pass whose input is an earlier `catalog_dd.txt` and whose `io.restart_shift_file` is the matching shift sidecar. Stage 2 always runs unless this is a table-only invocation. |
| `run.overwrite` | `true` | Allow primary catalog files in `output_dir` to be replaced. `false` stops before overwriting an existing primary result. It does not delete unrelated old files from the directory. |
| `run.build_travel_times_only` | `false` | Prepare/validate the travel-time model and return without loading theta data or relocating. The command-line `--build-tt-only` option provides the same practical workflow. |

## 5. Catalog columns: `[catalog]`

Column numbers are one-based. A negative number counts
backward from the end: `-1` is the last column and `-2` is the next-to-last.

| Setting | Default | Meaning |
| --- | --- | --- |
| `catalog.latitude_column` | `7` | Latitude column in decimal degrees. |
| `catalog.longitude_column` | `8` | Longitude column in decimal degrees. |
| `catalog.depth_column` | `9` | Event depth column in kilometres, interpreted according to `coordinates.event_vertical`. |
| `catalog.event_id_column` | `-1` | Persistent positive integer serial `EventID`. IDs must be unique, but need not be consecutive or equal to row numbers. |
| `catalog.origin_time_columns` | `[1, 2, 3, 4, 5, 6]` | Year, month, day, hour, minute, and second columns. GraphSplit adds the solved origin-time correction with full calendar rollover. Use `[]` when the catalog has no usable absolute origin times; the fields then remain unchanged and `dt0_s` is still available in the shift sidecar. |
| `catalog.origin_time_decimals` | `6` | Decimal places written in the corrected seconds field. This controls text precision only, not the internal solve. Allowed range: 0–12. |

The six origin-time columns, latitude, longitude, depth, and ID columns must be
distinct. GraphSplit replaces the configured location and origin-time fields
and preserves all other columns. If an active event has invalid calendar fields,
GraphSplit leaves those fields unchanged, emits a warning, and retains its
correction in `catalog_*_dxdydzt0.txt`.

## 6. Starting hypocenters: `[initialization]`

This section controls only the in-memory starting locations. It does not alter
the input file and does not add location constraints to the inversion.

| Setting | Default | Meaning |
| --- | --- | --- |
| `initialization.mode` | `"catalog"` | `catalog` uses every input latitude, longitude, and depth as its event's seed. `common_centroid` replaces all seeds by the mean input x, y, and depth. `common_manual` replaces all seeds by the location in the next three settings. Both common modes require Stage 1. |
| `initialization.latitude` | `0.0` | Common starting latitude in decimal degrees, used only by `common_manual`. It is an initial guess, not a fixed solution. |
| `initialization.longitude` | `0.0` | Common starting longitude in decimal degrees, used only by `common_manual`. |
| `initialization.depth_km` | `0.0` | Common starting depth in kilometres, used only by `common_manual` and interpreted using `coordinates.event_vertical`. |

The common modes are useful when DDSync has been run but a conventional
single-event locator has not produced a trustworthy input catalog. Stage 1
separates the initially coincident events using station-phase theta data;
Stage 2 then constructs its event graph from the Stage-1 locations. GraphSplit
rejects a common mode with `run.prelocation = false`, because nearest neighbors
among coincident points would be arbitrary.

`common_centroid` still uses the catalog only to select a reasonable shared
starting point. `common_manual` avoids using its location geometry altogether:

```toml
[initialization]
mode = "common_manual"
latitude = 64.0200
longitude = -21.2100
depth_km = 3.0

[run]
prelocation = true
```

If a common-start run struggles initially, inspect Stage-1 residuals and
coverage. A larger `prelocation.huber_k` can be a useful sensitivity test when
the first residuals are many sigma, but it should not be selected solely to
force a desired location.

## 7. Internal coordinates: `[coordinates]`

Most datasets should keep every value in this section at its default. These
settings describe how geographic inputs are converted into internal metres.

| Setting | Default | Meaning |
| --- | --- | --- |
| `coordinates.reference` | `"stations_mean"` | Origin of the local x/y coordinate frame. `stations_mean` uses the mean station latitude and circular-mean longitude; `catalog_mean` uses events; `manual` uses the next two settings. |
| `coordinates.reference_latitude` | `0.0` | Manual origin latitude in degrees; used only when `reference = "manual"`. |
| `coordinates.reference_longitude` | `0.0` | Manual origin longitude in degrees; used only when `reference = "manual"`. |
| `coordinates.earth_radius_m` | `6371000.0` | Radius used for geographic-to-local x/y conversion. This is separate from the optional travel-time Earth radius below. |
| `coordinates.station_vertical` | `"depth_from_elevation"` | `depth_from_elevation` converts elevation to depth using `z = -elevation_m`, so a station above sea level has negative depth. `elevation` uses the input number directly as internal z and is only for a nonstandard station file already using the desired sign. |
| `coordinates.event_vertical` | `"positive_depth"` | `positive_depth` interprets catalog depth as positive downward. `negative_depth` reverses its sign. `positive_depth_plus_z0` adds `event_z0_m` internally. |
| `coordinates.event_z0_m` | `0.0` | Vertical offset in metres, used only by `positive_depth_plus_z0`; it is removed again when writing output. |

## 8. Travel-time source and lifecycle: `[travel_time]`

| Setting | Default | Meaning |
| --- | --- | --- |
| `travel_time.type` | `"lookup"` | `lookup` uses the native layered-model table. `constant` bypasses table building and is intended mainly for synthetic tests. Experimental `3d` uses `[grid3d]`; see the [3D guide](THREE_DIMENSIONAL_MODELS.md). |
| `travel_time.table_file` | `"lookuptable/graphsplit_tt.gstt"` | Native memory-mapped table path. The parent directory is created automatically. Only the native `.gstt` format is supported. |
| `travel_time.velocity_model_file` | `"vm.txt"` | Layered P/S velocity model used to build or validate a lookup table. Its SHA-256 digest is stored in the table. |
| `travel_time.geometry` | `"auto"` | Geometry required when opening a table: `auto`, `cartesian`, or `radial`. `auto` reads an existing table's geometry. If a table must be built, `auto` delegates to `build_geometry`. An explicit geometry rejects/rebuilds a table of the other type. |
| `travel_time.build_geometry` | `"cartesian"` | Builder used when `geometry = "auto"` and no compatible table exists. `cartesian` is flat local range; `radial` uses spherical central angle. |
| `travel_time.auto_build` | `true` | Build the table if it is missing. When `false`, a missing table is an error. |
| `travel_time.auto_rebuild` | `true` | Rebuild an existing table when its model hash, geometry, Earth radius, spacing, or catalog/initialized-seed coverage is incompatible. When `false`, incompatibility is an error. |
| `travel_time.clamp_to_grid` | `true` | If a nonlinear update moves an event just outside the table, evaluate at the nearest boundary instead of returning a non-finite time. This avoids an abrupt failure, but repeated boundary locations mean the table should be enlarged and rebuilt. |
| `travel_time.vp_ms` | `6000.0` | P velocity in m/s, used only when `type = "constant"`. |
| `travel_time.vs_ms` | `3464.0` | S velocity in m/s, used only when `type = "constant"`. |
| `travel_time.earth_radius_m` | `0.0` | Radius used by radial travel times. Zero selects the radius inferred from a four-column model, or falls back to `coordinates.earth_radius_m`. A positive value overrides both. |

Example for a spherical table:

```toml
[travel_time]
geometry = "auto"
build_geometry = "radial"
```

Once created, the `.gstt` header remembers that it is radial.

## 9. Native lookup-table grid: `[lookup]`

Smaller spacing improves table resolution but increases build time and disk
size. Approximate stored size is
`2 phases × 4 bytes × Nrange × Ndepth × Nstation_depth`.

| Setting | Default | Meaning |
| --- | --- | --- |
| `lookup.horizontal_step_m` | `125.0` | Horizontal range spacing in metres. In radial mode this is a surface-arc spacing converted internally to angle. Changing it makes an existing table incompatible. |
| `lookup.depth_step_m` | `50.0` | Event-depth grid spacing in metres. It also sets the finest available station-depth spacing. |
| `lookup.station_depth_step_m` | `100.0` | Requested maximum spacing between station-depth slices. It must be at least `depth_step_m`; endpoint alignment may make actual spacing finer. More slices improve elevation interpolation but cost proportionally more time and disk. |
| `lookup.maximum_distance_km` | `0.0` | Table range. Zero derives it from the largest event-to-station distance across both the input catalog and initialized seed, plus a margin. A positive explicit value must cover both. This is unrelated to `graph.maximum_distance_km`. |
| `lookup.minimum_depth_km` | `0.0` | Shallow table limit. Zero automatically includes the shallowest model node, station, and seed event. A nonzero value is explicit; negative values are useful for stations above sea level. |
| `lookup.maximum_depth_km` | `0.0` | Deep table limit. Zero chooses an automatic limit; a positive value is explicit and must contain all seed events and stations. |
| `lookup.distance_margin_km` | `1.0` | Extra range used only when `maximum_distance_km = 0`. The builder uses at least two horizontal grid cells even if this margin is smaller. |
| `lookup.depth_margin_km` | `50.0` | Extra depth below the deepest seed event used by automatic `maximum_depth_km`. The automatic depth also considers range and velocity-model extent. Reduce cautiously if table size is excessive. |
| `lookup.maximum_sweeps` | `32` | Maximum fast-sweeping cycles for each phase and station-depth slice. Increase only if the builder reports non-converged slices. |
| `lookup.sweep_tolerance_s` | `1.0e-7` | Maximum travel-time change, in seconds, used to declare a slice converged. Smaller is stricter and may require more sweeps. |
| `lookup.threaded` | `false` | Build independent station-depth slices in parallel. Start Julia with multiple threads, for example `julia -t auto ...`; otherwise this switch has no speed effect. Parallel slices can increase peak memory and make progress messages less orderly. |
| `lookup.maximum_table_gib` | `8.0` | Safety limit on the estimated P+S cube size in GiB. |
| `lookup.allow_large_table` | `false` | Permit a table above `maximum_table_gib`. Set only after checking available disk, build time, and memory pressure. |

## 10. Stage 1 and Stage 2 solver settings

`[prelocation]` controls the optional theta-reference Stage 1.
`[relocation]` controls the final sparse event-pair Stage 2. The same setting
names are accepted in both sections; defaults that differ are shown separately.

| Setting | Prelocation default | Relocation default | Meaning |
| --- | ---: | ---: | --- |
| `prelocation.max_outer_iterations` / `relocation.max_outer_iterations` | 20 | 60 | Maximum nonlinear location iterations. Reaching this limit is not automatically a failure; inspect residual and step histories. |
| `prelocation.min_outer_iterations` / `relocation.min_outer_iterations` | 3 | 3 | Minimum iterations before a stopping rule is allowed. |
| `prelocation.linear_solver` / `relocation.linear_solver` | `"pcg"` | `"pcg"` | `pcg` is the matrix-free production solver. `direct` explicitly forms a dense operator and is limited to 4,000 free parameters; use it only for small tests. |
| `prelocation.preconditioner` / `relocation.preconditioner` | `"block_jacobi"` | `"block_jacobi"` | PCG acceleration: `block_jacobi` uses one coupled x/y/z/time 4×4 block per event and is recommended; `diagonal` is cheaper but weaker; `none` is primarily diagnostic. |
| `prelocation.inner_tolerance` / `relocation.inner_tolerance` | `1.0e-5` | `1.0e-4` | PCG relative residual tolerance for each linearized step. Smaller solves each step more accurately but may not improve the nonlinear result. |
| `prelocation.inner_max_iterations` / `relocation.inner_max_iterations` | 250 | 300 | Maximum PCG iterations per nonlinear iteration. A warning means the best available step was used. Repeated warnings deserve investigation. |
| `prelocation.huber_k` / `relocation.huber_k` | 1.345 | 1.345 | Huber threshold in units of each observation's sigma. Residuals within `k × sigma` retain full weight; larger residuals are progressively downweighted. Smaller values reject discordant data more aggressively. |
| `prelocation.min_sigma_s` / `relocation.min_sigma_s` | 0.002 | 0.002 | Lower uncertainty floor in seconds. It prevents exceptionally small or missing theta uncertainties from receiving unbounded weight. With no thetaStd, a pair receives `sqrt(2) × min_sigma_s`. |
| `prelocation.step_damping` / `relocation.step_damping` | 0.75 | 0.75 | Fraction of the computed update applied each nonlinear iteration. `1.0` takes the full step; smaller values move more cautiously. This slows all updates and is different from regularization below. |
| `prelocation.damping_lambda` / `relocation.damping_lambda` | `1.0e-3` | `5.0e-4` | Diagonal regularization added to the linearized normal equations. Larger values keep weakly resolved locations closer to their current positions and usually reduce drift, but too much damping preserves seed bias. Tune by factors of 2–10 and compare solutions. Zero selects only the internal numerical floor of `1.0e-8`, rather than mathematically exact zero damping. Because GraphSplit stores metres and seconds together internally, this scalar acts much more strongly on spatial motion than on origin-time corrections. |
| `prelocation.stop_step_rms_m` / `relocation.stop_step_rms_m` | 0.10 | 0.10 | Spatial RMS step threshold in metres. Small-step convergence requires this and the time threshold below simultaneously. |
| `prelocation.stop_step_rms_s` / `relocation.stop_step_rms_s` | `1.0e-4` | `1.0e-4` | Origin-time RMS step threshold in seconds. |
| `prelocation.stop_rms_improvement_s` / `relocation.stop_rms_improvement_s` | `1.0e-5` | `1.0e-5` | A robust-RMS change whose absolute value is at or below this many seconds counts as a stalled iteration. |
| `prelocation.stop_stall_iterations` / `relocation.stop_stall_iterations` | 5 | 5 | Consecutive stalled iterations required to stop. Zero disables stall stopping. |
| `prelocation.max_event_step_m` / `relocation.max_event_step_m` | 0.0 | 0.0 | Maximum 3-D movement of one event in one iteration. Zero disables clipping. This is an emergency safeguard, not a substitute for damping or a better seed. |
| `prelocation.max_origin_step_s` / `relocation.max_origin_step_s` | 0.0 | 0.0 | Maximum absolute origin-time change per event per iteration. Zero disables clipping. |
| `prelocation.verbose` / `relocation.verbose` | `true` | `true` | Print one convergence line per nonlinear iteration. It does not change results. |

A conservative response to obvious drift is, for example:

```toml
[prelocation]
damping_lambda = 0.005

[relocation]
damping_lambda = 0.0025
```

These values are five times the defaults. They are a diagnostic trial, not a
universal setting: compare movement, residual, and benchmark/independent
location evidence before accepting the more strongly damped result.

## 11. Observation selection: `[observations]`

Here, one “observation” means one station–phase contribution, such as P at
station ABC for a particular event pair. It does not mean one earthquake.

| Setting | Default | Meaning |
| --- | --- | --- |
| `observations.minimum_theta_degree` | 6 | Minimum DDSync graph degree recorded for an event in thetaStd. It filters poorly synchronized theta entries. If a thetaStd file is absent or has no degree column, this filter cannot be applied to that station–phase group. Zero disables it. |
| `observations.maximum_sigma_s` | 0.0 | Maximum combined theta uncertainty in seconds. For a pair it is `hypot(sigma_i, sigma_j)`. Zero disables the upper limit. |
| `observations.minimum_observations_per_pair` | 6 | Minimum number of surviving station–phase contributions required to keep an event pair in Stage 2. For example, four P stations plus two S stations give six contributions. |
| `observations.minimum_observations_per_event` | 0 | Minimum retained Stage-2 contributions touching each event. When positive, low-support events and their observations are removed iteratively because removing one event can weaken another. Zero disables this peel. |
| `observations.minimum_star_observations_per_event` | 0 | Minimum Stage-1 theta-reference contributions for the non-reference event. Observations belonging to events below the threshold are removed; zero disables it. This is Stage 1 only. |
| `observations.minimum_component_size` | 0 | Reject Stage-2 observations for events lying in geometric graph components smaller than this event count. Zero disables it. This uses the graph before measurement-support filtering. |

Start with the defaults. Relaxing filters increases coverage but may introduce
poor DDSync constraints; tightening them improves selectivity but can leave
events unchanged or split the usable catalog.

## 12. Reference frame and pinned events: `[gauge]`

Relative times alone do not define every possible common shift of the model.
The gauge tells GraphSplit how to choose a stable reference frame. You do not
need linear algebra to use it: keep the default unless you see coherent drift
or have trusted event locations.

| Setting | Default | Meaning |
| --- | --- | --- |
| `gauge.mode` | `"zero_mean"` | `zero_mean` adds a weak/common-shift reference described below. `pin` holds selected event fields exactly fixed. |
| `gauge.zero_mean` | `"origin_time"` | Used in zero-mean mode. `origin_time` removes only the common relative-time shift. `xyzt0_scaled` also discourages a common spatial translation after converting metres to seconds with `reference_velocity_ms`; it is useful for a weakly anchored first pass. `xyzt0` is unscaled and should normally be avoided. |
| `gauge.constraint_weight` | 10.0 | Strength of the zero-mean constraint relative to observations. It is not the same as damping. Larger values enforce the selected mean update more strongly. The origin-time mean is also recentered exactly after each iteration. It has no effect on fields removed by hard pinning. |
| `gauge.reference_velocity_ms` | 5000.0 | Metres-to-seconds scaling for `zero_mean = "xyzt0_scaled"`. It has no effect in other gauge choices. |
| `gauge.pin_event_ids` | `[]` | Persistent serial `EventID` values to hold fixed in pin mode, for example `[3, 50, 67]`. These are IDs from the configured catalog ID column—not row numbers. At least one is required when `mode = "pin"`. |
| `gauge.pin_fields` | `"xyz"` | Fields fixed for each listed event: any combination of `x`, `y`, `z`, and `t0`, such as `"xy"`, `"z"`, or `"xyzt0"`. `xyz` fixes the location but still solves relative origin time. |
| `gauge.pin_reference_catalog` | `""` | Optional catalog supplying trusted latitude/longitude/depth for pinned IDs. Empty means pin at the in-memory seed: normally `io.catalog_file`, or the common seed selected by `[initialization]`. The trusted file may contain only the pinned events and may use any row order, but it must use the same column layout configured under `[catalog]` and contain every listed serial ID. Trusted catalog date/time fields are not imported; a pinned `t0` is held at zero in an ordinary run or at its restored cumulative value in a DD-only pass. |

### Trusted catalog example

Suppose the working catalog has 1,000 events and trusted locations exist only
for serial IDs 3, 50, and 67. This is sufficient:

```toml
[gauge]
mode = "pin"
pin_event_ids = [3, 50, 67]
pin_fields = "xyz"
pin_reference_catalog = "catalog_trusted.txt"
```

`catalog_trusted.txt` may contain exactly three standard catalog rows:

```text
2025 01 01 00 00 03.000 35.410000 -117.910000 4.2000 1.0 3
2025 01 01 00 00 50.000 35.420000 -117.900000 5.1000 1.0 50
2025 01 01 00 01 07.000 35.430000 -117.890000 6.0000 1.0 67
```

Only columns 7, 8, 9, and the final serial ID are used here. The rows do not
need to occupy positions 3, 50, and 67, and the other 997 events do not need to
appear. “Same ID” means that the final values `3`, `50`, and `67` match the
working catalog.

### Why a displaced pin does not translate the catalog

Hard pins are exact: Huber weighting cannot move a pinned field. However,
moving only one event or a small subset to a distant absolute position and
pinning it is not a reliable way to shift the remaining swarm. The links from
the displaced events acquire large residuals, robust weighting can reduce
their leverage severely, and the event graph built after Stage 1 may place the
pins in a separate component. Increasing damping cannot transmit an absolute
shift across downweighted or missing connections.

Use `common_centroid` or `common_manual` when the objective is a full location
without individual input hypocenters. Reserve pins for trusted constraints in
a consistently initialized, connected catalog. If a displaced-pin sensitivity
test is attempted anyway, use `pin_fields = "xyz"`, set both
`prelocation.huber_k` and `relocation.huber_k` to a very large value, and verify
in `catalog_dd_graphmeta.csv` that the pins belong to the main retained
component. A very large Huber threshold disables most robust downweighting and
is appropriate here only as a diagnostic.

## 13. Physical and fixed depths: `[constraints]`

These settings are independent of `[gauge]`. A gauge selects the reference
frame of the relative problem; a depth constraint supplies external physical
information about particular z parameters. Both constraints apply in Stage 1,
Stage 2, and every nonlinear bootstrap replicate.

All configured depths use the catalog convention and kilometres. With the
default `coordinates.event_vertical = "positive_depth"`, `-0.8` km means an
elevation of 800 m above the datum.

### No-cross physical depth: `[constraints.depth_bound]`

| Setting | Default | Meaning |
| --- | --- | --- |
| `constraints.depth_bound.enabled` | `false` | Enable a physical minimum catalog depth. This is a model constraint, unlike `lookup.minimum_depth_km`, which only sizes the travel-time table. |
| `constraints.depth_bound.minimum_depth_km` | 0.0 | Shallowest permitted catalog depth in kilometres. In the default positive-down convention, every selected event must satisfy `depth >= minimum_depth_km`. |
| `constraints.depth_bound.scope` | `"all"` | `all` applies to every event and requires an empty ID list. `event_ids` applies only to the explicit serial IDs below and requires a nonempty list. |
| `constraints.depth_bound.event_ids` | `[]` | Persistent serial IDs used only with `scope = "event_ids"`. Missing, duplicate, or nonpositive IDs are errors. |
| `constraints.depth_bound.reflected_prelocation_restart` | `false` | When true, run an unconstrained Stage-1 pilot, reflect forbidden pilot depths, rerun bounded Stage 1, and build the Stage-2 graph only from that bounded solution. It requires `run.prelocation = true`. When false, the ordinary Stage-1 solve is bounded directly and all starting depths must already be feasible. |
| `constraints.depth_bound.pilot_shallow_margin_km` | 10.0 | Extra forbidden-side depth range reserved in an automatically built lookup table for the unconstrained pilot. It is used only by the reflected restart. Clamping is disabled during the pilot, so an explicitly sized or reused table that is still too small produces an error instead of a zero depth gradient. Increase this value and rebuild if the pilot leaves the table. |
| `constraints.depth_bound.mirror_plane` | `"stations_median"` | Plane for the explicit restart: `stations_median` uses the median internal station depth; `manual` uses `mirror_depth_km`. This has no effect when the reflected restart is off. |
| `constraints.depth_bound.mirror_depth_km` | 0.0 | Manual mirror-plane depth in catalog kilometres, used only with `mirror_plane = "manual"`. Reflected candidates must land inside the physical bound or GraphSplit stops with an error. |

The active-set solver first computes a trial update. If a selected event would
cross the physical bound, its z update is placed exactly on the boundary and
the coupled x, y, and t0 update is re-solved. At the next nonlinear iteration,
z is released automatically if its trial direction points into the admissible
interior. This avoids using table clamping as an accidental geological
constraint.

An event may still finish exactly at the bound. That means the constrained
objective prefers the forbidden side; it does not mean the boundary depth was
measured. Check `depth_constraint_status.csv` and, for bootstrap runs,
`bootstrap_depth_bound_status.txt`.

For sparse stations at similar elevation, the reflected restart searches the
deeper branch before the event graph is built:

```toml
[constraints.depth_bound]
enabled = true
minimum_depth_km = -0.8
scope = "all"
event_ids = []
reflected_prelocation_restart = true
pilot_shallow_margin_km = 10.0
mirror_plane = "stations_median"
mirror_depth_km = -0.8 # ignored for stations_median
```

The lookup-table compatibility check includes the physical boundary, the
configured pilot margin, and the reflection of the seed depth range. During
the pilot, table clamping is disabled. If the nonlinear search still leaves
that range, GraphSplit stops; increase `pilot_shallow_margin_km` and rebuild
the table rather than accepting zero out-of-grid depth gradients.

### Exact known depths: `[constraints.fixed_depth]`

| Setting | Default | Meaning |
| --- | --- | --- |
| `constraints.fixed_depth.enabled` | `false` | Freeze selected z parameters exactly while continuing to solve x, y, and relative origin time. |
| `constraints.fixed_depth.scope` | `"all"` | `all` fixes every event and requires an empty ID list. `event_ids` fixes only the listed serial IDs and requires at least one. |
| `constraints.fixed_depth.event_ids` | `[]` | Persistent serial IDs used with `scope = "event_ids"`. These are never interpreted as row numbers. |
| `constraints.fixed_depth.depth_km` | 0.0 | One fixed catalog depth in kilometres for every selected event. Ignored when `reference_catalog` is nonempty. |
| `constraints.fixed_depth.reference_catalog` | `""` | Optional same-format catalog supplying an individual depth for every selected ID. It may contain only the selected rows and may use any order. Latitude, longitude, time, and auxiliary columns are ignored. |

For explosions known to occur on a bench at 800 m elevation:

```toml
[constraints.fixed_depth]
enabled = true
scope = "all"
event_ids = []
depth_km = -0.8
reference_catalog = ""
```

Do not constrain the same z parameter through both fixed depth and a gauge pin;
GraphSplit rejects that ambiguity. Remove `z` from `gauge.pin_fields` when x/y
pinning and an independent fixed depth are both intended. A fixed-depth zero in
`linerrxyz.txt` is conditional on the imposed equality, not an estimated zero
geological error.

## 14. Choosing nearby event pairs: `[graph]`

GraphSplit does not compare every event with every other event. For each event
it finds a limited number of nearby candidates, producing a sparse list of
pairs. This is the main reason very large problems remain tractable.

| Setting | Default | Meaning |
| --- | --- | --- |
| `graph.neighbors` | 20 | Number of nearest candidate events requested for each event before mutual/radius/degree filtering. Larger values give more possible pairs and cost more memory and time. |
| `graph.maximum_degree` | 30 | Maximum number of base-graph connections retained for one event. A value at or below zero removes this cap. Connectivity repair and augmentation occur later and may exceed it. |
| `graph.metric` | `"xyz_scaled"` | Coordinates used to decide what is “near”: `xy` uses epicentral position only; `xyz` uses full 3-D metres; `xyz_scaled` multiplies depth difference by `depth_scale` before calculating distance. |
| `graph.depth_scale` | 0.5 | Depth multiplier used only by `xyz_scaled`. At 0.5, a 2 km depth difference counts like 1 km horizontally when ranking neighbors. Smaller values make horizontal proximity more important. |
| `graph.mutual` | `true` | With `true`, keep a candidate only if each event selects the other among its nearest neighbors. This prevents dense hubs but can isolate events. `false` keeps a pair if either event selects the other. |
| `graph.maximum_distance_km` | 0.0 | Optional maximum neighbor distance; zero disables it. Distance is measured in the selected graph metric, so depth is ignored for `xy` and scaled for `xyz_scaled`. This is unrelated to lookup-table range. |
| `graph.minimum_degree` | 0 | Attempt to restore additional candidate edges for events with fewer connections than this value, while respecting `maximum_degree`. Zero disables the restoration pass. It cannot invent edges outside the existing neighbor/radius/mutual candidate pool. |
| `graph.ensure_connected` | `false` | Try to join separate components using unused edges from the same base candidate pool. It may exceed `maximum_degree`. It cannot bridge components if mutual selection, `neighbors`, or the radius excluded every connecting pair. |
| `graph.maximum_bridge_distance_km` | 0.0 | Maximum metric distance for edges added by `ensure_connected`; zero gives no additional bridge limit beyond the base graph radius/candidate pool. It has no effect when connectivity repair is off. |

Practical interpretation of the run message:

- `components=1` means every event is connected through some chain of candidate
  pairs; it does **not** mean every pair has enough DDSync observations.
- `degree min=0` means at least one event has no geometric pair and will not
  participate in Stage 2.
- Use `catalog_dd_graphmeta.csv` to distinguish geometric connections from
  event pairs that survived measurement filters.

## 15. Optional graph augmentation: `[graph.augmentation]`

Augmentation is an advanced, normally disabled step. It looks through a broader
nearby-pair pool and adds a small number of connections intended to strengthen
weakly linked portions of the graph. It does not create new measurements: an
added pair still disappears from the inversion if too few station–phase theta
differences survive.

| Setting | Default | Meaning |
| --- | --- | --- |
| `graph.augmentation.enabled` | `false` | Enable the augmentation pass. When true, exactly one of the next two edge targets must be positive. |
| `graph.augmentation.candidate_neighbors` | 40 | Number of nearest events searched for possible added edges. This should normally exceed `graph.neighbors`; it still obeys the graph radius and `mutual` rule. |
| `graph.augmentation.add_edges` | 0 | Exact catalog-wide target number of new undirected edges. A positive value takes the exact-count form; set `add_edges_per_event = 0.0`. |
| `graph.augmentation.add_edges_per_event` | 0.0 | Proportional target expressed as average **added degree** per event. It is a float so fractional averages are possible. The edge target is `ceil(N × value / 2)` because each undirected edge contributes degree to two events. Set `add_edges = 0` when using it. |
| `graph.augmentation.maximum_added_per_event` | 3 | Maximum augmentation edges touching any one event. This limits only newly added edges, not total graph degree. The requested catalog-wide target may not be reached if this cap or the candidate pool is restrictive. |
| `graph.augmentation.power_iterations` | 30 | Iterations used to estimate which weak graph direction should be strengthened. More iterations refine ranking but do not add more edges. The calculation is deterministic. |

Therefore, `add_edges = 0` does **not** disable `enabled = true`: it selects the
proportional alternative. If `add_edges_per_event` is also zero, GraphSplit now
stops with a clear configuration error instead of silently adding nothing.

Two valid examples for a 1,000-event catalog are:

```toml
# Add up to exactly 100 edges across the whole catalog.
[graph.augmentation]
enabled = true
candidate_neighbors = 40
add_edges = 100
add_edges_per_event = 0.0
maximum_added_per_event = 3
```

```toml
# Request average added degree 0.2: ceil(1000 × 0.2 / 2) = 100 edges.
[graph.augmentation]
enabled = true
candidate_neighbors = 40
add_edges = 0
add_edges_per_event = 0.2
maximum_added_per_event = 3
```

Leave augmentation off unless ordinary graph settings and diagnostics show a
specific connectivity weakness. It is not a general “improve locations” switch.

## 16. Output switches: `[output]`

Every completed relocation writes the two full catalogs and their
`catalog_preloc_dxdydzt0.txt` and `catalog_dd_dxdydzt0.txt` sidecars. The
sidecars use metres for x/y/z, seconds for `t0`, and contain every input event
ID; their format and cumulative semantics are described in
[FILE_FORMATS.md](FILE_FORMATS.md#relocation-shift-sidecars).

| Setting | Default | Meaning |
| --- | --- | --- |
| `output.write_filtered_catalogs` | `true` | Write `catalog_preloc_filt.txt` and `catalog_dd_filt.txt`, containing only events that participated in retained observations for the relevant stage. The full catalogs are always written. |
| `output.write_graph_metadata` | `true` | Write `catalog_dd_graphmeta.csv` with component, geometric degree, retained pair degree, observation counts, and activity flags. |
| `output.write_solver_history` | `true` | Write `solver_history.csv` with robust RMS, spatial/time step RMS, and PCG iteration count for every nonlinear iteration. |
| `output.write_run_summary` | `true` | Write `run_summary.toml` with input counts, travel-time type, graph size, observation counts, iterations, and final residuals. |

These switches affect files only, not the location solution.

## 17. Error and stability estimates: `[uncertainty]`

Uncertainty calculations run only after the ordinary all-data relocation has
finished. They write separate, serial-ID-keyed files and never add columns to
`catalog_dd.txt` or `catalog_dd_filt.txt`.

| Setting | Default | Meaning |
| --- | --- | --- |
| `uncertainty.method` | `"none"` | `none` writes no uncertainty products; `linearized` writes the inexpensive regularized inverse-Hessian estimate; `bootstrap` runs station-phase resampling; `both` does both. The location catalogs are identical for all four choices. |

### Regularized linearized estimate: `[uncertainty.linearized]`

This calculation freezes the final robust weights and travel-time derivatives,
then estimates each event's spatial block of the inverse damped Hessian. It is
conditional on the velocity model, selected event graph, final linearization,
thetaStd weights, gauge, pins, and damping. It is therefore a useful formal
resolution/error measure, not a complete statement of hypocentral uncertainty.

| Setting | Default | Meaning |
| --- | --- | --- |
| `uncertainty.linearized.probes` | 12 | Random sign probes per spatial column. Each probe costs three PCG solves—one each for x, y, and z. More probes reduce random noise in the estimated 3×3 blocks but cost proportionally more. |
| `uncertainty.linearized.seed` | 24680 | Random seed for reproducible probes. It does not affect the relocated catalog. |
| `uncertainty.linearized.inner_tolerance` | 1.0e-3 | Relative PCG tolerance for uncertainty solves. This can be looser than relocation because randomized estimation already has sampling error. |
| `uncertainty.linearized.inner_max_iterations` | 150 | Maximum PCG iterations per randomized solve. Nonconverged solves are omitted and the converged count is printed. |

The output `linerrxyz.txt` contains only filtered Stage-2 serial IDs. Columns
are `EventID`, x/y/z standard deviations in metres, and the six unique terms of
the local east/north/depth covariance matrix in square metres. Spatially pinned
coordinates have zero formal spread because they were fixed, not because their
true locations are known without error.

### Block bootstrap: `[uncertainty.bootstrap]`

The bootstrap samples complete DDSync theta groups with replacement and reruns
the nonlinear Stage-2 relocation from the all-data solution. It keeps the
nominal event graph and nominal observation-support filtering fixed. This
isolates sensitivity to station-phase sampling without allowing duplicated
blocks to pass a support threshold by being counted as independent geometry.

| Setting | Default | Meaning |
| --- | --- | --- |
| `uncertainty.bootstrap.replicates` | 100 | Number of nonlinear bootstrap relocations. Values around 100 are useful for initial assessment; stable tail percentiles may require more. Runtime is roughly proportional to this number. |
| `uncertainty.bootstrap.resampling_unit` | `"station_phase"` | `station_phase` treats one complete station–phase theta file as a block and is recommended. `station` keeps P and S together, but dropping whole stations can leave sparse networks with too little geometry and unstable replicates. |
| `uncertainty.bootstrap.seed` | 12345 | Random seed for reproducible block multiplicities. |
| `uncertainty.bootstrap.summary_method` | `"percentile"` | `percentile` reports an axis-wise central interval; `standard_deviation` reports the bootstrap mean plus/minus a configurable multiple of sample standard deviation. These summaries do not replace the saved samples. Aliases `std` and `2std` are accepted for `standard_deviation`. |
| `uncertainty.bootstrap.confidence_level` | 0.95 | Central probability used by `percentile`; 0.95 gives the 2.5th and 97.5th percentiles. Ignored by `standard_deviation`. |
| `uncertainty.bootstrap.standard_deviation_multiplier` | 2.0 | Multiplier used by `standard_deviation`; 2.0 gives mean ± 2 sample standard deviations. Ignored by `percentile`. |
| `uncertainty.bootstrap.write_samples` | `true` | Write wide longitude, latitude, depth, and internal relative-origin-time tables. Each row begins with `EventID` and the all-data value, followed by every bootstrap realization. |
| `uncertainty.bootstrap.write_catalogs` | `false` | Also write one ordinary filtered catalog per converged replicate under `bootstrap_catalogs/`. This can create many large files and is normally unnecessary because the wide tables retain every sampled location. |

`booterrxyz.txt` contains axis-wise offset summaries and the empirical local
x/y/z covariance. Its offsets are relative to the all-data `catalog_dd` result.
The four `bootstrap_samples_*.txt` files are the primary uncertainty product:
they preserve skewed, irregular, or multimodal sample clouds that an ellipse or
one number per axis would conceal. An event with no resampled observations in a
replicate is written as `NaN`, rather than being left at the all-data location
and incorrectly appearing to have zero uncertainty.
`bootstrap_block_counts.txt` records the integer multiplicity of every named
station-phase (or station) block in every replicate, so the resampling itself
can be audited independently of the random seed.

Samples are held in a temporary disk-backed array rather than resident RAM, so
bootstrap storage does not scale up the solver's memory use. Temporary binary
space is approximately `32 × filtered_events × replicates` bytes for x, y, z,
and t0, in addition to the final text files. For example, one million filtered
events and 100 replicates require about 3.2 GB of temporary disk. Setting
`write_samples = false` suppresses the final wide tables but the temporary
store is still needed to calculate percentile summaries.

The bootstrap is conditional on the theta products supplied to GraphSplit. It
does not represent velocity-model uncertainty, and it does not reproduce the
uncertainty introduced while DDSync estimated theta. A fully end-to-end study
would resample the original differential-time data, rerun DDSync, and then rerun
GraphSplit.

## 18. Experimental theta-bias diagnostic: `[experimental.bias]`

This entire section is deliberately absent from the minimal TOML. The feature
has not shown a consistent location benefit and should not be part of a default
scientific workflow.

| Setting | Default | Meaning |
| --- | --- | --- |
| `experimental.bias.enabled` | `false` | After relocation, fit and write `theta_bias_report.csv`. The scan itself does not change the current solution. |
| `experimental.bias.apply_model_file` | `""` | Path to a reviewed `theta_bias_report.csv` from an earlier run. Nonempty applies rows marked applicable to theta values before relocation. This can be used even when `enabled = false`. |
| `experimental.bias.fit_dimensions` | `"xy"` | Fit coherent residual gradients in `xy` or `xyz`. The model has no intercept and is fitted separately for each station–phase group. |
| `experimental.bias.minimum_observations` | 30 | Minimum theta-reference differences required to fit one station–phase group. |
| `experimental.bias.huber_k` | 1.345 | Huber threshold in sigma units for the robust bias fit. |
| `experimental.bias.minimum_sigma_s` | 0.002 | Uncertainty floor in seconds used only by the bias fit. |
| `experimental.bias.maximum_irls_iterations` | 6 | Maximum robust reweighting iterations for each group fit. |
| `experimental.bias.minimum_improvement_fraction` | 0.05 | Minimum fractional robust-RMS reduction required for a report row to be marked applicable. `0.05` means at least 5%. |

Always retain and report an uncorrected baseline if this experiment is used.

## 19. Complete worked configurations

### Basic run with visible damping

Start from `config/graphsplit_template.toml`. Its two damping values are shown
explicitly so a user never has to discover them in source code.

### Iterated run

```toml
[io]
catalog_file = "../pass1/catalog_dd.txt"
restart_shift_file = "../pass1/catalog_dd_dxdydzt0.txt"
output_dir = "pass2"

[run]
prelocation = false
```

The catalog and shift file must come from the same pass. You may tighten graph
radius or degree after a broad first graph has improved a poor sparse-network
seed. The restored cumulative `t0` keeps the new theta residuals consistent;
only the new time increment is added to the already-corrected calendar.

### Trusted partial pins

Use the three-line example in Section 12. The trusted file does not need to be
a duplicate of the complete working catalog.

### Reproducibility

For any published result, retain the actual TOML alongside the input catalog,
stations, velocity model, theta/thetaStd folders, `run_summary.toml`, graph
metadata, and solver history. The complete/default TOML documents software
defaults, but only the run-specific TOML records which values you changed.


## Experimental fixed 3D grids: `[grid3d]`

These settings apply only when `travel_time.type = "3d"`. Read the
[3D workflow](THREE_DIMENSIONAL_MODELS.md) before choosing a grid. The ordinary
`[lookup]` spacing settings and `travel_time.clamp_to_grid` do not control 3D
volumes. Queries outside the valid 3D domain are always rejected.

| Setting | Default | Meaning |
| --- | --- | --- |
| `grid3d.model_format` | `"nll_velocity"` | Build from velocity grids, or use `"nll_time"` for precomputed station TIME grids. |
| `grid3d.vp_file` | `"vp.mod.hdr"` | P model header; matching `.buf` is required. |
| `grid3d.vs_file` | `"vs.mod.hdr"` | S model header; required only when S is included in `phases`. |
| `grid3d.time_root` | `"time/model"` | Imported TIME filename prefix, followed by `.P.STA.time.hdr` or `.S.STA.time.hdr`. |
| `grid3d.coordinate_system` | `"header"` | Use the header's SIMPLE geographic transform. `"local"` requires a missing/NONE transform and a manual GraphSplit reference. |
| `grid3d.byte_order` | `"little"` | Binary input byte order: `little`, `big` or `native`. The header does not encode byte order. |
| `grid3d.spacing_m` | `[500.0,500.0,500.0]` | Maximum x/y/z spacing when building from velocity, metres. Imported TIME grids retain their native spacing. |
| `grid3d.bounds_km` | `[]` | Empty uses the model box. Otherwise give `[xmin,xmax,ymin,ymax,zmin,zmax]` in model coordinates. Must lie inside the velocity model; not supported for TIME imports. |
| `grid3d.model_interpolation` | `"nearest"` | Transfer input velocity by nearest node, or interpolate reciprocal velocity with `"slowness_linear"`. Only used for velocity models. |
| `grid3d.surface_file` | `""` | Optional regular `x_km y_km elevation_m` surface grid. Above-surface nodes cannot be used. Empty treats the whole box as material. |
| `grid3d.phases` | `["P","S"]` | Build/import these phases for every station in `stations.txt`. Use `["P"]` for a P-only dataset; omit unused stations to reduce storage. |
| `grid3d.accuracy_order` | `2` | FMM upwind order, 1 or 2; second order falls back where unavailable. Ignored for imported TIME grids. |
| `grid3d.cache_dir` | `"lookuptable/3d"` | One content-identified cache directory per model/settings/station combination. Interrupted or damaged caches are rebuilt when allowed. |
| `grid3d.maximum_memory_gib` | `4.0` | Refuse a larger travel-time build estimate before grid allocation. Excludes location-solver memory and the Julia runtime. |
| `grid3d.warning_memory_gib` | `1.0` | Warn at this estimate. Must not exceed `maximum_memory_gib`. |

Velocity and time volumes must share their grid and projection. Times use seconds;
velocity grids may use NonLinLoc `VELOCITY` (km/s), `VELOCITY_METERS` (m/s),
`SLOWNESS` (s/km), or `SLOW_LEN` (seconds across one x cell; equal native cell
spacing required). FLOAT and DOUBLE buffers are supported. The reader converts
units before use. Invalid values, incompatible projections, missing phases and
source-coordinate mismatches are errors.
