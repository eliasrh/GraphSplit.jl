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
2. **Starting catalog and prelocation:** `io.catalog_file` and
   `run.prelocation` determine the seed supplied to the final relocation.
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
| PCG repeatedly reaches its iteration limit | Graph support, damping, preconditioner | Increase damping, retain `block_jacobi`, or increase `inner_max_iterations`; do not begin by making `inner_tolerance` tighter |
| Locations pile up at a table boundary | Lookup depth/range and clamping | Rebuild a larger table; clamping prevents a crash but does not make boundary locations reliable |
| Run is too slow or table is too large | Lookup spacings and station-depth slices | Coarsen the table carefully, especially `station_depth_step_m`; compare locations before accepting reduced resolution |

The numerical values below are defaults, not universal recommendations.

## 3. Files and directories: `[io]`

All five paths may be relative to `graphsplit.toml` or absolute.

| Setting | Default | Meaning |
| --- | --- | --- |
| `io.catalog_file` | `"catalog.txt"` | Starting earthquake catalog. It supplies seed latitude, longitude, and depth plus the persistent serial event ID used to join every input and output. |
| `io.stations_file` | `"stations.txt"` | Station coordinates and elevations. Station codes must exactly match the codes in theta filenames. |
| `io.theta_dir` | `"theta"` | Directory containing DDSync `theta_<STA>_<P\|S>.txt` files. |
| `io.thetastd_dir` | `"thetastd"` | Directory containing matching DDSync `std_theta_<STA>_<P\|S>.txt` uncertainty/degree files. Missing individual thetaStd files are allowed, but their observations then use the sigma floor and cannot be degree-filtered. |
| `io.output_dir` | `"graphsplit_output"` | Directory receiving catalogs and diagnostics. It is created when necessary. |

See [FILE_FORMATS.md](FILE_FORMATS.md) for column-by-column examples.

## 4. Run control: `[run]`

| Setting | Default | Meaning |
| --- | --- | --- |
| `run.prelocation` | `true` | Run Stage 1 before the sparse pair relocation. Set to `false` for a second pass whose input is an earlier `catalog_dd.txt`. Stage 2 always runs unless this is a table-only invocation. |
| `run.overwrite` | `true` | Allow primary catalog files in `output_dir` to be replaced. `false` stops before overwriting an existing primary result. It does not delete unrelated old files from the directory. |
| `run.build_travel_times_only` | `false` | Prepare/validate the travel-time model and return without loading theta data or relocating. The command-line `--build-tt-only` option provides the same practical workflow. |

## 5. Catalog columns: `[catalog]`

Column numbers are one-based, as in Julia and MATLAB. A negative number counts
backward from the end: `-1` is the last column and `-2` is the next-to-last.

| Setting | Default | Meaning |
| --- | --- | --- |
| `catalog.latitude_column` | `7` | Latitude column in decimal degrees. |
| `catalog.longitude_column` | `8` | Longitude column in decimal degrees. |
| `catalog.depth_column` | `9` | Event depth column in kilometres, interpreted according to `coordinates.event_vertical`. |
| `catalog.event_id_column` | `-1` | Persistent positive integer serial `EventID`. IDs must be unique, but need not be consecutive or equal to row numbers. |

The latitude, longitude, depth, and ID columns must be distinct. GraphSplit
preserves all other columns when writing catalogs.

## 6. Internal coordinates: `[coordinates]`

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

## 7. Travel-time source and lifecycle: `[travel_time]`

| Setting | Default | Meaning |
| --- | --- | --- |
| `travel_time.type` | `"lookup"` | `lookup` uses the native layered-model table. `constant` bypasses table building and is intended mainly for synthetic tests. |
| `travel_time.table_file` | `"lookuptable/graphsplit_tt.gstt"` | Native memory-mapped table path. The parent directory is created automatically. MATLAB `.mat` tables are not accepted. |
| `travel_time.velocity_model_file` | `"vm.txt"` | Layered P/S velocity model used to build or validate a lookup table. Its SHA-256 digest is stored in the table. |
| `travel_time.geometry` | `"auto"` | Geometry required when opening a table: `auto`, `cartesian`, or `radial`. `auto` reads an existing table's geometry. If a table must be built, `auto` delegates to `build_geometry`. An explicit geometry rejects/rebuilds a table of the other type. |
| `travel_time.build_geometry` | `"cartesian"` | Builder used when `geometry = "auto"` and no compatible table exists. `cartesian` is flat local range; `radial` uses spherical central angle. |
| `travel_time.auto_build` | `true` | Build the table if it is missing. When `false`, a missing table is an error. |
| `travel_time.auto_rebuild` | `true` | Rebuild an existing table when its model hash, geometry, Earth radius, spacing, or seed-catalog coverage is incompatible. When `false`, incompatibility is an error. |
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

## 8. Native lookup-table grid: `[lookup]`

Smaller spacing improves table resolution but increases build time and disk
size. Approximate stored size is
`2 phases × 4 bytes × Nrange × Ndepth × Nstation_depth`.

| Setting | Default | Meaning |
| --- | --- | --- |
| `lookup.horizontal_step_m` | `125.0` | Horizontal range spacing in metres. In radial mode this is a surface-arc spacing converted internally to angle. Changing it makes an existing table incompatible. |
| `lookup.depth_step_m` | `50.0` | Event-depth grid spacing in metres. It also sets the finest available station-depth spacing. |
| `lookup.station_depth_step_m` | `100.0` | Requested maximum spacing between station-depth slices. It must be at least `depth_step_m`; endpoint alignment may make actual spacing finer. More slices improve elevation interpolation but cost proportionally more time and disk. |
| `lookup.maximum_distance_km` | `0.0` | Table range. Zero derives it from the largest seed-event-to-station distance plus a margin. A positive explicit value must cover every seed event/station pair. This is unrelated to `graph.maximum_distance_km`. |
| `lookup.minimum_depth_km` | `0.0` | Shallow table limit. Zero automatically includes the shallowest model node, station, and seed event. A nonzero value is explicit; negative values are useful for stations above sea level. |
| `lookup.maximum_depth_km` | `0.0` | Deep table limit. Zero chooses an automatic limit; a positive value is explicit and must contain all seed events and stations. |
| `lookup.distance_margin_km` | `1.0` | Extra range used only when `maximum_distance_km = 0`. The builder uses at least two horizontal grid cells even if this margin is smaller. |
| `lookup.depth_margin_km` | `50.0` | Extra depth below the deepest seed event used by automatic `maximum_depth_km`. The automatic depth also considers range and velocity-model extent. Reduce cautiously if table size is excessive. |
| `lookup.maximum_sweeps` | `32` | Maximum fast-sweeping cycles for each phase and station-depth slice. Increase only if the builder reports non-converged slices. |
| `lookup.sweep_tolerance_s` | `1.0e-7` | Maximum travel-time change, in seconds, used to declare a slice converged. Smaller is stricter and may require more sweeps. |
| `lookup.threaded` | `false` | Build independent station-depth slices in parallel. Start Julia with multiple threads, for example `julia -t auto ...`; otherwise this switch has no speed effect. Parallel slices can increase peak memory and make progress messages less orderly. |
| `lookup.maximum_table_gib` | `8.0` | Safety limit on the estimated P+S cube size in GiB. |
| `lookup.allow_large_table` | `false` | Permit a table above `maximum_table_gib`. Set only after checking available disk, build time, and memory pressure. |

## 9. Stage 1 and Stage 2 solver settings

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

## 10. Observation selection: `[observations]`

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

## 11. Reference frame and pinned events: `[gauge]`

Relative times alone do not define every possible common shift of the model.
The gauge tells GraphSplit how to choose a stable reference frame. You do not
need linear algebra to use it: keep the default unless you see coherent drift
or have trusted event locations.

| Setting | Default | Meaning |
| --- | --- | --- |
| `gauge.mode` | `"zero_mean"` | `zero_mean` adds a weak/common-shift reference described below. `pin` holds selected event fields exactly fixed. |
| `gauge.zero_mean` | `"origin_time"` | Used in zero-mean mode. `origin_time` removes only the common relative-time shift. `xyzt0_scaled` also discourages a common spatial translation after converting metres to seconds with `reference_velocity_ms`; it is useful for a weakly anchored first pass. `xyzt0` is unscaled and should normally be avoided. |
| `gauge.constraint_weight` | 10.0 | Strength of the zero-mean constraint relative to observations. It is not the same as damping. Larger values enforce the selected mean update more strongly. The origin-time mean is also recentered exactly after each iteration. |
| `gauge.reference_velocity_ms` | 5000.0 | Metres-to-seconds scaling for `zero_mean = "xyzt0_scaled"`. It has no effect in other gauge choices. |
| `gauge.pin_event_ids` | `[]` | Persistent serial `EventID` values to hold fixed in pin mode, for example `[3, 50, 67]`. These are IDs from the configured catalog ID column—not row numbers. At least one is required when `mode = "pin"`. |
| `gauge.pin_fields` | `"xyz"` | Fields fixed for each listed event: any combination of `x`, `y`, `z`, and `t0`, such as `"xy"`, `"z"`, or `"xyzt0"`. `xyz` fixes the location but still solves relative origin time. |
| `gauge.pin_reference_catalog` | `""` | Optional catalog supplying trusted latitude/longitude/depth for pinned IDs. Empty means pin at the locations in `io.catalog_file`. The trusted file may contain only the pinned events and may use any row order, but it must use the same column layout configured under `[catalog]` and contain every listed serial ID. Trusted catalog date/time fields are not imported; a pinned `t0` is held at GraphSplit's initial relative value of zero. |

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

## 12. Choosing nearby event pairs: `[graph]`

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

## 13. Optional graph augmentation: `[graph.augmentation]`

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

## 14. Output switches: `[output]`

| Setting | Default | Meaning |
| --- | --- | --- |
| `output.write_filtered_catalogs` | `true` | Write `catalog_preloc_filt.txt` and `catalog_dd_filt.txt`, containing only events that participated in retained observations for the relevant stage. The full catalogs are always written. |
| `output.write_graph_metadata` | `true` | Write `catalog_dd_graphmeta.csv` with component, geometric degree, retained pair degree, observation counts, and activity flags. |
| `output.write_solver_history` | `true` | Write `solver_history.csv` with robust RMS, spatial/time step RMS, and PCG iteration count for every nonlinear iteration. |
| `output.write_run_summary` | `true` | Write `run_summary.toml` with input counts, travel-time type, graph size, observation counts, iterations, and final residuals. |

These switches affect files only, not the location solution.

## 15. Experimental theta-bias diagnostic: `[experimental.bias]`

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

## 16. Complete worked configurations

### Basic run with visible damping

Start from `config/graphsplit_template.toml`. Its two damping values are shown
explicitly so a user never has to discover them in source code.

### Iterated run

```toml
[io]
catalog_file = "../pass1/catalog_dd.txt"
output_dir = "pass2"

[run]
prelocation = false
```

You may also tighten graph radius/degree in the second pass. This is simply a
new run using the first output as its seed.

### Trusted partial pins

Use the three-line example in Section 11. The trusted file does not need to be
a duplicate of the complete working catalog.

### Reproducibility

For any published result, retain the actual TOML alongside the input catalog,
stations, velocity model, theta/thetaStd folders, `run_summary.toml`, graph
metadata, and solver history. The complete/default TOML documents software
defaults, but only the run-specific TOML records which values you changed.
