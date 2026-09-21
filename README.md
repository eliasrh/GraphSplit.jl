# GraphSplit.jl

GraphSplit relocates earthquake catalogs using the synchronized `theta` and
`thetaStd` files produced by [DDSync](https://github.com/eliasrh/DDSync). It is
designed for large catalogs: the event graph is sparse, the nonlinear solve is
matrix-free, and native travel-time tables are memory-mapped.

GraphSplit uses only Julia standard libraries, supports Cartesian and radial
layered-Earth travel times, and is controlled by one TOML file.

## Manuscript

**GraphSplit: Sparse Double-Difference Earthquake Relocation From Synchronized
Differential-Time Graphs** (tentative title). Elías Rafn Heimisson.
Submission year: **2026**.

Software citation metadata are in [CITATION.cff](CITATION.cff). The manuscript
citation will be updated when publication details are available.

## Start here

The repository documentation is organized by task:

1. [Manual](docs/MANUAL.md) — the complete workflow, beginning with DDSync and
   ending with interpretation and reproducibility.
2. [Configuration reference](docs/CONFIGURATION_REFERENCE.md) — every TOML
   setting, its units, interactions, and when to change it.
3. [Input and output file formats](docs/FILE_FORMATS.md) — literal catalog,
   station, theta, velocity-model, constraint, and uncertainty formats.
4. [Configuration recipes](examples/README.md) — copy-ready TOML fragments for
   common seeds, second passes, pins, fixed depths, depth bounds, radial tables,
   and uncertainty estimates.
5. [Validation guide](docs/VALIDATION.md) — tests and checks to perform before
   interpreting a field catalog.

New users should read the first three in that order. The minimal
[`graphsplit_template.toml`](config/graphsplit_template.toml) contains the
settings most likely to be changed. The commented
[`graphsplit_complete.toml`](config/graphsplit_complete.toml) exposes every
available option.

## Before GraphSplit: run DDSync

GraphSplit does not read `dt.cc` directly. First run DDSync and retain:

```text
theta/
thetastd/
```

The catalog supplied to GraphSplit must preserve DDSync's persistent serial
`EventID` values. Matching is by ID, never by row number or row order.

## Five-minute setup

Requirements: Julia 1.10 or later and a completed DDSync run.

1. Copy the minimal configuration beside the dataset:

   ```bash
   cp /path/to/GraphSplit.jl/config/graphsplit_template.toml graphsplit.toml
   ```

2. Arrange or point the TOML to these inputs:

   ```text
   my_run/
     graphsplit.toml
     catalog.txt
     stations.txt
     vm.txt
     theta/
     thetastd/
   ```

3. In `graphsplit.toml`, check at least:

   | Question | Setting |
   | --- | --- |
   | Are the input and output paths correct? | `[io]` |
   | Are individual catalog locations usable? | `initialization.mode` |
   | Is the network local or geographically broad? | `travel_time.build_geometry` |
   | Is the lookup spacing adequate? | `[lookup]` |
   | How local should Stage-2 pairs be? | `[graph]` |
   | How much damping is appropriate? | `prelocation.damping_lambda`, `relocation.damping_lambda` |
   | Is a physical depth limit or known depth required? | `[constraints.depth_bound]`, `[constraints.fixed_depth]` |
   | Are uncertainty samples needed? | `uncertainty.method` |

4. Run GraphSplit:

   ```bash
   julia --project=/path/to/GraphSplit.jl \
     /path/to/GraphSplit.jl/run_graphsplit.jl \
     /path/to/my_run/graphsplit.toml
   ```

5. Inspect `solver_history.csv`, `catalog_dd_graphmeta.csv`, and
   `run_summary.toml` before interpreting `catalog_dd_filt.txt`.

No package installation is needed. `--project` selects this repository's Julia
environment. Paths inside the TOML are resolved relative to the TOML file, not
the terminal's current directory.

## Required inputs

- `catalog.txt`: numeric catalog. Defaults use latitude, longitude, and depth
  in columns 7, 8, and 9 and the positive serial `EventID` in the final column.
- `stations.txt`: `STA latitude longitude elevation_m`, or
  `NET STA latitude longitude elevation_m`.
- `theta/theta_<STA>_<P|S>.txt`: synchronized DDSync theta values.
- `thetastd/std_theta_<STA>_<P|S>.txt`: theta uncertainty and normally DDSync
  graph degree.
- `vm.txt`: `depth_km Vp_km_s Vs_km_s`, or the radial
  `depth_km radius_km Vp_km_s Vs_km_s` form.

See [FILE_FORMATS.md](docs/FILE_FORMATS.md) before adapting a catalog or
velocity model.

## What the two stages do

- Stage 1 uses theta-reference star observations to improve or construct the
  starting catalog.
- Stage 2 builds a sparse local event graph from the Stage-1 solution and
  performs double-difference relocation on its retained event pairs.

The default `initialization.mode = "catalog"` uses individual input
hypocenters. `common_centroid` and `common_manual` start every event together;
they require Stage 1 so the event graph is never constructed from coincident
points.

On the first lookup-model run, GraphSplit builds a native `.gstt` travel-time
table. Later runs validate its velocity-model hash, geometry, spacing, and
coverage before reuse. To build or validate only the table:

```bash
julia --project=. build_travel_times.jl graphsplit.toml
julia --project=. build_travel_times.jl graphsplit.toml --force
```

## Depth constraints

Depth constraints are independent of the gauge and apply in both relocation
stages and bootstrap replicates.

To prevent events from becoming shallower than a physical surface:

```toml
[constraints.depth_bound]
enabled = true
minimum_depth_km = -0.8
scope = "all"
reflected_prelocation_restart = true
pilot_shallow_margin_km = 10.0
mirror_plane = "stations_median"
```

With the normal positive-down convention, `-0.8` km represents 800 m elevation.
The solver uses an active bound rather than travel-time-table clamping. The
optional reflected restart first searches for the sparse-network mirror branch,
reflects forbidden Stage-1 solutions about the selected plane, reruns bounded
Stage 1, and only then constructs the event graph. If the pilot exceeds the
reserved lookup range, increase `pilot_shallow_margin_km` and rebuild the table.

To locate explosions horizontally and in origin time while keeping their depth
known exactly:

```toml
[constraints.fixed_depth]
enabled = true
scope = "all"                # or "event_ids"
event_ids = []                # required only for event_ids scope
depth_km = -0.8
```

Constraint status is written to `depth_constraint_status.csv`; the scientific
catalog columns remain unchanged. A boundary-active event has unresolved
one-sided depth, not a measured location exactly on the boundary.

## Principal outputs

The output directory contains:

- `catalog_preloc.txt` and `catalog_preloc_filt.txt`;
- `catalog_dd.txt` and `catalog_dd_filt.txt`;
- `catalog_preloc_dxdydzt0.txt` and `catalog_dd_dxdydzt0.txt`, with cumulative
  east/north/vertical/time shifts for every serial ID;
- `catalog_dd_graphmeta.csv`, with graph and observation support per event;
- `solver_history.csv`, with nonlinear and inner-solver convergence;
- `run_summary.toml`, with run settings and core counts;
- `depth_constraint_status.csv` when a physical or fixed depth is enabled.

For the default HypoDD-style catalog, the relocated files apply the solved
origin-time correction to year/month/day/hour/minute/second, including calendar
rollover. Positive `dt0_s` moves the origin time later. Auxiliary columns are
preserved. If a catalog has no usable calendar time, set
`catalog.origin_time_columns = []`; its time-like columns remain unchanged while
the shift sidecars still report `dt0_s`.

The shift files contain `dx_m dy_m dz_m dt0_s EventID`. Events that have never
participated in a retained relocation have four zeros. The values are cumulative
from the first-pass input, which makes them useful for displacement analysis and
also preserves the relative-time state required by an optional DD-only second
pass. See [file formats](docs/FILE_FORMATS.md#relocation-shift-sidecars).

With `uncertainty.method = "linearized"`, `"bootstrap"`, or `"both"`, GraphSplit
writes separate sidecars without changing any catalog format. Bootstrap output
includes every longitude, latitude, depth, and relative-origin-time sample in
wide serial-ID-keyed tables. See the [manual](docs/MANUAL.md#10-uncertainty-estimates)
for interpretation.

## Important workflow choices

- Do not move one or a few events far from the swarm and expect hard pins to
  translate the catalog. Large residuals can be Huber-downweighted, and the
  displaced pins can become graph-isolated.
- Use `common_centroid` or `common_manual` when individual seed locations are
  unavailable or intentionally ignored.
- A DD-only second pass is an advanced option for a poor, sparsely connected
  seed: pair a previous `catalog_dd.txt` with its
  `catalog_dd_dxdydzt0.txt`, then rebuild a tighter graph after the broad first
  graph has improved event proximity.
- Use fixed depth when depth is known externally; do not imitate it with very
  large depth damping.
- A cluster at the physical depth bound indicates unresolved depth. Inspect
  `depth_constraint_status.csv` and bootstrap bound-active fractions.
- `lookup.minimum_depth_km` controls table coverage only. It is not a physical
  event bound.

## Testing

Run:

```bash
julia --project=. -e 'using Pkg; Pkg.test()'
```

GitHub Actions tests Julia 1.10 and the current stable Julia release. Scientific
validation remains dataset-specific; follow [VALIDATION.md](docs/VALIDATION.md).

## Additional tools

The [benchmark comparator](benchmark/yifan2025/README.md) reports
horizontal/depth accuracy, local neighbor-pair precision, and point-cloud
Chamfer distance. Its README defines the metrics, including the Chamfer
normalization, and explains how to obtain the external true locations.
The coherent theta-bias scan remains experimental and disabled by default.

## Benchmark data attribution

The example catalog, station coordinates and velocity model in
`benchmark/yifan2025/input/` come from the synthetic experiment of:

Yu, Y., Ellsworth, W. L., and Beroza, G. C. (2025). *Accuracy and Precision of
Earthquake Location Programs: Insights from a Synthetic Controlled Experiment*.
**Seismological Research Letters, 96**(3), 1860–1874.
[https://doi.org/10.1785/0220240354](https://doi.org/10.1785/0220240354).

Please cite that paper when using these example data. The authors' code and
source data are available from
[Yu's benchmark repository](https://github.com/YuYifan2000/comparison_hypoDD_GrowClust).
True locations are not distributed here; obtain the matching realization
through that repository's instructions. The example inputs do not include the
DDSync potentials needed to run a relocation.

## Citation and license

Cite the GraphSplit manuscript when publication details are available, and
identify the software version or commit used. Software metadata are in
[CITATION.cff](CITATION.cff). GraphSplit is distributed under the
[GraphSplit Non-Commercial License v1.0](LICENSE). Benchmark source-data
attribution is separate from the software license; retain the Yu et al.
citation when reusing those inputs.
