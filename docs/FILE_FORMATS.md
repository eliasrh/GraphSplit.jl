# Input file formats

A normal run needs five kinds of input: a catalog, a station file, a P/S
velocity model, a directory of DDSync theta files, and normally a directory of
DDSync thetaStd files. The TOML tells GraphSplit where each one is located.

```text
my_run/
  graphsplit.toml
  catalog.txt
  stations.txt
  vm.txt
  theta/
    theta_STA1_P.txt
    theta_STA1_S.txt
  thetastd/
    std_theta_STA1_P.txt
    std_theta_STA1_S.txt
```

All text inputs are whitespace-separated. Catalog and numeric readers also
accept commas. Blank lines and trailing `#` or `%` comments are ignored. A
single textual header before the numeric data is accepted in numeric files.

The persistent serial `EventID` is the join key. GraphSplit does not assume
that row 50 is event 50: it explicitly looks up the ID written in each file.
This allows catalogs and trusted pin files to be reordered or filtered.

## Catalog

The default DDSync/HypoDD-style layout is:

```text
YYYY MM DD HH MM SS.sss latitude longitude depth_km magnitude ... EventID
```

Defaults select columns 7/8/9 for latitude/longitude/depth and the final column
for `EventID`. Override these one-based indices under `[catalog]`. IDs must be
unique positive integers, but they need not be contiguous or match row order.
They must not exceed 9,007,199,254,740,991 (`2^53 - 1`), the supported exact
integer range of the numeric file reader. GraphSplit rejects larger IDs
instead of rounding them. This applies to catalogs, theta event and reference
IDs, theta uncertainty event IDs, and restart-shift event IDs.

If DDSyncJulia reported automatic ID mapping, use its `catalog_seq.txt` and
the matching theta directories. The external IDs in the original catalog
and `dt_sync.cc` then belong to a different numbering system. Keep
`event_id_map.csv` for translating GraphSplit results back to those external
IDs. Reordering rows with their IDs attached is allowed; independently
renumbering the catalog after synchronization is not. When preparing theta
files outside DDSync, apply any ID mapping consistently to every input,
including reference IDs and event-specific constraint settings.

For example:

```text
2025 01 03 12 14 05.230 35.410000 -117.910000 4.2000 1.2 3
```

With the defaults, this event has latitude 35.410000°, longitude −117.910000°,
depth 4.2 km, and serial ID 3. GraphSplit does not treat the input origin time as
an absolute constraint, but it adds the solved relative `t0` correction to
columns 1–6 when writing the relocated catalog. Positive `dt0_s` moves the
origin later; year, month, and day rollover are handled as Gregorian calendar
time. Magnitude and all other auxiliary columns are preserved.

`catalog.origin_time_columns` identifies year/month/day/hour/minute/second and
defaults to `[1, 2, 3, 4, 5, 6]`. Set it to `[]` for a point catalog or any input
without usable calendar timing. Those fields are then copied unchanged, while
the shift sidecar remains the authoritative `dt0_s` output. An active row with
invalid calendar fields is likewise preserved with a warning. Output seconds
use `catalog.origin_time_decimals` digits after the decimal (six by default).

Output preserves the original column count and order, replacing the three
location columns and, when configured, the six origin-time columns. Filtered
output preserves the original serial IDs.

With `initialization.mode = "common_centroid"` or `"common_manual"`, the
latitude, longitude, and depth columns remain syntactically required but are
overridden as individual in-memory seeds. The serial ID and all preserved
auxiliary columns are still taken from this file. Output location columns
contain the relocated results, not the common starting point.

Uncertainty estimates never add columns to a catalog. They are written as
separate text sidecars keyed by the same persistent `EventID`; their formats
are described below.

## Stations

Either form is accepted:

```text
STA latitude longitude elevation_m
NET STA latitude longitude elevation_m
```

Station codes must match the `<STA>` part of theta filenames exactly.

Examples:

```text
STA1 35.5000 -117.8000 820.0
CI STA2 35.5500 -117.7500 1010.0
```

Elevation is in metres, positive upward in the normal convention. With
`coordinates.station_vertical = "depth_from_elevation"`, 820 m elevation
becomes an internal station depth of −820 m. Network code is accepted for
compatibility but station matching uses `STA1` or `STA2`, not `CI.STA2`.

## DDSync theta

Filename:

```text
theta_<STA>_<P|S>.txt
```

Rows:

```text
EventID theta_s refEventID
```

Non-finite rows are skipped. Values are stored sparsely by serial ID.

For example:

```text
3  0.18423  1
50 0.39210  1
67 0.00000 67
```

The first row says that event ID 3 has synchronized potential 0.18423 s in
this station–phase group and uses event ID 1 as its DDSync reference. The
reference may differ between connected DDSync groups. In Stage 2, GraphSplit
forms a pair only when both events have theta values tied to the same reference.

The station and phase come from the filename. A station name may contain
underscores; the final underscore-separated token must be `P` or `S`.

## DDSync thetaStd

Filename:

```text
std_theta_<STA>_<P|S>.txt
```

Base rows:

```text
EventID std_theta_s refEventID degree
```

DDSync's optional fifth and sixth weight columns are accepted and ignored by
the core relocation. Matching is by `EventID`, not theta row.

For example:

```text
3 0.0031 1 12
```

means event 3 has theta standard deviation 0.0031 s, reference ID 1, and
DDSync degree 12. GraphSplit combines the two event uncertainties with
`hypot(sigma_i, sigma_j)` before weighting a difference.

A missing thetaStd file does not stop a run. That entire station–phase group
then uses the configured sigma floor, and `minimum_theta_degree` cannot be
applied because no degree values exist. Supplying thetaStd is therefore
strongly recommended for scientific runs.

## Velocity model

Three-column form:

```text
depth_km Vp_km_s Vs_km_s
```

Four-column radial-model form:

```text
depth_km radius_km Vp_km_s Vs_km_s
```

For the four-column form, `depth + radius` estimates Earth radius. Repeated
depths represent a discontinuity: the first velocity is the value approaching
from above and the final repeated value is used at and below the interface.

Example three-column model:

```text
0.0  5.50 3.20
2.0  5.80 3.35
8.0  6.20 3.55
20.0 6.60 3.80
```

All depths are kilometres positive downward; velocities are km/s. The model
must cover the table depth range sensibly. GraphSplit extends the sampled end
velocities outside the listed nodes, so an accidentally shallow model may run
without representing the intended deeper structure.

## Trusted pin catalog

`gauge.pin_reference_catalog` uses the same format and configured column
numbers as the working catalog. It may be a complete catalog, but it normally
needs only one row for each ID in `gauge.pin_event_ids`.

For example, with:

```toml
[gauge]
mode = "pin"
pin_event_ids = [3, 50, 67]
pin_fields = "xyz"
pin_reference_catalog = "catalog_trusted.txt"
```

the trusted file may contain only three rows, in any order. Each final serial
ID must be exactly 3, 50, or 67. GraphSplit copies the requested latitude,
longitude, and/or depth into the starting state and then holds those fields
fixed. It does not import trusted catalog origin time into the internal `t0`.

Without `pin_reference_catalog`, pins remain at the initialized in-memory seed.
That is normally their `io.catalog_file` location, but under a common
initialization mode it is the shared centroid or manual point.

## Fixed-depth reference catalog

`constraints.fixed_depth.reference_catalog` uses the same numeric catalog
format and configured depth/ID columns. It may contain only the selected IDs
and may use any row order. GraphSplit reads only the serial ID and depth; its
latitude, longitude, date/time, and auxiliary values are ignored.

For example, with `scope = "event_ids"` and `event_ids = [101, 202]`, a
two-row file is sufficient. Every selected ID must appear exactly once in the
reference catalog. When this path is empty, `constraints.fixed_depth.depth_km`
supplies one common depth.

## Native `.gstt` table

The binary format is internal but stable within format version 1. It contains:

- 16-byte magic and version;
- geometry plus a completion/scalar marker;
- travel-time and coordinate Earth radii;
- axis lengths and the SHA-256 velocity-model digest;
- Float64 range/angle, event-depth, and station-depth axes;
- contiguous Float32 P and S cubes in `(q, z, zs)` order.

Do not edit the file. GraphSplit rejects unsupported or truncated headers and
rebuilds incompatible tables when configured to do so. The header remains
marked incomplete until both P and S cubes have been flushed, so an interrupted
build cannot be mistaken for a valid table on the next run.

## Relocation-shift sidecars

Every run writes:

```text
catalog_preloc_dxdydzt0.txt
catalog_dd_dxdydzt0.txt
```

Each contains every input event in input-catalog order:

```text
# dx_m dy_m dz_m dt0_s EventID
12.450000 -3.200000 8.750000 0.014230000000 3
0.000000 0.000000 0.000000 0.000000000000 50
```

`dx_m` is local east and `dy_m` is local north. `dz_m` follows the configured
internal `coordinates.event_vertical` convention; under the default
`positive_depth`, positive `dz_m` means deeper. `dt0_s` is the origin-time
correction, with positive values meaning later. The coordinate reference and
vertical convention are recorded in `run_summary.toml`.

The prelocation file describes the Stage-1 state and the DD file describes the
final state. Shifts are cumulative relative to the catalog supplied to the
first pass. An event that has never participated in retained observations is
written as `0 0 0 0 EventID`. On a later DD-only pass, an event without new
support retains its earlier cumulative shift rather than being mislabeled as
never relocated. As with the relocation itself, a common time shift is fixed by
the selected gauge, so `dt0_s` is meaningful under that recorded gauge.

For an iterated DD graph, set `io.catalog_file` to a prior `catalog_dd.txt` and
`io.restart_shift_file` to the `catalog_dd_dxdydzt0.txt` from the same run. The
restart reader joins by `EventID`, so row order may differ and the shift file may
contain additional IDs, but every event in the working catalog must be present
exactly once. GraphSplit restores cumulative `t0` for its theta prediction and
adds only the newly solved time increment to the already-corrected input
calendar. Use the same stations and coordinate convention if cumulative spatial
shifts are to remain directly comparable across passes.

## Depth-constraint sidecars

When either depth constraint is enabled, `depth_constraint_status.csv` contains:

```text
EventID,bound_applies,fixed_depth,bound_depth_km,fixed_depth_km,reflected_prelocation,bound_hits_prelocation,bound_active_prelocation,bound_hits_relocation,bound_active_relocation,final_depth_km
```

It contains the union of bound-selected and fixed-depth events. Boolean fields
are written as 0/1. `bound_hits_*` counts nonlinear iterations in which the
active-set solver had to prevent a crossing. `bound_active_* = 1` means the
final step ended exactly on the physical boundary. `NaN` in an inapplicable
configured-depth field means that constraint does not apply to the event.

This file is deliberately separate from all four scientific catalogs.

## Uncertainty sidecars

`linerrxyz.txt` is whitespace-separated and contains:

```text
EventID std_x_m std_y_m std_z_m cov_xx_m2 cov_xy_m2 cov_xz_m2 cov_yy_m2 cov_yz_m2 cov_zz_m2
```

The x/y axes are local east/north and z follows the configured internal event
vertical convention. Only events in `catalog_dd_filt.txt` are written.
An exactly fixed z has zero conditional variance. For an event whose final
solution is active at the one-sided physical bound, z standard deviation and
z covariance terms are `NaN`; the x/y block remains the covariance conditional
on the active bound.

`booterrxyz.txt` contains one row per filtered ID:

```text
EventID n_valid valid_fraction center_dx_m lower_dx_m upper_dx_m ... cov_zz_m2
```

The d-values are local offsets from the all-data `catalog_dd` solution. The
meaning of lower/upper is recorded in `bootstrap_metadata.toml` and is either a
central percentile interval or mean plus/minus a selected standard-deviation
multiple. The six covariance fields are empirical sample covariances in m².

Each wide sample table has this form:

```text
EventID all_data sample_0001 sample_0002 ... sample_NNNN
```

The four files hold longitude (degrees), latitude (degrees), depth (km), and
internal relative origin-time adjustment (s). `NaN` records an inactive event
or a nonconverged replicate. The first `all_data` value is the ordinary solution
using every retained station-phase group, not a bootstrap draw.

`bootstrap_block_counts.txt` has one row per resampling block and one integer
multiplicity per replicate. With station-phase resampling, its label is the
theta filename; with station resampling, it is the station code. This small file
records exactly which information was omitted, retained once, or repeated in
every bootstrap realization.

With a physical depth bound, `bootstrap_depth_bound_status.txt` contains:

```text
EventID n_valid n_bound_active bound_active_fraction
```

The denominator includes only converged realizations in which that event has
resampled support. A high active fraction means the saved depth cloud is
strongly truncated by the physical prior.

## Benchmark truth

`benchmark/yifan2025/compare_catalogs.jl` expects:

```text
latitude longitude depth_km [optional columns...]
```

Catalog serial ID `k` selects truth row `k`. The benchmark utility checks that
every requested row exists.

The true-location file is not distributed with GraphSplit. See the
[benchmark data source and citation](../benchmark/yifan2025/README.md#data-source-and-citation)
for retrieval instructions.

## Experimental 3D velocity and travel-time grids

The [3D guide](THREE_DIMENSIONAL_MODELS.md) defines the NonLinLoc model and TIME
inputs, the four-column ASCII velocity-grid option, coordinate conventions,
surface-file format and cache lifecycle. A
3D run replaces `vm.txt` with the files selected under `[grid3d]`. Catalogs,
stations and theta inputs keep their existing formats.
