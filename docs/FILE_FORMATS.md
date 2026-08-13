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

For example:

```text
2025 01 03 12 14 05.230 35.410000 -117.910000 4.2000 1.2 3
```

With the defaults, this event has latitude 35.410000°, longitude −117.910000°,
depth 4.2 km, and serial ID 3. Columns 1–6, magnitude, and any other auxiliary
columns are preserved, but GraphSplit does not use catalog origin times as
absolute constraints. Its relative origin-time corrections are internal
nuisance parameters and are not written back to columns 1–6.

Output preserves the original column count and order, replacing only the three
location columns. Filtered output preserves the original serial IDs.

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

## Benchmark truth

`benchmark/yifan2025/compare_catalogs.jl` expects:

```text
latitude longitude depth_km [optional columns...]
```

Catalog serial ID `k` selects truth row `k`. The benchmark utility checks that
every requested row exists.
