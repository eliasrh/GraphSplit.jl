# Configuration recipes

These are TOML fragments to merge into a copy of
`config/graphsplit_template.toml`. The runner and code are identical for every
case.

## Basic

Set the five input/output paths, choose `build_geometry`, and run. The minimal
template already selects Stage 1 followed by Stage 2.

## Common centroid when individual catalog locations are not trusted

```toml
[initialization]
mode = "common_centroid"

[run]
prelocation = true
```

All events begin at the mean input x, y, and depth. Stage 1 must separate them
before Stage 2 constructs the event graph.

## Manually specified common starting location

```toml
[initialization]
mode = "common_manual"
latitude = 64.0200
longitude = -21.2100
depth_km = 3.0

[run]
prelocation = true
```

The manual point is an initial guess, not a hard constraint. The input catalog
is still used for serial IDs and preserved columns, but its individual
hypocenters are overridden in memory.

## Second pass

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
```

This is the complete iterative mechanism: the previous DD output is the next
seed, and Stage 1 is disabled.

## Pin exact serial IDs at their seed locations

```toml
[gauge]
mode = "pin"
pin_event_ids = [101, 202, 303]
pin_fields = "xyz"
```

## Pin exact serial IDs to a trusted catalog

```toml
[gauge]
mode = "pin"
pin_event_ids = [101, 202, 303]
pin_fields = "xyz"
pin_reference_catalog = "trusted_catalog.txt"
```

The trusted catalog must contain those same IDs in its configured ID column,
but it does not need to contain any unpinned events. For example, it may be a
three-line standard catalog containing only IDs 101, 202, and 303. Matching is
by serial ID, not row order or row number; its location columns replace the
seed locations before those fields are fixed.

Do not move only one or a few events far from the swarm and expect these pins
to translate the complete catalog. Although the pins remain exact, Huber can
strongly downweight their inconsistent links and the Stage-2 graph can isolate
them. For a diagnostic test of that construction, use `pin_fields = "xyz"`,
set both Stage-1 and Stage-2 `huber_k` values very large, and check graph
components. Common initialization is the cleaner workflow when individual
seed locations are unavailable.

See the [complete configuration reference](../docs/CONFIGURATION_REFERENCE.md)
for a literal three-line trusted-catalog example and an explanation of every
setting.

## Radial table

```toml
[travel_time]
geometry = "auto"
build_geometry = "radial"
```

Once built, the table's radial geometry is detected from its header.

## Bootstrap location clouds

```toml
[uncertainty]
method = "bootstrap"

[uncertainty.bootstrap]
replicates = 100
resampling_unit = "station_phase"
summary_method = "percentile"
confidence_level = 0.95
write_samples = true
write_catalogs = false
```

The normal catalogs are unchanged. The four `bootstrap_samples_*.txt` tables
contain the all-data location plus every resampled realization for each
filtered serial ID. Set `write_catalogs = true` only if complete catalog-format
files are also needed for every replicate.

## Fast formal uncertainty for a large catalog

```toml
[uncertainty]
method = "linearized"

[uncertainty.linearized]
probes = 12
```

This writes `linerrxyz.txt`. Use `method = "both"` to produce it together with
the bootstrap products.
