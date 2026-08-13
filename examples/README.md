# Configuration recipes

These are TOML fragments to merge into a copy of
`config/graphsplit_template.toml`. The runner and code are identical for every
case.

## Basic

Set the five input/output paths, choose `build_geometry`, and run. The minimal
template already selects Stage 1 followed by Stage 2.

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
