# Experimental relocation with a fixed 3D model

This branch adds relocation in a prescribed three-dimensional velocity model.
It does not estimate velocity. The submitted GraphSplit study uses the existing
1D travel-time calculation; its configurations and results are unchanged.
The 3D implementation has small synthetic and numerical tests, but has not been
validated on large field catalogs. Start with a small subset and check grid
refinement before interpreting a fault geometry.

## What changes in a relocation?

A 1D model has the same velocity at a given depth everywhere. Its travel times
can be stored in a compact table indexed by event–station distance and depth.
A 3D model allows velocity to vary horizontally as well. GraphSplit instead
needs a volume of travel times from each station, separately for P and S.

The embedded fast-marching method solves the isotropic eikonal equation on a
regular Cartesian grid. It propagates the earliest arrival from a station
through the supplied velocity model. Reciprocity gives the event-to-station
travel time. The default uses second-order upwind differences where accepted
neighbors permit them and first-order differences elsewhere. It is a different
forward calculation from the 1D fast-sweeping lookup builder. No NonLinLoc,
GrowClust3D or TomoSplit installation is required.

At each location step, GraphSplit interpolates the time and its three spatial
derivatives at each event. These derivatives describe how the predicted time
changes when the earthquake moves east, north or down. They are still needed
for the linearized relocation and are stored in the observation arrays as in
1D. The code does **not** retain derivatives with respect to velocity at every
model node; those would be needed for tomography. The event graph, robust
weights, reference conditions and location solver otherwise work as before.

## Two ways to supply travel times

1. **Velocity grids:** give P and S NonLinLoc `.hdr`/`.buf` model pairs. GraphSplit
   samples them onto its forward grid and computes the travel-time volumes.
2. **Precomputed NonLinLoc `TIME` grids:** supply one `.hdr`/`.buf` pair per
   station and phase. This is the travel-time-grid route used by
   [GrowClust3D](https://github.com/dttrugman/GrowClust3D.jl/wiki/Documentation).
   GraphSplit reads these times directly, retaining their native spacing and
   extent. It does not run FMM again or reproduce NonLinLoc's model builder.

The input formats follow the [NonLinLoc implementation](https://github.com/ut-beg-texnet/NonLinLoc).
This first version supports its **SIMPLE** geographic projection, including
rotation, or explicitly declared local Cartesian coordinates. Other projections
are rejected, rather than treated as equivalent. `TIME2D`, anisotropy, converted
phases, reflected arrivals and spherical propagation are not supported.

## Start with a velocity model

Copy [the 3D template](../config/graphsplit_3d.toml) beside a dataset that already
has a catalog, stations and the matching DDSync theta directories. Set:

```toml
[travel_time]
type = "3d"

[grid3d]
model_format = "nll_velocity"
vp_file = "model.P.mod.hdr"
vs_file = "model.S.mod.hdr"
spacing_m = [500.0, 500.0, 500.0]
```

Then use the usual commands:

```bash
julia --project=/path/to/GraphSplit.jl /path/to/GraphSplit.jl/build_travel_times.jl graphsplit_3d.toml
julia --project=/path/to/GraphSplit.jl /path/to/GraphSplit.jl/run_graphsplit.jl graphsplit_3d.toml
```

Paths are relative to the TOML file. The optional builder lets you check model
coverage and memory before relocation. `--force` rebuilds the times. A completed
cache records the model, surface, station coordinates, settings and implementation
hashes; input changes create a different cache. Cached fields are checked for
size and checksum before reuse. Incomplete fields are not used.

The forward grid covers the input-model box unless `bounds_km` selects a smaller
box, in **model coordinates**. Each axis is divided into equal intervals no
larger than its requested spacing. The exact spacing and dimensions are reported
in `run_summary.toml`. Coarsening is explicit; it never extends the model.

## Coordinates, elevation and topography

Catalog depths and the grid's z coordinate are positive downward in kilometres.
Station elevations are positive upward in metres, relative to the **same
vertical datum**. GraphSplit requires `event_vertical = "positive_depth"` and
`station_vertical = "depth_from_elevation"` in this mode. A station at 800 m
elevation has z = −0.8 km, so the model must extend above zero depth.

For a header with `TRANSFORM SIMPLE`, the reader uses the NonLinLoc geographic
conversion, including its latitude-dependent longitude scale and rotation.
It also transforms the derivatives back into GraphSplit's local east/north
coordinates. Using the same origin alone is not sufficient to equate the two
coordinate formulas.

A header with no transform or `TRANSFORM NONE` is accepted only with
`grid3d.coordinate_system = "local"` and `coordinates.reference = "manual"`.
Here the model axes must be east, north and down in GraphSplit's local frame:
`x = R cos(latitude_origin) Δlongitude`, `y = R Δlatitude`, with angles in
radians. State the reference latitude and longitude explicitly. Do not ignore
a declared projection to make an incompatible model fit.

Station elevation alone does **not** stop a ray from passing through air. To
exclude air, set `surface_file` to a text file with one row per surface node:

```text
# x_km y_km elevation_m
-5.0 -5.0 200.0
-5.0 -4.5 220.0
```

The full file must contain a complete regular x/y grid covering the model box,
in the same horizontal coordinates and vertical datum. Elevation is bilinearly
interpolated onto the forward grid. Nodes above this surface are excluded from
propagation. Short station-to-grid segments are checked against that surface.
A station may lie on or below it, including a buried receiver. A station above
it is rejected: check the DEM, elevation precision and vertical datum. The code
does not move stations to make them fit. For imported TIME grids, a surface can
mask interpolation cells, but cannot repair rays previously calculated through
air; the original travel-time calculation must already use a suitable model.

This is a regular-grid approximation to terrain, not a mesh that follows the
surface exactly. A usable event interpolation cell needs eight finite solid
corners. Very shallow events and steep slopes therefore require finer spacing.
Positive velocity must be supplied throughout the input box; values in air are
not propagated when a surface is supplied. Without a surface file, the whole
box is treated as material, including its negative-depth part.

All stations and starting events must lie within the model. The inversion never
clamps an event to a grid edge or extrapolates a travel time. If an update would
leave the usable domain, it is shortened along the same direction and a warning
is printed. A severely restricted step stops with an error. Boundary warnings
are a reason to enlarge/refine the grid and rerun, not evidence of convergence.
The reflected-depth prelocation restart is disabled in 3D because its pilot can
leave the fixed model domain.

## Layers and grid resolution

NonLinLoc model files contain values at grid nodes; a sharp layer is represented
at the resolution of those nodes. `model_interpolation = "nearest"` is the
default when transferring a model to a different grid. It avoids smoothing a
velocity jump during this transfer. `"slowness_linear"` interpolates reciprocal
velocity and is suitable for a smooth model. Neither option preserves a thin
layer that the forward grid fails to sample. A discontinuity may move by part
of a cell when grids differ, and the finite-difference solution introduces its
own discretization error. This is not an interface-fitted ray tracer.

Use a grid that samples important layers and terrain, then halve the spacing
for a smaller test. Compare travel times **and relocated locations**. A coarse
3D grid can add errors comparable to the differential-time uncertainty, even
when the input velocity model is accurate. The default 500 m spacing is a
memory-conscious starting point, not a claim that 500 m is adequate for a
particular sequence.

## Precomputed times

```toml
[travel_time]
type = "3d"

[grid3d]
model_format = "nll_time"
time_root = "time/model"
byte_order = "little"
phases = ["P", "S"]
```

For station `ABC`, the expected pairs are `time/model.P.ABC.time.hdr/.buf`
and `time/model.S.ABC.time.hdr/.buf`. The header station name and source
coordinates must agree with `stations.txt` (coordinate tolerance 1 cm).
All grids must share their dimensions, projection, origin and spacing.
`spacing_m`, `model_interpolation` and `accuracy_order` do not change imported
times; `bounds_km` must be empty. Use `phases = ["P"]` for P-only theta inputs.
A missing station or phase field is an error.

## Memory and disk space

For N grid nodes and G station–phase fields, retained Float32 times occupy
approximately **4 N G bytes**, both on disk and as addressable mapped data.
Memory mapping does not make these volumes free of RAM cost. FMM builds one
field at a time, using additional arrays proportional to N. Its trial-node
queue has at most N entries; it does not grow with duplicate queue records.

For a 100 × 100 × 30 km box and 50 stations with both phases:

| Spacing | Nodes | Retained times alone |
| --- | ---: | ---: |
| 500 m | 2,464,461 | 0.92 GiB |
| 250 m | 19,456,921 | 7.25 GiB |

Halving all three spacings increases volume storage by roughly eight times.
More stations or phases increase it linearly. This differs from the compact 1D
table, which can be shared by stations with the same elevation sampling.

Before allocating the forward grid, GraphSplit reports retained-field storage
and a conservative build estimate: `4 N G + 80 N + largest input buffer +
16 × surface-file bytes`. The extra terms allow for serial FMM work, model
sampling, masks, conversion and surface parsing. The estimate is deliberately
conservative for precomputed TIME imports. It excludes Julia's runtime,
catalogs, observations, location-solver arrays and other applications; it is
not a process-wide RAM guarantee.

`maximum_memory_gib = 4.0` refuses a larger estimate before loading arrays;
`warning_memory_gib = 1.0` prints an earlier warning. Increase limits only after
checking available RAM. Also allow disk space for the full retained fields and
one temporary field. The cache directory is safe to remove when no run is
using it, since it can be rebuilt from the input model. Old cache versions are
kept when inputs change and can be deleted selectively.

## Reproduce the small test

```bash
julia --project=. examples/three_dimensional/run_synthetic.jl
```

This writes inputs, two run TOMLs, relocated catalogs and error metrics under
`examples/three_dimensional/generated/`. See the [example guide](../examples/three_dimensional/README.md)
for the model, noise and test results. This is a known-model recovery test,
not a comparison of field performance with GrowClust3D or NonLinLoc.
