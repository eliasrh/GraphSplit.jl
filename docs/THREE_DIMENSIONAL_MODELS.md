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
forward calculation from the 1D fast-sweeping lookup builder. No NonLinLoc or
GrowClust3D installation is required.

At each location step, GraphSplit interpolates the time and its three spatial
derivatives at each event. These derivatives describe how the predicted time
changes when the earthquake moves east, north or down. They are still needed
for the linearized relocation and are stored in the observation arrays as in
1D. The code does **not** retain derivatives with respect to velocity at every
model node; those would be needed for tomography. The event graph, robust
weights, reference conditions and location solver otherwise work as before.

## Input choices

1. **Velocity grids:** give P and S NonLinLoc `.hdr`/`.buf` model pairs. GraphSplit
   samples them onto its forward grid and computes the travel-time volumes.
2. **Precomputed NonLinLoc `TIME` grids:** supply one `.hdr`/`.buf` pair per
   station and phase. This is the travel-time-grid route used by
   [GrowClust3D](https://github.com/dttrugman/GrowClust3D.jl/wiki/Documentation).
   GraphSplit reads these times directly, retaining their native spacing and
   extent. It does not run FMM again or reproduce NonLinLoc's model builder.
3. **Plain-text velocity grids:** supply four-column ASCII files with
   latitude/longitude/depth/velocity, or local x/y/depth/velocity. Rows may be
   unordered, but must form a complete regularly spaced grid. The next section
   gives an example.

The binary input formats follow the [NonLinLoc implementation](https://github.com/ut-beg-texnet/NonLinLoc).
This first version supports its **SIMPLE** geographic projection, including
rotation, or explicitly declared local Cartesian coordinates. Other projections
are rejected, rather than treated as equivalent. `TIME2D`, anisotropy, converted
phases, reflected arrivals and spherical propagation are not supported.

## What is in a NonLinLoc grid file?

A NonLinLoc model is **not an ASCII point cloud**. Each pair has:

- A `.hdr` text file with node counts, x/y/z origin, grid spacing, value type,
  numeric precision and geographic projection. TIME headers also give the
  station's name and coordinates.
- A `.buf` binary array, with depth changing fastest, then y, then x. TIME
  values belong to nodes. Standard NLL velocity values belong to cell centers,
  with an unused final plane on each axis; see [Layers and grid resolution](#layers-and-grid-resolution).

The regular grid uses projected Cartesian x/y coordinates and positive-down
z, in kilometres. It does not store latitude and longitude beside every value.
The geographic origin and rotation in a SIMPLE header tell GraphSplit how to
place stations and earthquakes on that grid. For example:

```text
25 23 17 -3.0 -2.75 -0.5 0.25 0.25 0.25 VELOCITY FLOAT
TRANSFORM SIMPLE LatOrig 64.0 LongOrig -20.0 RotCW 27.0
```

This describes a 25 × 23 × 17 grid with 250 m spacing, starting at x = −3 km,
y = −2.75 km and depth = −0.5 km in the rotated model frame. The values are
P or S velocities in km/s. Imported TIME values are in seconds. P and S
velocity grids must have the same spatial definition; all TIME grids must also
agree on that definition.

## A plain-text latitude/longitude velocity model

Use [the ASCII template](../config/graphsplit_3d_ascii.toml). Supply one file per
phase, with exactly these four columns and explicit units:

```text
# latitude_deg longitude_deg depth_km velocity_km_s
64.000 -20.000 -0.5 4.8
64.000 -20.000  0.5 5.1
64.000 -19.990 -0.5 4.9
64.000 -19.990  0.5 5.2
64.010 -20.000 -0.5 4.8
64.010 -20.000  0.5 5.1
64.010 -19.990 -0.5 4.9
64.010 -19.990  0.5 5.2
```

This tiny example includes all eight combinations of two latitudes, two
longitudes and two depths. Real models will have many more nodes. Blank lines
and `#` comments are allowed; row order does not matter. Set:

```toml
[travel_time]
type = "3d"

[grid3d]
model_format = "ascii_velocity"
coordinate_system = "geographic"
vp_file = "vp.txt"
vs_file = "vs.txt"

[coordinates]
reference = "manual"
reference_latitude = 64.0
reference_longitude = -20.0
```

Choose an origin near the study area and keep it fixed. Geographic rows are
converted to GraphSplit's local east/north frame using that origin. This is
**not** an implicit NonLinLoc SIMPLE conversion: the two formulas have different
longitude scaling. A catalog and stations in geographic coordinates are
converted through the same GraphSplit frame. Run metadata records the reference
and input coordinate convention.

For text rows already in local coordinates, use `coordinate_system = "local"`
and columns `x_km y_km depth_km velocity_km_s`. The manual reference must describe
those axes. Depth is always positive down relative to the station-elevation
datum, not depth below the local land surface. `bounds_km` and a surface file
remain in model x/y coordinates even when the velocity input is geographic.

The initial ASCII reader requires a **complete regular grid**, with at least
two nodes per axis. It rejects missing nodes, duplicates, irregular spacing and
nonpositive velocities. A genuinely scattered point cloud needs a separately
chosen interpolation method and a coverage rule. Automatically choosing nearest
neighbors or inverse-distance averaging could fill unsupported areas or smear a
velocity interface. Grid those data explicitly first; this version does not
make that scientific choice silently.

The file is scanned for its axes before allocating a dense velocity array.
The ordinary 3D memory guard then includes the native ASCII model and the
resampled forward grid. It does not retain a second full list of text rows.

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
conversion, including its 6371.0087714 km sphere radius, latitude-dependent
longitude scale and rotation. This matches current NonLinLoc; files built with
older/custom geographic constants need to be regenerated or converted explicitly.
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

NonLinLoc's `Vel2Grid` and `Vel2Grid3D` place velocity values at **cell
centers**, half a spacing beyond the header origin. Their last array plane on
each axis is unused padding. Travel times instead belong to **grid nodes**.
The default `model_sampling = "cell_centers"` follows this convention and
ignores the padding. If another exporter writes velocities at the header's
node coordinates, set `model_sampling = "nodes"` explicitly. The NLL header
does not record this distinction. ASCII velocities always belong to their
listed coordinates, and imported TIME grids are always nodal.

`model_interpolation = "nearest"` transfers the closest velocity sample to
the FMM grid without blending across a layer boundary. At an exact cell
boundary, the cell in the positive coordinate direction is selected.
`"slowness_linear"` interpolates reciprocal velocity and is suitable for a
smooth model. For cell-centered input, the boundary half-cell retains the
nearest interior value. Neither method extends outside the model box or uses
the padding. Neither preserves a thin layer that the FMM grid fails to sample.
A discontinuity may move by part of a cell when grids differ, and the two
solvers have different discretization errors. Import existing TIME grids if
the intention is to retain travel times calculated by NonLinLoc.

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
16 × surface-file bytes`. For ASCII input, the input-buffer allowance is
16 bytes per native model node, including parsing and occupancy overhead. The extra terms allow for serial FMM work, model
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


## Import and geographic-reference checks

The tests include small P/S velocity and TIME grids generated by NonLinLoc's
own `Vel2Grid` and `Grid2Time`, using a SIMPLE origin of 64°N, 20°W and a 27°
rotation. NonLinLoc itself supplies independent projected coordinates and
interpolated travel times for five geographic test points. GraphSplit imports
both file types while deliberately using a different internal geographic origin.
The TIME-import coordinates and interpolated values are checked against those
independent queries, as are GraphSplit's spatial derivatives. Station-coordinate
mismatches are rejected.

This verifies actual NonLinLoc files and the supported projection, rather than
only files written by GraphSplit's test helpers. It does not establish support
for the other NonLinLoc projections or constitute an end-to-end GrowClust3D run.
The [fixture notes](../test/fixtures/nll_simple/README.md) identify the source
version and reproduction commands. Separate tests cover regular ASCII models in
both geographic and local coordinates.
