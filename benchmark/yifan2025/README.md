# Yu–Ellsworth–Beroza benchmark utility

This standalone comparator evaluates a relocated catalog against externally
obtained true locations. It does not run GraphSplit or DDSync.

## Data source and citation

The example files in `input/` are derived from the synthetic benchmark of:

Yu, Y., Ellsworth, W. L., and Beroza, G. C. (2025). *Accuracy and Precision of
Earthquake Location Programs: Insights from a Synthetic Controlled Experiment*.
**Seismological Research Letters, 96**(3), 1860–1874.
[https://doi.org/10.1785/0220240354](https://doi.org/10.1785/0220240354).

Please cite this paper when using the example data. The authors' code and
source data are at
[YuYifan2000/comparison_hypoDD_GrowClust](https://github.com/YuYifan2000/comparison_hypoDD_GrowClust).
The included `catalog.txt`, `stations.txt` and `vm.txt` provide the starting
catalog, station geometry and layered velocity model. Their formatting follows
[GraphSplit's input conventions](../../docs/FILE_FORMATS.md).

`truelocs.txt` is **not included**. Follow the source repository's instructions
to obtain the matching true locations from the authors. The expected file has
latitude, longitude and depth in kilometers, with benchmark event ID `k` in
row `k`. These example inputs also omit synchronized theta files; running a
relocation additionally requires DDSync processing products.

## Run the comparison

From this directory:

```bash
julia compare_catalogs.jl relocated_catalog.txt /path/to/truelocs.txt [output_prefix]
```

The command writes `<output_prefix>_metrics.csv` and
`<output_prefix>_event_errors.csv`. Catalog IDs select truth rows, allowing
filtered catalogs to be evaluated without renumbering. Compare methods on the
same event population; filtering changes the evaluated neighbor pairs.

## Metric definitions

The calculation follows the definitions in the authors' public
`util/error_analysis.py`, with the Chamfer normalization noted below:

- Horizontal accuracy is the WGS84 geodesic distance from each event to its
  true location; depth accuracy is the absolute depth difference.
- Truth neighbors have both horizontal and depth separation below 2 km.
- Each event's precision is the RMS error in its neighbor separations, followed
  by a mean across events.
- This utility reports **half the sum** of the two directed mean nearest-neighbor
  distances after WGS84 geodetic-to-ENU conversion. Yu's published Python
  evaluation and the GraphSplit manuscript use the **sum**. Multiply this
  utility's `chamfer_distance_km` by two for that convention; the accuracy and
  precision values need no conversion.

For the matching 1,000-event realization, comparing `input/catalog.txt` to the
external truth should give approximately 18,200 unique neighbor pairs, mean
accuracy `0.917479 / 1.124088 km` (horizontal/depth), mean precision
`0.551643 / 0.711940 km`, and half-sum Chamfer `0.856395 km`.
Small last-digit differences between WGS84 inverse implementations are
acceptable. The repository's automated smoke test uses generated test data and
does not require the external true locations.
