# Yu–Ellsworth–Beroza benchmark utility

This directory is separate from GraphSplit. It evaluates a relocated standard
catalog against the supplied `truelocs.txt` and writes two CSV files.

```bash
julia compare_catalogs.jl relocated_catalog.txt input/truelocs.txt [output_prefix]
```

The implementation follows the authors' public `util/error_analysis.py`:

- horizontal accuracy is WGS84 geodesic distance from each event to its truth;
- depth accuracy is absolute depth difference;
- truth neighbors have both horizontal and vertical separation below 2 km;
- each event's precision is the RMS misfit of its neighbor distances, followed
  by a mean across events;
- Chamfer distance follows Point Cloud Utils: half the sum of the two mean
  nearest-neighbor Euclidean distances after WGS84 geodetic-to-ENU conversion
  about the same reference point used by the authors.

Serial catalog ID `k` selects truth row `k`, allowing filtered GraphSplit
catalogs to be evaluated without renumbering.

As a smoke-test, comparing the supplied seed `input/catalog.txt` to its truth
should report approximately 18,200 unique neighbor pairs, mean accuracy
`0.917479 / 1.124088 km` (horizontal/depth), mean precision
`0.551643 / 0.711940 km`, and Chamfer distance `0.856395 km`. Small last-digit
differences between WGS84 inverse implementations are acceptable.

Reference: Y. Yu, W. L. Ellsworth, and G. C. Beroza, “Accuracy and Precision of
Earthquake Location Programs: Insights from a Synthetic Controlled
Experiment,” *Seismological Research Letters* 96, 1860–1874 (2025),
doi:10.1785/0220240354. The benchmark code and data cited by the paper are at
<https://github.com/YuYifan2000/comparison_hypoDD_GrowClust/>.
