"""Optional illustration of the Julia example; requires NumPy and Matplotlib."""
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

here = Path(__file__).resolve().parent
source = here / "generated"
a = np.genfromtxt(source / "one_dimensional_errors.csv", delimiter=",", names=True)
b = np.genfromtxt(source / "three_dimensional_errors.csv", delimiter=",", names=True)
plt.rcParams.update({"font.size": 10, "axes.spines.top": False, "axes.spines.right": False})
fig, axes = plt.subplots(2, 3, figsize=(10.4, 6.6), constrained_layout=True)
for col, (label, data, prefix) in enumerate((("Starting locations", a, "initial"), ("1D relocation", a, "relocated"), ("True 3D model", b, "relocated"))):
    for row, (horizontal, horizontal_label) in enumerate((("x", "East (km)"), ("y", "North (km)"))):
        ax = axes[row, col]
        tx = data[f"truth_{horizontal}_m"] / 1000
        tz = data["truth_z_m"] / 1000
        px = data[f"{prefix}_{horizontal}_m"] / 1000
        pz = data[f"{prefix}_z_m"] / 1000
        for xi, zi, xr, zr in zip(tx, tz, px, pz):
            ax.plot([xi, xr], [zi, zr], color="#aaaaaa", lw=.55, zorder=1)
        ax.scatter(tx, tz, marker="+", color="#222222", s=21, lw=.8, label="True location", zorder=3)
        ax.scatter(px, pz, color="#2076a3", s=11, alpha=.85, label="Estimated location", zorder=2)
        ax.set(xlim=(-3.6, 3.6), ylim=(6.5, .5), xlabel=horizontal_label)
        ax.set_aspect("equal", adjustable="box")
        ax.grid(alpha=.14)
        if col == 0:
            ax.set_ylabel("Depth (km)")
        else:
            ax.tick_params(labelleft=False)
        if row == 0:
            ax.set_title(label)
axes[0, 0].legend(loc="lower left", frameon=True, framealpha=.95, fontsize=8)
fig.savefig(here / "comparison.png", dpi=180)
