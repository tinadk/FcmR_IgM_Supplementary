#!/usr/bin/env python3
"""
Fig. 7A: RMSD, Rg, and SASA time series (combined 1×3 panel).
All parameters pulled from config.py.  Data expected as:
    data/{system}/rmsd.xvg
    data/{system}/rg.xvg
    data/{system}/sasa.xvg
or as specified in the script's DATA_PATH section.
"""

import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from config import (SYSTEMS, COLORS, DATA_DIR, OUT_DIR,
                    DPI, MPL_RCPARAMS, TIME_MAX_NS, EQUILIBRATION_NS)

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

# Apply global plot style
plt.rcParams.update(MPL_RCPARAMS)

# ------------------- Output directory -------------------
OUTPUT_DIR = OUT_DIR / "fig7a_output"
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

# ---------- Per‑figure customization ----------
LINE_WIDTH = 1.8
ALPHA = 0.9
SMOOTH_WINDOW = 50          # rolling window for visual smoothing


# ---------- Helper: load and smooth ----------

# ---------- Plotting ----------
fig, axes = plt.subplots(1, 3, figsize=(24, 6))

plot_configs = [
    ('rmsd.xvg', 'RMSD (nm)', (0, None)),
    ('rg.xvg', 'Rg (nm)', (0, None)),
    ('sasa.xvg', 'SASA (nm²)', (110, None))
]

for ax, (fname, ylabel, title, ylim) in zip(axes, plot_configs):
    for sys in SYSTEMS:
        c = COLORS[sys]
        t, y = load_and_smooth(sys, fname)
        if t is None:
            continue
        ax.plot(t, y, color=c, linewidth=LINE_WIDTH, alpha=ALPHA, label=sys)
    ax.axvline(EQUILIBRATION_NS, color='gray', linestyle='--',
               alpha=0.5, linewidth=1)
    ax.set(xlabel='Time (ns)', ylabel=ylabel, title=title)
    ax.legend(loc='lower right')
    if ylim[0] is not None or ylim[1] is not None:
        ax.set_ylim(*ylim)
    # For SASA, expand y top margin slightly
    if 'SASA' in title:
        ymin, ymax = ax.get_ylim()
        ax.set_ylim(ymin, ymax * 1.05)

plt.tight_layout()

# Save in three formats
basename = 'S30_rmsd_rg_sasa_combine'
for fmt in ('png', 'jpg', 'pdf'):
    out_path = OUT_DIR / f"{basename}.{fmt}"
    fig.savefig(out_path, dpi=DPI, bbox_inches='tight', facecolor='white')
    print(f"Saved: {out_path}")

plt.close()
print("Combined RMSD/Rg/SASA figure saved.")