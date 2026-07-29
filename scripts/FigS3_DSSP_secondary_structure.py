#!/usr/bin/env python3
"""
Figure S3: DSSP secondary structure heatmaps for all systems.
Uses config.py.  No overall title; system names retained as subplot labels.
Exports: figure (png, jpg, pdf), combined CSV, per‑system Excel sheets.
"""

import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import ListedColormap, BoundaryNorm

from config import DATA_DIR, OUT_DIR, SYSTEMS, COLORS, DPI, MPL_RCPARAMS

plt.rcParams.update(MPL_RCPARAMS)

# ------------------- Output directory -------------------
OUTPUT_DIR = OUT_DIR / "FigS3_DSSP_heatmap"
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

# ------------------- Constants -------------------
FRAME_PS = 10.0
SS_NAMES = ['Helix', 'Strand', 'Turn', 'Coil']
SS_COLORS = ['#1f77b4', '#d62728', '#ff7f0e', '#7f7f7f']

CODE_MAP = {
    'H': 0, 'G': 0, 'I': 0,
    'E': 1, 'B': 1,
    'T': 2, 'S': 2,
    'C': 3, ' ': 3
}

# ------------------- Plotting -------------------
fig, axes = plt.subplots(2, 2, figsize=(16, 10), constrained_layout=True)
axes = axes.flatten()
all_data = []
im = None

for ax, sys in zip(axes, SYSTEMS):
    f = DATA_DIR / sys / 'ss.dat'
    if not f.exists():
        ax.text(0.5, 0.5, 'No data', ha='center', va='center')
        ax.set_title(sys, fontweight='bold', color=COLORS[sys])
        continue

    data = parse_ss_file(f)
    if data is None or data.size == 0:
        ax.text(0.5, 0.5, 'Empty', ha='center', va='center')
        ax.set_title(sys, fontweight='bold', color=COLORS[sys])
        continue

    n_frames, n_res = data.shape
    time_ns = np.arange(n_frames) * FRAME_PS / 1000.0
    residues = np.arange(1, n_res + 1)

    cmap = ListedColormap(SS_COLORS)
    bounds = [0, 1, 2, 3, 4]
    norm = BoundaryNorm(bounds, cmap.N)

    extent = [time_ns[0], time_ns[-1], residues[0], residues[-1]]
    im = ax.imshow(data.T, origin='lower', cmap=cmap, norm=norm,
                   aspect='auto', extent=extent, interpolation='none')
    ax.set_xlabel('Time (ns)')
    ax.set_ylabel('Residue number')
    # Keep system name as subplot identifier (no DSSP keyword here)
    ax.set_title(sys, fontweight='bold', color=COLORS[sys])

    if n_res > 100:
        step = max(1, n_res // 20)
        ax.set_yticks(np.arange(1, n_res + 1, step))

    for i, t in enumerate(time_ns):
        for j, code in enumerate(data[i]):
            all_data.append([sys, t, j + 1, code])

# Add single colorbar for all subplots
if im is not None:
    cbar = fig.colorbar(im, ax=axes, location='right', shrink=0.7, pad=0.02)
    cbar.set_ticks([0.5, 1.5, 2.5, 3.5])
    cbar.set_ticklabels(SS_NAMES)
    cbar.set_label('Secondary structure')

# No overall suptitle (removed)

# ---------- Save figure ----------
for ext in ('pdf', 'jpg', 'png'):
    out_path = OUTPUT_DIR / f'FigS3_DSSP_heatmap.{ext}'
    plt.savefig(out_path, dpi=DPI, facecolor='white')
    print(f"Saved: {out_path}")
plt.close()

# ---------- Export data ----------
if all_data:
    df = pd.DataFrame(all_data, columns=['System', 'Time_ns', 'Residue', 'SS_code'])
    csv_path = OUTPUT_DIR / 'FigS3_DSSP.csv'
    df.to_csv(csv_path, index=False)
    print(f"CSV saved to {csv_path}")

    xlsx_path = OUTPUT_DIR / 'FigS3_DSSP_data.xlsx'
    try:
        with pd.ExcelWriter(xlsx_path, engine='openpyxl') as writer:
            for sys in SYSTEMS:
                sys_df = df[df['System'] == sys]
                if not sys_df.empty:
                    sys_df.to_excel(writer, sheet_name=sys, index=False)
        print(f"Excel saved to {xlsx_path}")
    except Exception as e:
        print(f"Excel export failed: {e}")

print("FigS3 generation complete.")