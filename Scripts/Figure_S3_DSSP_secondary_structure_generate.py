#!/usr/bin/env python3
"""
Figure S3: DSSP secondary structure heatmaps for all systems.
Uses config.py.  No overall title; system names retained as subplot labels.
Exports: figure (pdf, png, jpg), combined CSV, per-system Excel sheets.
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
OUTPUT_DIR = OUT_DIR / "Figure_S3_DSSP_heatmap"
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

# ------------------- Parser function -------------------
def parse_ss_file(fpath):
    """
    Read DSSP file, return numpy array of shape (n_frames, n_residues)
    with integer codes for each residue.
    """
    with open(fpath) as f:
        lines = [line.strip() for line in f if line.strip() and not line.startswith('#')]
    if not lines:
        return None
    max_len = max(len(line) for line in lines)
    padded = [line + 'C' * (max_len - len(line)) for line in lines]
    encoded = [[CODE_MAP.get(ch, 3) for ch in line] for line in padded]
    return np.array(encoded, dtype=int)

# ------------------- Plotting -------------------
# Use subplots without constrained_layout to avoid warning with subplots_adjust
fig, axes = plt.subplots(2, 2, figsize=(16, 10))
axes = axes.flatten()
plt.subplots_adjust(wspace=0.25, hspace=0.25, left=0.06, right=0.98, top=0.94, bottom=0.06)

all_data = []
im = None
labels = ['A', 'B', 'C', 'D']

for idx, (ax, sys) in enumerate(zip(axes.flat, SYSTEMS)):
    ax.text(-0.12, 1.06, labels[idx], transform=ax.transAxes,
            fontweight='bold', fontsize=9, va='bottom', ha='center', color='black')
    ax.set_title(sys, fontweight='normal', fontsize=9, color=COLORS[sys])

    f = DATA_DIR / sys / 'ss.dat'
    if not f.exists():
        ax.text(0.5, 0.5, 'No data', ha='center', va='center')
        continue

    data = parse_ss_file(f)
    if data is None or data.size == 0:
        ax.text(0.5, 0.5, 'Empty', ha='center', va='center')
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

    if n_res > 100:
        step = max(1, n_res // 20)
        ax.set_yticks(np.arange(1, n_res + 1, step))

    # Collect data for export
    for i, t in enumerate(time_ns):
        for j, code in enumerate(data[i]):
            all_data.append([sys, t, j + 1, code])

# Add single colorbar for all subplots
if im is not None:
    cbar = fig.colorbar(im, ax=axes, location='right', shrink=1.0, pad=0.02)
    cbar.set_ticks([0.5, 1.5, 2.5, 3.5])
    cbar.set_ticklabels(SS_NAMES)
    cbar.set_label('Secondary structure')

# ---------- Save figure in all three formats ----------
for ext in ('png', 'jpg', 'pdf'):      # <-- 修改此处：同时输出三种格式
    out_path = OUTPUT_DIR / f'Figure_S3_DSSP_heatmap.{ext}'
    plt.savefig(out_path, dpi=DPI, facecolor='white', bbox_inches='tight')
    print(f"Saved: {out_path}")

plt.close()

# ---------- Export data ----------
if all_data:
    df = pd.DataFrame(all_data, columns=['System', 'Time_ns', 'Residue', 'SS_code'])
    csv_path = OUTPUT_DIR / 'Figure_S3_DSSP.csv'
    df.to_csv(csv_path, index=False)
    print(f"CSV saved to {csv_path}")

    xlsx_path = OUTPUT_DIR / 'Figure_S3_DSSP_data.xlsx'
    try:
        with pd.ExcelWriter(xlsx_path, engine='openpyxl') as writer:
            for sys in SYSTEMS:
                sys_df = df[df['System'] == sys]
                if not sys_df.empty:
                    sys_df.to_excel(writer, sheet_name=sys, index=False)
        print(f"Excel saved to {xlsx_path}")
    except Exception as e:
        print(f"Excel export failed: {e}")

print("Figure_S3 generation complete.")