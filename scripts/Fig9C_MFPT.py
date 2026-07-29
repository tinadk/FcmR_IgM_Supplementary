#!/usr/bin/env python3
"""
Fig. 9C – Mean First Passage Time (MFPT) matrices for all four systems.
All raw MFPT values shown directly in the heatmap (log‑scale colormap 'plasma',
NO white: bad values set to colormap minimum).
Exports:
  - 2x2 figure (log‑scale plasma, exact values labeled)
  - Combined CSV (all systems, raw values)
  - Combined XLSX (four sheets, raw values)
No per‑system CSVs. Uses config.py for settings.
"""

import sys
from pathlib import Path
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from config import SYSTEMS, OUT_DIR, DPI, MPL_RCPARAMS

import numpy as np
import matplotlib.pyplot as plt
import matplotlib
from matplotlib.colors import LogNorm
import pandas as pd
import warnings

warnings.filterwarnings('ignore')
plt.rcParams.update(MPL_RCPARAMS)

# ---------- Parameters ----------
LAG_TIME = 10          # frames per step
FRAME_TO_NS = 0.1      # ns per frame

# ---------- Paths ----------
MSM_DIR = ROOT / 'MSM_final' / 'archive' / 'complete' / 'data' / 'msm'
TICA_DIR = ROOT / 'MSM_final' / 'archive' / 'complete' / 'data' / 'tica'

OUTPUT_DIR = OUT_DIR / 'fig9c_mfpt'
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

# ---------- Load common data ----------
pcca_membership = np.load(MSM_DIR / 'pcca_membership.npy')
dtraj_all = np.load(MSM_DIR / 'dtraj_k150.npy')
system_labels = np.load(TICA_DIR / 'system_labels.npy')
macro_assign = pcca_membership[dtraj_all]
n_macro = pcca_membership.max() + 1
macro_labels = [f'M{i+1}' for i in range(n_macro)]
sys_to_id = {name: i for i, name in enumerate(SYSTEMS)}

# ---------- Helper functions ----------

# ---------- Compute MFPT ----------
mfpt = {}
for sys_name in SYSTEMS:
    mfpt[sys_name] = mfpt_for_system(sys_to_id[sys_name])

# ---------- Plot 2x2 with plasma colormap (bad values = colormap minimum) ----------
fig, axes = plt.subplots(2, 2, figsize=(14, 12))
axes = axes.flatten()

all_vals = np.concatenate([m[m > 1e-6] for m in mfpt.values()])
vmin = max(1e-3, np.nanmin(all_vals))
vmax = np.nanmax(all_vals)

cmap = plt.cm.plasma.copy()
cmap.set_bad(cmap(0))   # set NaN/masked cells to the darkest color (deep purple)

for ax, sys_name in zip(axes, SYSTEMS):
    data = mfpt[sys_name].copy()
    masked = np.ma.array(data, mask=np.isnan(data))

    im = ax.imshow(masked, cmap=cmap, aspect='auto', origin='lower',
                   norm=LogNorm(vmin=vmin, vmax=vmax))

    cbar = plt.colorbar(im, ax=ax, label='MFPT (ns)', fraction=0.046, pad=0.04)

    # Label each cell with its original value (white text)
    for i in range(n_macro):
        for j in range(n_macro):
            val = data[i, j]
            if i == j or np.isnan(val):
                continue
            if val >= 1e6:
                txt = f'{val:.1e}'
            else:
                txt = f'{val:.0f}'
            ax.text(j, i, txt, ha='center', va='center', fontsize=6.5,
                    color='white', fontweight='bold')

    ax.set_xticks(range(n_macro))
    ax.set_xticklabels(macro_labels)
    ax.set_yticks(range(n_macro))
    ax.set_yticklabels(macro_labels)
    ax.set_xlabel('Target macrostate')
    ax.set_ylabel('Source macrostate')
    ax.set_title(sys_name)

plt.tight_layout()

for ext in ('png', 'jpg', 'pdf'):
    plt.savefig(OUTPUT_DIR / f'Fig9C_MFPT_matrices.{ext}',
                dpi=DPI, bbox_inches='tight', facecolor='white')
plt.close()

# ---------- Export combined CSV (all systems, raw values) ----------
combined_rows = []
for sys_name in SYSTEMS:
    mat = mfpt[sys_name]
    for i, source in enumerate(macro_labels):
        for j, target in enumerate(macro_labels):
            if i != j and not np.isnan(mat[i, j]):
                combined_rows.append([sys_name, source, target, mat[i, j]])

df_combined = pd.DataFrame(combined_rows, columns=['System', 'Source', 'Target', 'MFPT_ns'])
df_combined.to_csv(OUTPUT_DIR / 'Fig9C_MFPT_all_systems.csv', index=False, float_format='%.6f')
print(f'Combined CSV saved to {OUTPUT_DIR / "Fig9C_MFPT_all_systems.csv"}')

# ---------- Export combined XLSX (one sheet per system, raw values) ----------
xlsx_path = OUTPUT_DIR / 'Fig9C_MFPT_data.xlsx'
with pd.ExcelWriter(xlsx_path, engine='openpyxl') as writer:
    for sys_name in SYSTEMS:
        df = pd.DataFrame(mfpt[sys_name], index=macro_labels, columns=macro_labels)
        df.index.name = 'Source'
        df.columns.name = 'Target'
        df.to_excel(writer, sheet_name=sys_name)
print(f'Combined XLSX saved to {xlsx_path}')

print('All files ready. No per‑system CSVs generated.')