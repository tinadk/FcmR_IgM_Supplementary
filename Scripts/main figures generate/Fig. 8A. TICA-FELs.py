#!/usr/bin/env python3
"""
Fig. 8A – TICA-FEL overlay (colored by system) and 2x2 free energy landscapes.
2x2 colormap: plasma (purple‑yellow, no black), contour lines in white (thin, semi‑transparent).
Exports: overlay figure, 2x2 figure, combined CSV, combined XLSX (single sheet).
Uses config.py for settings.
"""

import sys
from pathlib import Path
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from config import SYSTEMS, COLORS, OUT_DIR, DPI, MPL_RCPARAMS

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import pandas as pd
import warnings

warnings.filterwarnings('ignore')
plt.rcParams.update(MPL_RCPARAMS)

TICA_DIR = ROOT / 'MSM_final' / 'archive' / 'complete' / 'data' / 'tica'
OUTPUT_DIR = OUT_DIR / 'fig8a_tica_fels_output'
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

kT = 2.494
BINS = 80
COLOR_LIST = [COLORS[s] for s in SYSTEMS]

def load_tic_data():
    data = {}
    for sys in SYSTEMS:
        fpath = TICA_DIR / f'{sys}_tica.npy'
        if fpath.exists():
            data[sys] = np.load(fpath)[:, :2]
        else:
            print(f'  WARNING: {fpath} not found, skipping {sys}.')
            data[sys] = None
    return data

DATA = load_tic_data()

# ---------- Combine all data ----------
all_rows = []
for sys, arr in DATA.items():
    if arr is not None:
        for row in arr:
            all_rows.append([sys, row[0], row[1]])

if all_rows:
    combined_csv = OUTPUT_DIR / 'Fig8A_overlay_data.csv'
    with open(combined_csv, 'w') as f:
        f.write('System,TIC1,TIC2\n')
        for r in all_rows:
            f.write(f'{r[0]},{r[1]:.6f},{r[2]:.6f}\n')
    print(f'  Saved {combined_csv}')

    df_all = pd.DataFrame(all_rows, columns=['System', 'TIC1', 'TIC2'])
    xlsx_path = OUTPUT_DIR / 'Fig8A_FEL_data.xlsx'
    df_all.to_excel(xlsx_path, sheet_name='All_Systems', index=False)
    print(f'  Saved {xlsx_path}')
else:
    print('  No data to export.')

# ===================== Overlay (colored by system) =====================
print('1. Plotting TICA scatter overlay...')
fig, ax = plt.subplots(figsize=(10, 8))
for sys, col in zip(SYSTEMS, COLOR_LIST):
    arr = DATA.get(sys)
    if arr is None:
        continue
    ax.scatter(arr[:, 0], arr[:, 1], c=col, s=4, alpha=0.4, label=sys, rasterized=True)

ax.set_xlabel('TIC 1')
ax.set_ylabel('TIC 2')
ax.set_title('TICA projection)')
ax.legend(markerscale=4, loc='lower right', ncol=1, fontsize=8,
          framealpha=0.7, edgecolor='none', labelspacing=0.5, handletextpad=0.4)
plt.tight_layout()
for ext in ('png', 'jpg', 'pdf'):
    plt.savefig(OUTPUT_DIR / f'Fig8A_TICA_overlay.{ext}', dpi=DPI,
                bbox_inches='tight', facecolor='white')
plt.close()
print('  Overlay saved.')

# ===================== 2x2 FEL (plasma + white contours) =====================
print('2. Plotting TICA-FEL 2x2 (plasma, white contours)...')
fig, axes = plt.subplots(2, 2, figsize=(12, 10))
for ax, sys in zip(axes.flat, SYSTEMS):
    arr = DATA.get(sys)
    if arr is None:
        ax.text(0.5, 0.5, 'No data', transform=ax.transAxes, ha='center')
        continue
    H, xe, ye = np.histogram2d(arr[:, 0], arr[:, 1], bins=BINS, density=True)
    F = -kT * np.log(H.T + 1e-12)
    F -= F.min()
    im = ax.imshow(F, origin='lower', extent=[xe[0], xe[-1], ye[0], ye[-1]],
                   aspect='auto', cmap='plasma',
                   vmin=0, vmax=np.percentile(F, 95))
    levels = np.arange(0, F.max(), 2)
    ax.contour(F, levels, extent=[xe[0], xe[-1], ye[0], ye[-1]],
               origin='lower', colors='white', linewidths=0.4, alpha=0.4)
    plt.colorbar(im, ax=ax, label='Free energy (kJ/mol)')
    ax.set_title(sys, fontweight='bold')
    ax.set_xlabel('TIC 1')
    ax.set_ylabel('TIC 2')
plt.suptitle('TICA free energy landscapes', fontweight='bold', fontsize=14)
plt.tight_layout()
for ext in ('png', 'jpg', 'pdf'):
    plt.savefig(OUTPUT_DIR / f'Fig8A_TICA_FEL_2x2.{ext}', dpi=DPI,
                bbox_inches='tight', facecolor='white')
plt.close()
print('  2x2 saved.')

print(f'All Fig.8A outputs in {OUTPUT_DIR}')