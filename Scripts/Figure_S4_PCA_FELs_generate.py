#!/usr/bin/env python3
"""
Figure S4 – PCA-FEL 2x2 free energy landscapes.
Reads PCA data from DATA_DIR/<system>/pca_2d.xvg.
Plasma colormap with white contour lines.
Shared colorbar placed at the right of the figure.
Exports: figure (png/jpg/pdf), combined CSV, combined XLSX.
"""

import sys
from pathlib import Path
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from config import SYSTEMS, COLORS, OUT_DIR, DPI, MPL_RCPARAMS, DATA_DIR

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import pandas as pd
import warnings

def load_pca_data():
    data = {}
    for sys in SYSTEMS:
        f = DATA_DIR / sys / 'pca_2d.xvg'
        if f.exists():
            arr = np.loadtxt(f, comments=['@', '#'])
            data[sys] = arr
        else:
            print(f'[WARNING] {sys}: pca_2d.xvg not found')
    return data

warnings.filterwarnings('ignore')
plt.rcParams.update(MPL_RCPARAMS)

OUTPUT_DIR = OUT_DIR / 'Figure_S4_PCA_output'
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

kT = 2.494
BINS = 80
LABEL_FONTSIZE = 9
PANEL_LABEL_POS = (-0.12, 1.06)

DATA = load_pca_data()

all_rows = []
for sys, arr in DATA.items():
    if arr is not None:
        for row in arr:
            all_rows.append([sys, row[0], row[1]])

if all_rows:
    # Combined CSV (unchanged)
    combined_csv = OUTPUT_DIR / 'Figure_S4_PCA_FEL_2x2.csv'
    with open(combined_csv, 'w') as f:
        f.write('System,PC1,PC2\n')
        for r in all_rows:
            f.write(f'{r[0]},{r[1]:.6f},{r[2]:.6f}\n')
    print(f'  Saved {combined_csv}')

    # Build full DataFrame
    df_all = pd.DataFrame(all_rows, columns=['System', 'PC1', 'PC2'])

    # Excel: only per-system sheets, no summary sheet
    xlsx_path = OUTPUT_DIR / 'Figure_S4_PCA_FEL_data.xlsx'
    with pd.ExcelWriter(xlsx_path, engine='openpyxl') as writer:
        for sys in SYSTEMS:
            sys_df = df_all[df_all['System'] == sys]
            if not sys_df.empty:
                sys_df.to_excel(writer, sheet_name=sys, index=False)
    print(f'  Saved {xlsx_path}')
else:
    print('  No data to export.')

fig, axes = plt.subplots(2, 2, figsize=(14, 12))

plt.subplots_adjust(wspace=0.25, hspace=0.25, left=0.06, right=0.98, top=0.94, bottom=0.06)

labels = ['A', 'B', 'C', 'D']

for idx, (ax, sys) in enumerate(zip(axes.flat, SYSTEMS)):
    ax.text(PANEL_LABEL_POS[0], PANEL_LABEL_POS[1], labels[idx],
            transform=ax.transAxes, fontweight='bold', fontsize=LABEL_FONTSIZE,
            va='bottom', ha='center', color='black')

    arr = DATA.get(sys)
    if arr is None:
        ax.text(0.5, 0.5, 'No data', transform=ax.transAxes, ha='center', va='center')
        continue

    H, xe, ye = np.histogram2d(arr[:,0], arr[:,1], bins=BINS, density=True)
    F = -kT * np.log(H.T + 1e-12)
    F -= F.min()
    im = ax.imshow(F, origin='lower', extent=[xe[0], xe[-1], ye[0], ye[-1]],
                   aspect='auto', cmap='plasma',
                   vmin=0, vmax=np.percentile(F, 95))
    levels = np.arange(0, F.max(), 2)
    ax.contour(F, levels, extent=[xe[0], xe[-1], ye[0], ye[-1]],
               origin='lower', colors='white', linewidths=0.4, alpha=0.4)
    ax.set_title(sys, fontweight='normal', fontsize=LABEL_FONTSIZE, color=COLORS[sys])
    ax.set_xlabel('PC1')
    ax.set_ylabel('PC2')

cbar = fig.colorbar(im, ax=axes.ravel().tolist(), location='right', shrink=1.0, pad=0.02)
cbar.set_label('Free energy (kJ/mol)', fontsize=LABEL_FONTSIZE)


for ext in ('png', 'jpg', 'pdf'):
    fig.savefig(OUTPUT_DIR / f'Figure_S4_PCA_FEL_2x2.{ext}', dpi=DPI,
                bbox_inches='tight', facecolor='white')
plt.close()
print('  2x2 saved with shared colorbar.')
print(f'Figure_S4 outputs saved to {OUTPUT_DIR}')