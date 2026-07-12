#!/usr/bin/env python3
"""
Extracted plotting code with correct TICA path.
1. TICA-FEL scatter overlay (dot size=4)
2. TICA-FEL 2x2 per system
3. FEL scatter overlay (using TICA data, generic labels)
4. FEL 2x2 (using TICA data, generic labels)
(Removed PCA scatter overlay and PCA-FEL 2x2 as requested)
Outputs: CSV files per system, combined CSV for each figure, Arial font, DPI 1600.
"""

import numpy as np
import matplotlib
matplotlib.use('Agg')          # for headless environments
import matplotlib.pyplot as plt
from pathlib import Path
import warnings

warnings.filterwarnings('ignore')

# ================= CONFIGURATION =================
TICA_DIR = Path('MSM_rebuild_final_part3/LEGACY_ARCHIVE/ARCHIVE_COMPLETE/data/tica')
PCA_DIR  = Path('data')        # kept for possible future use, not used here
OUT_DIR  = Path('use_MSM_publication_figures')
OUT_DIR.mkdir(exist_ok=True)

# DPI settings: all 1600 for publication quality
SAVE_DPI = {'png': 1600, 'jpg': 1600, 'pdf': 1600}

# Set Arial as default font
matplotlib.rcParams['font.family'] = 'Arial'

kT = 2.494
SYSTEMS = ['WT', 'D111A', 'T110A_D111A', 'T60A_S63A']
COLORS  = {'WT': '#2C3E50', 'D111A': '#E74C3C', 'T110A_D111A': '#27AE60', 'T60A_S63A': '#2980B9'}
COLOR_LIST = [COLORS[s] for s in SYSTEMS]

def save_fig(fig, name):
    """Save figure in multiple formats with appropriate DPI."""
    print(f'  Saving {name}...')
    for ext, dpi in SAVE_DPI.items():
        fig.savefig(OUT_DIR / f'{name}.{ext}', dpi=dpi,
                    bbox_inches='tight', facecolor='white')
    plt.close(fig)
    print(f'  {name} saved.')

def make_fel(ax, x, y, title, xlab='TIC 1', ylab='TIC 2'):
    """Generate a free energy landscape using histogram2d."""
    H, xe, ye = np.histogram2d(x, y, bins=60, density=True)
    F = -kT * np.log(H.T + 1e-12)
    F -= F.min()
    im = ax.imshow(F, origin='lower', extent=[xe[0], xe[-1], ye[0], ye[-1]],
                   aspect='auto', cmap='coolwarm')
    plt.colorbar(im, ax=ax, label='Free energy (kJ/mol)')
    ax.set_title(title, fontweight='bold')
    ax.set_xlabel(xlab)
    ax.set_ylabel(ylab)

def load_tic_data():
    """Load TIC data for all systems, returns dict {system: (N,2) array} or None if missing."""
    data = {}
    for sys in SYSTEMS:
        fpath = TICA_DIR / f'{sys}_tica.npy'
        if fpath.exists():
            data[sys] = np.load(fpath)[:, :2]
        else:
            print(f'  WARNING: {fpath} not found, skipping {sys}.')
            data[sys] = None
    return data

def save_combined_csv(data_dict, csv_name):
    """Save combined TIC data with a System column."""
    rows = []
    for sys, arr in data_dict.items():
        if arr is not None:
            for row in arr:
                rows.append([sys, row[0], row[1]])
    if rows:
        import csv
        csv_path = OUT_DIR / csv_name
        with open(csv_path, 'w', newline='') as f:
            writer = csv.writer(f)
            writer.writerow(['System', 'TIC1', 'TIC2'])
            writer.writerows(rows)
        print(f'    Saved {csv_path}')
    else:
        print(f'    WARNING: No data for {csv_name}')

# ================= SAVE INDIVIDUAL TICA DATA AS CSV =================
print('0. Saving individual TICA coordinates as CSV...')
for sys in SYSTEMS:
    fpath = TICA_DIR / f'{sys}_tica.npy'
    if not fpath.exists():
        print(f'  WARNING: {fpath} not found, skipping CSV.')
        continue
    X = np.load(fpath)[:, :2]          # first two TICs
    csv_path = OUT_DIR / f'{sys}_tica.csv'
    np.savetxt(csv_path, X, delimiter=',', header='TIC1,TIC2', comments='')
    print(f'    Saved {csv_path}')

# Pre-load all TIC data once for both plotting and combined CSV
DATA = load_tic_data()

# ================= 1. TICA-FEL SCATTER OVERLAY =================
print('1. Plotting TICA-FEL scatter overlay...')
fig, ax = plt.subplots(figsize=(10, 8))
for sys, col in zip(SYSTEMS, COLOR_LIST):
    arr = DATA.get(sys)
    if arr is None:
        continue
    ax.scatter(arr[:, 0], arr[:, 1], c=col, s=4, alpha=0.4, label=sys, rasterized=True)
ax.set_xlabel('TIC 1')
ax.set_ylabel('TIC 2')
ax.set_title('TICA Free Energy Landscape (scatter)')
ax.legend(markerscale=6, loc='lower right')
save_fig(fig, 'TICA_FEL_overlay')
# Combined CSV for this figure
save_combined_csv(DATA, 'TICA_FEL_overlay_data.csv')

# ================= 2. TICA-FEL 2x2 =================
print('2. Plotting TICA-FEL 2x2...')
fig, axes = plt.subplots(2, 2, figsize=(12, 10))
for ax, sys in zip(axes.flat, SYSTEMS):
    arr = DATA.get(sys)
    if arr is None:
        ax.text(0.5, 0.5, 'No data', transform=ax.transAxes, ha='center')
        continue
    make_fel(ax, arr[:, 0], arr[:, 1], sys)
plt.suptitle('TICA Free Energy Landscapes', fontweight='bold', fontsize=14)
plt.tight_layout()
save_fig(fig, 'TICA_FEL_2x2')
save_combined_csv(DATA, 'TICA_FEL_2x2_data.csv')

# ================= 3. GENERIC FEL SCATTER OVERLAY =================
print('3. Plotting FEL scatter overlay (generic)...')
fig, ax = plt.subplots(figsize=(10, 8))
for sys, col in zip(SYSTEMS, COLOR_LIST):
    arr = DATA.get(sys)
    if arr is None:
        continue
    ax.scatter(arr[:, 0], arr[:, 1], c=col, s=4, alpha=0.4, label=sys, rasterized=True)
ax.set_xlabel('Reaction coordinate 1')
ax.set_ylabel('Reaction coordinate 2')
ax.set_title('Free Energy Landscape')
ax.legend(markerscale=6, loc='lower right')
save_fig(fig, 'FEL_overlay')
save_combined_csv(DATA, 'FEL_overlay_data.csv')

# ================= 4. GENERIC FEL 2x2 =================
print('4. Plotting FEL 2x2 (generic)...')
fig, axes = plt.subplots(2, 2, figsize=(12, 10))
for ax, sys in zip(axes.flat, SYSTEMS):
    arr = DATA.get(sys)
    if arr is None:
        ax.text(0.5, 0.5, 'No data', transform=ax.transAxes, ha='center')
        continue
    make_fel(ax, arr[:, 0], arr[:, 1], sys, xlab='RC1', ylab='RC2')
plt.suptitle('Free Energy Landscapes', fontweight='bold', fontsize=14)
plt.tight_layout()
save_fig(fig, 'FEL_2x2')
save_combined_csv(DATA, 'FEL_2x2_data.csv')

print('All done! Outputs saved in', OUT_DIR)