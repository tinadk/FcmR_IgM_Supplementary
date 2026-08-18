#!/usr/bin/env python3
"""
Figure 1A – Cα RMSD heatmap among FcμR structures.
Reads PDB files from config.PDB_DIR, saves figure (png, jpg, pdf) and data (CSV, XLSX)
to config output folder. No top title. Colorbar height matches heatmap exactly.
"""

import sys
from pathlib import Path
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from config import PDB_DIR, OUT_DIR, DPI, MPL_RCPARAMS

import numpy as np
import MDAnalysis as mda
from MDAnalysis.analysis import rms
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd

plt.rcParams.update(MPL_RCPARAMS)

# ---------- Parameters ----------
structures = [
    "7YTE.pdb",
    "7YSG.pdb",
    "7YTC.pdb",
    "7YTD.pdb",
    "8BPE.pdb",
    "8BPF.pdb",
    "8BPG.pdb",
]

chain_ids = ["C", "R", "R", "R", "I", "I", "A"]
residue_range = "resid 18-124 and name CA"

OUTPUT_DIR = OUT_DIR / 'fig1a_output'
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

# ---------- Calculate the RMSD matrix ----------
resid_sets = []
atoms_list = []

for pdb_file, chain_id in zip(structures, chain_ids):
    pdb_path = PDB_DIR / pdb_file
    u = mda.Universe(str(pdb_path))
    sel = f"chainID {chain_id} and {residue_range}"
    atoms = u.select_atoms(sel)
    if len(atoms) == 0:
        raise ValueError(f"Chain {chain_id} residues 18-124 not found in {pdb_path}")
    resid_sets.append(set(atoms.resids))
    atoms_list.append(atoms)
    print(f"{chain_id} {pdb_file}: {atoms.positions.shape}, resids: {min(atoms.resids)}-{max(atoms.resids)}")

common_resids = sorted(set.intersection(*resid_sets))
print(f"\nCommon residues count: {len(common_resids)}")

coords_list = []
for atoms in atoms_list:
    resid_to_coord = {resid: pos for resid, pos in zip(atoms.resids, atoms.positions)}
    filtered_coords = np.array([resid_to_coord[resid] for resid in common_resids])
    coords_list.append(filtered_coords)
    print(filtered_coords.shape)

n = len(structures)
rmsd_matrix = np.zeros((n, n))
for i in range(n):
    for j in range(i, n):
        if i == j:
            rmsd_matrix[i, j] = 0.0
        else:
            rmsd_val = rms.rmsd(coords_list[i], coords_list[j], superposition=True)
            rmsd_matrix[i, j] = rmsd_val
            rmsd_matrix[j, i] = rmsd_val

print("\nRMSD matrix (Å):")
print(np.round(rmsd_matrix, 2))

# ---------- Draw heatmap ----------
labels = ["7YTE", "7YSG", "7YTC", "7YTD", "8BPE", "8BPF", "8BPG"]

fig, ax = plt.subplots(figsize=(6, 5))
sns.heatmap(
    rmsd_matrix,
    annot=True,
    fmt='.2f',
    cmap='viridis',
    xticklabels=labels,
    yticklabels=labels,
    square=True,
    linewidths=0.5,
    cbar_kws={'label': 'Cα RMSD (Å)', 'shrink': 1.0},
    annot_kws={'size': 9},
    ax=ax
)
ax.set_xticklabels(labels, rotation=45, ha='right', fontsize=8)
ax.set_yticklabels(labels, rotation=0, fontsize=8)
ax.set_xlabel('PDB ID', fontsize=8)
ax.set_ylabel('PDB ID', fontsize=7)

plt.tight_layout()

cbar = ax.collections[0].colorbar
cax = cbar.ax
pos_ax = ax.get_position()
pos_cax = cax.get_position()
cax.set_position([pos_cax.x0, pos_ax.y0, pos_cax.width, pos_ax.height])

for ext in ('png', 'jpg', 'pdf'):
    plt.savefig(OUTPUT_DIR / f'Figure_1A_rmsd_heatmap.{ext}', dpi=DPI, bbox_inches='tight')
plt.close()
print(f"Heatmap saved to {OUTPUT_DIR}")

# ---------- Export data table ----------
df_rmsd = pd.DataFrame(rmsd_matrix, index=labels, columns=labels)
# Set index name for Excel export
df_rmsd.index.name = 'FcμR-D1 and IgM Structure PDB'

csv_path = OUTPUT_DIR / 'Figure_1A_RMSD_matrix.csv'
df_rmsd.to_csv(csv_path, float_format='%.2f')
print(f"CSV saved to {csv_path}")

xlsx_path = OUTPUT_DIR / 'Figure_1A_RMSD_matrix_data.xlsx'
df_rmsd.to_excel(xlsx_path, float_format='%.2f')
print(f"XLSX saved to {xlsx_path}")