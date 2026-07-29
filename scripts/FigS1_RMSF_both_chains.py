#!/usr/bin/env python3
"""
Fig S1 – RMSF of both chains: IgM-Fcμ-Cμ4 (chain B) and FcμR-D1 (chain C).
Uses the same per-residue averaging logic as Fig7B, concatenates chains left to right.
Reads paths from config.py.  No top title.  Exports CSV and XLSX.
"""

import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from config import OUT_DIR, DPI, MPL_RCPARAMS, PDB_DIR, DATA_DIR, SYSTEMS, COLORS

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

plt.rcParams.update(MPL_RCPARAMS)

OUTPUT_DIR = OUT_DIR / "FigS1_RMSF_both_chains_output"
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

fig, ax = plt.subplots(figsize=(16, 5))
max_len = 0
all_chainB_labels = None
all_chainC_labels = None
final_labels = None
all_data = []

for sys_name in SYSTEMS:
    rmsf_f = DATA_DIR / sys_name / 'rmsf.xvg'
    pdb_f = PDB_DIR / sys_name / 'protein.pdb'
    if not pdb_f.exists():
        pdb_f = DATA_DIR / sys_name / 'protein.pdb'
    if not rmsf_f.exists() or not pdb_f.exists():
        print(f'[WARNING] {sys_name}: missing files, skipping')
        continue

    resB = extract_chain_rmsf(rmsf_f, pdb_f, 'B')
    if resB is None:
        print(f'[WARNING] {sys_name}: chain B not found')
        continue
    labelsB, yB, _ = resB

    resC = extract_chain_rmsf(rmsf_f, pdb_f, 'C')
    if resC is None:
        print(f'[WARNING] {sys_name}: chain C not found')
        continue
    labelsC, yC, _ = resC

    y_vals = np.concatenate([yB, yC])
    labels = labelsB + labelsC

    y_smooth = pd.Series(y_vals).rolling(5, center=True, min_periods=1).mean()
    ax.plot(np.arange(len(y_smooth)), y_smooth,
            color=COLORS[sys_name], lw=1.5, alpha=0.9, label=sys_name)

    for i, lbl in enumerate(labels):
        all_data.append([sys_name, lbl, y_vals[i]])

    if len(y_smooth) > max_len:
        max_len = len(y_smooth)
        final_labels = labels
        all_chainB_labels = labelsB
        all_chainC_labels = labelsC
    elif len(y_smooth) == max_len:
        final_labels = labels
        all_chainB_labels = labelsB
        all_chainC_labels = labelsC

# ---------- Chain boundary ----------
if all_chainB_labels is not None and all_chainC_labels is not None:
    nB = len(all_chainB_labels)
    boundary = nB - 0.5
    ax.axvline(x=boundary, color='gray', linestyle='--', alpha=0.7, linewidth=1)

    mid_B = (0 + nB - 1) / 2
    mid_C = nB + (0 + len(all_chainC_labels) - 1) / 2
    ax.text(mid_B, -0.32, 'IgM-Fcμ-Cμ4', transform=ax.get_xaxis_transform(),
            ha='center', va='top', fontsize=9, color='gray')
    ax.text(mid_C, -0.32, 'FcμR-D1', transform=ax.get_xaxis_transform(),
            ha='center', va='top', fontsize=9, color='gray')

ax.set_xlabel('Residue', labelpad=40)

# ---------- Axis formatting ----------
if final_labels:
    step = max(1, max_len // 20)
    ax.set_xticks(np.arange(0, max_len, step))
    ax.set_xticklabels([final_labels[i] for i in np.arange(0, max_len, step)],
                       rotation=45, ha='right', fontsize=7)
ax.set_xlim(0, max_len - 1)
ax.set_ylabel('RMSF (nm)')
# No title – removed
ax.legend(fontsize=8, loc='upper right', framealpha=0.6, edgecolor='none')
ax.set_ylim(0, None)
ax.grid(False)
for spine in ax.spines.values():
    spine.set_visible(True)

plt.subplots_adjust(bottom=0.32)
plt.tight_layout(rect=[0, 0.22, 1, 1])

# ---------- Save figure ----------
for fmt in ('png', 'jpg', 'pdf'):
    plt.savefig(OUTPUT_DIR / f'FigS1_RMSF_both_chains.{fmt}',
                dpi=DPI, bbox_inches='tight', facecolor='white')
plt.close()

# ---------- Export data ----------
df = pd.DataFrame(all_data, columns=['System', 'Residue_Label', 'RMSF (nm)'])
csv_path = OUTPUT_DIR / 'RMSF_both_chains_all_residues.csv'
df.to_csv(csv_path, index=False)
print(f'CSV saved to {csv_path}')

xlsx_path = OUTPUT_DIR / 'RMSF_both_chains_all_residues.xlsx'
df.to_excel(xlsx_path, index=False)
print(f'XLSX saved to {xlsx_path}')

print(f'FigS1 and data files saved to {OUTPUT_DIR}')