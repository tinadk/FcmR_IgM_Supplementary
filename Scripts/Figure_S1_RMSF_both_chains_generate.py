#!/usr/bin/env python3
"""
Figure_S1 – RMSF of both chains: IgM-Fcμ-Cμ4 (chain B) and FcμR-D1 (chain C).
Uses the same per-residue averaging logic as Fig7B, concatenates chains left to right.
Reads paths from config.py.  No top title.  Exports CSV and XLSX.
"""

import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from config import OUT_DIR, DPI, MPL_RCPARAMS, PDB_DIR, DATA_DIR, SYSTEMS, COLORS

plt.rcParams.update(MPL_RCPARAMS)

# ---------- Helper functions (inline) ----------
def get_chain_residues(pdb_path, chain_id):
    """Return list of (resname, resid, start_idx, end_idx) for a chain."""
    ca_list = []
    with open(pdb_path) as f:
        for line in f:
            if line.startswith(('ATOM', 'HETATM')) and line[12:16].strip() == 'CA':
                if line[21].strip() == chain_id:
                    ca_list.append((line[17:20].strip(), int(line[22:26].strip())))
    if not ca_list:
        return []
    residues = []
    cur_res = None
    start = 0
    for i, (name, rid) in enumerate(ca_list):
        if (name, rid) != cur_res:
            if cur_res is not None:
                residues.append((cur_res[0], cur_res[1], start, i))
            cur_res = (name, rid)
            start = i
    residues.append((cur_res[0], cur_res[1], start, len(ca_list)))
    return residues

def extract_chain_rmsf(rmsf_file, pdb_path, chain_id):
    residues = get_chain_residues(pdb_path, chain_id)
    if not residues:
        return None
    rmsf_all = np.loadtxt(rmsf_file, comments=['@', '#'])
    ca_order = []
    with open(pdb_path) as f:
        for line in f:
            if line.startswith(('ATOM', 'HETATM')) and line[12:16].strip() == 'CA':
                ca_order.append(line[21].strip())
    start = ca_order.index(chain_id)
    n_ca = ca_order.count(chain_id)
    if start + n_ca > len(rmsf_all):
        y_ca = rmsf_all[:, 1]
    else:
        y_ca = rmsf_all[start:start + n_ca, 1]
    y_res = []
    labels = []
    for name, rid, s, e in residues:
        if s < len(y_ca) and e <= len(y_ca):
            y_res.append(np.mean(y_ca[s:e]))
            labels.append(f"{name}{rid}")
    return labels, np.array(y_res), y_ca

# ---------- Annotation parameters (same as Fig7B) ----------
MARKER_RESIDS = [60, 63, 110, 111]
MARKER_LABELS = {60: 'THR60', 63: 'SER63', 110: 'THR110', 111: 'ASP111'}
MARKER_COLOR = '#333333'

# ---------- Main plotting ----------
OUTPUT_DIR = OUT_DIR / "Figure_S1_RMSF_both_chains_output"
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

fig, ax = plt.subplots(figsize=(16, 5), constrained_layout=True)
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

    y_vals = np.concatenate([yB, [np.nan], yC])
    labels = labelsB + [''] + labelsC

    # Directly use raw data (no smoothing)
    ax.plot(np.arange(len(y_vals)), y_vals,
            color=COLORS[sys_name], lw=1.5, alpha=0.9, label=sys_name)
    for i, lbl in enumerate(labels):
        all_data.append([sys_name, lbl, y_vals[i]])

    if len(y_vals) > max_len:
        max_len = len(y_vals)
        final_labels = labels
        all_chainB_labels = labelsB
        all_chainC_labels = labelsC
    elif len(y_vals) == max_len:
        final_labels = labels
        all_chainB_labels = labelsB
        all_chainC_labels = labelsC

# ---------- Chain boundary ----------
if all_chainB_labels is not None and all_chainC_labels is not None:
    nB = len(all_chainB_labels)
    boundary = nB
    ax.axvline(x=boundary, color='#7F8C8D', linestyle='-.', alpha=0.8, linewidth=1.5)

    mid_B = (0 + nB - 1) / 2
    mid_C = nB + (0 + len(all_chainC_labels) - 1) / 2
    ax.text(mid_B, -0.14, 'IgM-Fcμ-Cμ4', transform=ax.get_xaxis_transform(),
            ha='center', va='top', fontsize=9, color='gray')
    ax.text(mid_C, -0.14, 'FcμR-D1', transform=ax.get_xaxis_transform(),
            ha='center', va='top', fontsize=9, color='gray')

# ---------- Annotations for key residues on FcμR-D1 (C chain) ----------
if all_chainB_labels is not None and all_chainC_labels is not None:
    nB = len(all_chainB_labels)
    # Get WT PDB to resolve residue positions (use WT as reference)
    wt_pdb = PDB_DIR / 'WT' / 'protein.pdb'
    if not wt_pdb.exists():
        wt_pdb = DATA_DIR / 'WT' / 'protein.pdb'
    if wt_pdb.exists():
        wt_residues = get_chain_residues(wt_pdb, 'C')
        if wt_residues:
            # Build a dict from residue label string to index in chain C
            # We'll just search in all_chainC_labels
            for marker_id in MARKER_RESIDS:
                label_str = MARKER_LABELS[marker_id]  # e.g., 'THR60'
                if label_str in all_chainC_labels:
                    c_idx = all_chainC_labels.index(label_str)
                    abs_idx = nB + c_idx
                    ax.axvline(x=abs_idx, color=MARKER_COLOR, linestyle='--',
                               alpha=0.7, linewidth=1.0)
                    # Position text: left for 60 and 110, right for 63 and 111
                    if marker_id in [60, 110]:
                        dx, ha = (-0.5, 'right')
                    else:
                        dx, ha = (0.5, 'left')
                    y_top = ax.get_ylim()[1] if ax.get_ylim()[1] else 0.5
                    ax.text(abs_idx + dx, y_top * 0.99,
                            MARKER_LABELS.get(marker_id, str(marker_id)),
                            rotation=45, ha=ha, va='top', fontsize=9,
                            color=MARKER_COLOR, alpha=0.9,
                            bbox=dict(boxstyle='round,pad=0.1', facecolor='white',
                                      alpha=0.7, edgecolor='none'))

ax.set_xlabel('Residue', labelpad=15, fontsize=9)

# ---------- Axis formatting ----------
if final_labels:
    step = max(1, max_len // 20)
    ax.set_xticks(np.arange(0, max_len, step))
    ax.set_xticklabels([final_labels[i] for i in np.arange(0, max_len, step)],
                       rotation=45, ha='right', fontsize=9)
ax.margins(x=0.03)
ax.set_ylabel('RMSF (nm)')
ax.legend(fontsize=9, loc='upper left', framealpha=0.6, edgecolor='none')
ax.set_ylim(0, None)
ax.grid(False)
for spine in ax.spines.values():
    spine.set_visible(True)

# ---------- Save figure ----------
for fmt in ('png', 'jpg', 'pdf'):
    plt.savefig(OUTPUT_DIR / f'Figure_S1_RMSF_both_chains.{fmt}',
                dpi=DPI, bbox_inches='tight', facecolor='white')
plt.close()

# ---------- Export data ----------
# First column remains "System" (no renaming)
df = pd.DataFrame(all_data, columns=['System', 'Residue_Label', 'RMSF (nm)'])

csv_path = OUTPUT_DIR / 'Figure_S1_RMSF_both_chains_all_residues.csv'
df.to_csv(csv_path, index=False)
print(f'CSV saved to {csv_path}')

xlsx_path = OUTPUT_DIR / 'Figure_S1_RMSF_both_chains_all_residues_data.xlsx'
df.to_excel(xlsx_path, index=False)
print(f'XLSX saved to {xlsx_path}')

print(f'Figure_S1 and data files saved to {OUTPUT_DIR}')