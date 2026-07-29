#!/usr/bin/env python3
"""
Fig7B – Cα RMSF of FcμR chain C with annotated key residues.
Reads paths from config.py, automatically generates residue_map_C.csv.
Handles per‑atom RMSF by averaging per residue.
"""

import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from config import OUT_DIR, DPI, MPL_RCPARAMS, PDB_DIR, DATA_DIR, SYSTEMS, COLORS

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

plt.rcParams.update(MPL_RCPARAMS)

OUTPUT_DIR = OUT_DIR / "fig7b_output"
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

MARKER_RESIDS = [60, 63, 110, 111]
MARKER_LABELS = {60: 'THR60', 63: 'SER63', 110: 'THR110', 111: 'ASP111'}
MARKER_COLOR = '#333333'

FcMuR_CHAIN = 'C'   # FcμR is chain C


def build_residue_map(sys_name):
    """Build residue map for FcμR chain (C). Returns (csv_path, residue_boundaries)."""
    map_file = DATA_DIR / sys_name / 'residue_map_C.csv'
    if map_file.exists():
        map_file.unlink()

    # Find PDB (try PDB_DIR first, then DATA_DIR)
    pdb_path = PDB_DIR / sys_name / 'protein.pdb'
    if not pdb_path.exists():
        pdb_path = DATA_DIR / sys_name / 'protein.pdb'
    if not pdb_path.exists():
        print(f'[WARNING] {sys_name}: PDB not found')
        return None, None

    ca_list = []
    with open(pdb_path) as f:
        for line in f:
            if line.startswith(('ATOM', 'HETATM')) and line[12:16].strip() == 'CA':
                if line[21].strip() == FcMuR_CHAIN:
                    ca_list.append((line[17:20].strip(), int(line[22:26].strip())))

    if not ca_list:
        print(f'[WARNING] {sys_name}: no chain {FcMuR_CHAIN} CA atoms')
        return None, None

    # Build unique residues with boundaries
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

    # Save CSV
    map_file.parent.mkdir(parents=True, exist_ok=True)
    with open(map_file, 'w') as f:
        f.write('index,chain,resname,resid\n')
        for idx, (name, rid, s, e) in enumerate(residues):
            f.write(f'{idx+1},{FcMuR_CHAIN},{name},{rid}\n')
    print(f'[Generated] {sys_name}: {len(residues)} residues -> {map_file}')
    return map_file, residues


def get_chain_start_and_n(pdb_path, chain):
    """Return (start_idx, n_CA) for a given chain in PDB."""
    ca_order = []
    with open(pdb_path) as f:
        for line in f:
            if line.startswith(('ATOM', 'HETATM')) and line[12:16].strip() == 'CA':
                ca_order.append(line[21].strip())
    start = ca_order.index(chain)
    n = ca_order.count(chain)
    return start, n


fig, ax = plt.subplots(figsize=(14, 5))
max_len = 0
final_labels = None
all_data = []
peaks = []

for sys_name in SYSTEMS:
    rmsf_f = DATA_DIR / sys_name / 'rmsf.xvg'
    if not rmsf_f.exists():
        print(f'[WARNING] {sys_name}: rmsf.xvg missing')
        continue

    map_f, res_info = build_residue_map(sys_name)
    if res_info is None:
        print(f'[WARNING] {sys_name}: residue map failed, skipping')
        continue

    pdb_f = PDB_DIR / sys_name / 'protein.pdb'
    if not pdb_f.exists():
        pdb_f = DATA_DIR / sys_name / 'protein.pdb'
    if not pdb_f.exists():
        print(f'[WARNING] {sys_name}: protein.pdb not found')
        continue

    # Load full RMSF
    rmsf_all = np.loadtxt(rmsf_f, comments=['@', '#'])
    start_C, n_C = get_chain_start_and_n(pdb_f, FcMuR_CHAIN)

    # If chain C indices exceed RMSF data, assume RMSF file contains only chain C
    if start_C + n_C > len(rmsf_all):
        print(f'[INFO] {sys_name}: RMSF file appears to contain only chain C (indices exceed). Using all data.')
        y_ca = rmsf_all[:, 1]  # take all rows
    else:
        y_ca = rmsf_all[start_C:start_C + n_C, 1]

    # Check if the number of CA atoms in PDB matches the extracted RMSF data
    n_ca_pdb = sum(1 for _, _, s, e in res_info for _ in range(s, e))  # total CA in residues
    if len(y_ca) != n_ca_pdb:
        print(f'[WARNING] {sys_name}: CA count mismatch (data: {len(y_ca)}, PDB: {n_ca_pdb}), trying to align...')
        # Simple trimming to the shorter one
        min_len = min(len(y_ca), n_ca_pdb)
        y_ca = y_ca[:min_len]
        # Adjust residue boundaries to fit min_len
        trimmed_res = []
        for name, rid, s, e in res_info:
            if s >= min_len:
                break
            e = min(e, min_len)
            if e > s:
                trimmed_res.append((name, rid, s, e))
        res_info = trimmed_res

    # Average RMSF per residue
    y_res = []
    labels = []
    for name, rid, s, e in res_info:
        if s < len(y_ca) and e <= len(y_ca):
            y_res.append(np.mean(y_ca[s:e]))
            labels.append(f"{name}{rid}")

    y_vals = np.array(y_res)
    y_smooth = pd.Series(y_vals).rolling(5, center=True, min_periods=1).mean()
    ax.plot(np.arange(len(y_smooth)), y_smooth,
            color=COLORS[sys_name], lw=1.5, alpha=0.9, label=sys_name)

    for i, val in enumerate(y_vals):
        all_data.append([sys_name, labels[i], val])

    mean_rmsf = np.mean(y_vals)
    std_rmsf = np.std(y_vals)
    threshold = mean_rmsf + 1.5 * std_rmsf
    for idx in np.where(y_vals > threshold)[0]:
        peaks.append([sys_name, labels[idx], y_vals[idx]])

    if len(y_smooth) > max_len:
        max_len = len(y_smooth)
        final_labels = labels

# ---------- Annotations (based on WT residue map) ----------
if final_labels:
    wt_map_f, wt_info = build_residue_map('WT')
    if wt_info:
        wt_resids = [r[1] for r in wt_info]   # list of resid
        y_top = ax.get_ylim()[1] if ax.get_ylim()[1] else 0.5

        for marker_id in MARKER_RESIDS:
            if marker_id in wt_resids:
                idx = wt_resids.index(marker_id)
                ax.axvline(x=idx, color=MARKER_COLOR, linestyle='--',
                           alpha=0.7, linewidth=1.0)
                dx, ha = (-0.5, 'right') if marker_id == 110 else (0.5, 'left')
                ax.text(idx + dx, y_top * 0.90,
                        MARKER_LABELS.get(marker_id, str(marker_id)),
                        rotation=45, ha=ha, va='top', fontsize=7,
                        color=MARKER_COLOR, alpha=0.9,
                        bbox=dict(boxstyle='round,pad=0.1', facecolor='white',
                                  alpha=0.7, edgecolor='none'))

# ---------- Axis formatting ----------
if final_labels:
    step = max(1, max_len // 15)
    ax.set_xticks(np.arange(0, max_len, step))
    ax.set_xticklabels([final_labels[i] for i in np.arange(0, max_len, step)],
                       rotation=45, ha='right', fontsize=8)
ax.set_xlim(0, max_len - 1 if max_len > 0 else 0)
ax.set_xlabel('FcμR-D1 residue')
ax.set_ylabel('Cα RMSF of FcμR-D1 chain (nm)')
if final_labels:
    ax.legend(fontsize=8)
ax.set_ylim(0, None)
ax.grid(False)
for spine in ax.spines.values():
    spine.set_visible(True)
plt.tight_layout()

# ---------- Save ----------
for fmt in ('png', 'jpg', 'pdf'):
    plt.savefig(OUTPUT_DIR / f'Fig7B_RMSF_chainC_annotated.{fmt}',
                dpi=DPI, bbox_inches='tight', facecolor='white')
plt.close()

pd.DataFrame(all_data, columns=['System', 'Residue_Label', 'RMSF (nm)']).to_csv(
    OUTPUT_DIR / 'RMSF_chainC_all_residues.csv', index=False)
pd.DataFrame(peaks, columns=['System', 'Residue_Label', 'RMSF (nm)']).to_csv(
    OUTPUT_DIR / 'RMSF_peak_residues.csv', index=False)

print(f'Figure and CSV files saved to {OUTPUT_DIR}')