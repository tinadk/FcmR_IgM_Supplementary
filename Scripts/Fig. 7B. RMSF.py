#!/usr/bin/env python3
"""
RMSF of FcμR chain (B) with annotated key residues.
- 60, 63, 111 labels placed right of the line.
- 110 label placed left of the line.
Output: 1200 dpi PNG/JPG/PDF + CSV files.
"""

import os, numpy as np, matplotlib.pyplot as plt, pandas as pd

# ================= CONFIGURATION =================
SYSTEMS = ['WT', 'D111A', 'T110A_D111A', 'T60A_S63A']
COLORS = {
    'WT': '#2C3E50',
    'D111A': '#E74C3C',
    'T110A_D111A': '#27AE60',
    'T60A_S63A': '#2980B9'
}

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
BASE_DIR = os.path.dirname(SCRIPT_DIR)
DATA_DIR = os.path.join(BASE_DIR, 'data')
OUT_DIR  = os.path.join(BASE_DIR, 'output')
PDB_BASE = os.path.join(BASE_DIR, 'data')

DPI = 1600
os.makedirs(OUT_DIR, exist_ok=True)

# ---------- Key residues to annotate (FcμR numbering) ----------
MARKER_RESIDS = [60, 63, 110, 111]                     # original residue IDs
MARKER_LABELS = {60: 'THR60', 63: 'SER63', 110: 'THR110', 111: 'ASP111'}
MARKER_COLOR = '#333333'

# ================= HELPER: find chain B indices =================

# ================= MAIN PLOT =================
fig, ax = plt.subplots(figsize=(14, 5))
max_len = 0
final_labels = None
all_data = []      # for RMSF_chainB_all_residues.csv
peaks = []         # for RMSF_peak_residues.csv

for sys in SYSTEMS:
    rmsf_f = os.path.join(DATA_DIR, sys, 'rmsf.xvg')
    map_f  = os.path.join(DATA_DIR, sys, 'residue_map_B.csv')
    pdb_f = os.path.join(PDB_BASE, sys, 'protein.pdb')

    if not all(map(os.path.exists, [rmsf_f, map_f, pdb_f])):
        print(f'[WARNING] {sys}: missing files, skipping')
        continue

    # 1) load full RMSF data (all chains)
    rmsf_all = np.loadtxt(rmsf_f, comments=['@','#'])
    # 2) get chain B start index
    start_B, n_B = get_chain_B_start_and_n(pdb_f)
    y_vals_B = rmsf_all[start_B:start_B + n_B, 1]

    # 3) load residue labels for chain B (CSV)
    res_df = pd.read_csv(map_f)
    if n_B != len(res_df):
        n = min(n_B, len(res_df))
        y_vals_B = y_vals_B[:n]
        res_df = res_df.iloc[:n]
    labels = [f"{row.resname}{row.resid}" for _, row in res_df.iterrows()]

    # 4) smooth and plot
    y_smooth = pd.Series(y_vals_B).rolling(5, center=True, min_periods=1).mean()
    ax.plot(np.arange(len(y_smooth)), y_smooth, color=COLORS[sys],
            lw=1.5, alpha=0.9, label=sys)

    # 5) store data for CSV
    for i, (_, row) in enumerate(res_df.iterrows()):
        all_data.append([sys, labels[i], y_vals_B[i]])

    # 6) identify peak residues (RMSF > mean + 1.5*std)
    mean_rmsf = np.mean(y_vals_B)
    std_rmsf = np.std(y_vals_B)
    threshold = mean_rmsf + 1.5 * std_rmsf
    peak_mask = y_vals_B > threshold
    for idx in np.where(peak_mask)[0]:
        peaks.append([sys, labels[idx], y_vals_B[idx]])

    # keep labels from the longest system (for x-axis)
    if len(y_smooth) > max_len:
        max_len = len(y_smooth)
        final_labels = labels

# ================= ANNOTATIONS =================
if final_labels:
    wt_map_file = os.path.join(DATA_DIR, 'WT', 'residue_map_B.csv')
    if os.path.exists(wt_map_file):
        wt_res_df = pd.read_csv(wt_map_file)
        wt_resids = wt_res_df['resid'].values
        y_top = ax.get_ylim()[1] if ax.get_ylim()[1] else 0.5

        for marker_id in MARKER_RESIDS:
            matches = np.where(wt_resids == marker_id)[0]
            if len(matches) == 0:
                continue
            idx = matches[0]   # index in chain B (0‑based)

            # draw vertical line on all systems
            ax.axvline(x=idx, color=MARKER_COLOR, linestyle='--',
                       alpha=0.7, linewidth=1.0)

            # Decide text placement: residue 110 goes left, others go right
            if marker_id == 110:
                dx = -0.5
                ha = 'right'
            else:
                dx = +0.5
                ha = 'left'

            label_text = MARKER_LABELS.get(marker_id, str(marker_id))
            ax.text(idx + dx, y_top * 0.90, label_text,
                    rotation=45, ha=ha, va='top',
                    fontsize=7, color=MARKER_COLOR, alpha=0.9,
                    bbox=dict(boxstyle='round,pad=0.1', facecolor='white',
                              alpha=0.7, edgecolor='none'))

# ================= AXIS FORMATTING =================
if final_labels:
    step = max(1, max_len // 15)
    ax.set_xticks(np.arange(0, max_len, step))
    ax.set_xticklabels([final_labels[i] for i in np.arange(0, max_len, step)],
                       rotation=45, ha='right', fontsize=8)
ax.set_xlim(0, max_len - 1)
ax.set_xlabel('FcμR-D1 residue')
ax.set_ylabel('Cα RMSF of FcμR-D1 chain (nm)')
# ax.set_title('Cα RMSF of FcμR chain (5‑residue smoothing)')
ax.legend(fontsize=8)
ax.set_ylim(0, None)
ax.grid(False)
for spine in ax.spines.values():
    spine.set_visible(True)

plt.tight_layout()

# ================= SAVE FIGURES =================
for fmt in ('png', 'jpg', 'pdf'):
    plt.savefig(os.path.join(OUT_DIR, f'S02_RMSF_chainB_annotated.{fmt}'),
                dpi=DPI, bbox_inches='tight', facecolor='white')
plt.close()

# ================= SAVE CSV TABLES =================
pd.DataFrame(all_data, columns=['System', 'Residue_Label', 'RMSF (nm)']).to_csv(
    os.path.join(OUT_DIR, 'RMSF_chainB_all_residues.csv'), index=False)
pd.DataFrame(peaks, columns=['System', 'Residue_Label', 'RMSF (nm)']).to_csv(
    os.path.join(OUT_DIR, 'RMSF_peak_residues.csv'), index=False)

print('Annotated figure and CSV tables saved. Check rmsf_true_residue/output/')