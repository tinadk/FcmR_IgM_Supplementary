#!/usr/bin/env python3
"""
Figure 6 – TICA projections colored by MSM macrostates, system overlay,
and relative free energy (ΔG) per macrostate.
Layout: 2x3 panels:
A–D: Individual systems colored by macrostate (M1–M5) with circles
E: Overlay of all systems colored by system with squares
F: Relative free energy (ΔG) per macrostate
Exports: per-system CSV, combined CSV, and combined Excel with sheets per panel.
"""

import sys
from pathlib import Path
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from matplotlib.lines import Line2D
import pandas as pd
from config import SYSTEMS, COLORS, OUT_DIR, DPI, MPL_RCPARAMS

plt.rcParams.update(MPL_RCPARAMS)

TICA_DIR = ROOT / 'MSM' / 'archive' / 'complete' / 'data' / 'tica'
MSM_DIR = ROOT / 'supplementary_data'
OUTPUT_DIR = OUT_DIR / 'fig6_output'
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

LABEL_FONTSIZE = 9
PANEL_LABEL_POS = (-0.12, 1.06)
STATE_COLORS = ['#E41A1C', '#377EB8', '#4DAF4A', '#984EA3', '#FF7F00']
STATE_LABELS = ['M1', 'M2', 'M3', 'M4', 'M5']

def load_data():
    all_tica, sys_lengths = [], []
    for sys in SYSTEMS:
        tica_path = TICA_DIR / f'{sys}_tica.npy'
        if not tica_path.exists():
            print(f'Warning: {tica_path} not found')
            return None, None
        tica = np.load(tica_path)[:, :2]
        all_tica.append(tica)
        sys_lengths.append(len(tica))

    dtraj_path = MSM_DIR / 'dtraj_k150.npy'
    membership_path = MSM_DIR / 'pcca_membership.npy'
    if not dtraj_path.exists() or not membership_path.exists():
        print(f'Warning: {dtraj_path} or {membership_path} not found')
        return None, None

    dtraj = np.load(dtraj_path)
    membership = np.load(membership_path)

    if membership.ndim == 2:
        microstate_macro = np.argmax(membership, axis=1) + 1
    elif membership.ndim == 1:
        microstate_macro = membership
    else:
        raise ValueError("membership array must be 1D or 2D")

    macro_labels = microstate_macro[dtraj]
    total = sum(sys_lengths)
    if len(macro_labels) != total:
        print(f'Warning: length mismatch: {len(macro_labels)} vs {total}')
        min_len = min(len(macro_labels), total)
        macro_labels = macro_labels[:min_len]

    split = []
    start = 0
    for l in sys_lengths:
        end = start + l
        split.append(macro_labels[start:end])
        start = end
        if start >= len(macro_labels):
            break
    while len(split) < len(SYSTEMS):
        split.append(np.array([], dtype=int))
    return all_tica, split

tica_data, macro_labels = load_data()
if tica_data is None or macro_labels is None:
    raise SystemExit('Data loading failed.')

# Populations for panel F
populations = []
for idx, sys in enumerate(SYSTEMS):
    labels = macro_labels[idx]
    counts = [np.sum(labels == s) for s in range(1, 6)]
    pop = np.array(counts) / len(labels) if len(labels) > 0 else np.zeros(5)
    populations.append(pop)

# ---------- Compute deltaG (with safe handling) ----------
pop_array = np.array(populations)          # shape: (4, 5)
pop_max = np.max(pop_array, axis=1, keepdims=True)   # (4, 1)
pop_max_safe = np.where(pop_max == 0, 1.0, pop_max)
ratio = pop_array / pop_max_safe
ratio[ratio == 0] = 1e-10
kT = 2.494
deltaG = -kT * np.log(ratio)
deltaG = np.nan_to_num(deltaG, nan=0.0, posinf=0.0, neginf=0.0)

# ============================================================
# Export data: per-system CSV, combined CSV, and combined Excel
# ============================================================

# 1) Per-system TICA + macrostate (CSV)
for idx, sys in enumerate(SYSTEMS):
    tica = tica_data[idx]
    labels = macro_labels[idx]
    min_len = min(len(tica), len(labels))
    if len(tica) != min_len:
        tica = tica[:min_len]
        labels = labels[:min_len]
    df = pd.DataFrame({'TIC1': tica[:,0], 'TIC2': tica[:,1], 'Macrostate': labels})
    df.to_csv(OUTPUT_DIR / f'{sys}_tica_macro.csv', index=False, float_format='%.6f')

# 2) Combined data (CSV)
all_rows = []
for idx, sys in enumerate(SYSTEMS):
    tica = tica_data[idx]
    labels = macro_labels[idx]
    min_len = min(len(tica), len(labels))
    for i in range(min_len):
        all_rows.append([sys, tica[i,0], tica[i,1], labels[i]])
df_combined = pd.DataFrame(all_rows, columns=['System', 'TIC1', 'TIC2', 'Macrostate'])
df_combined.to_csv(OUTPUT_DIR / 'Figure_6_TICA_data.csv', index=False, float_format='%.6f')

# 3) DeltaG matrix (CSV)
deltaG_df = pd.DataFrame(deltaG, index=SYSTEMS, columns=STATE_LABELS)
deltaG_df.index.name = 'System'
deltaG_df.to_csv(OUTPUT_DIR / 'Figure_6F_relative_free_energy.csv')

# ============================================================
# Combined Excel with descriptive sheet names (clear and concise)
# ============================================================
excel_path = OUTPUT_DIR / 'Figure_6.xlsx'
with pd.ExcelWriter(excel_path, engine='openpyxl') as writer:
    # 6A–6D: TICA projections for each system
    for idx, sys in enumerate(SYSTEMS):
        tica = tica_data[idx]
        labels = macro_labels[idx]
        min_len = min(len(tica), len(labels))
        df = pd.DataFrame({
            'TIC1': tica[:min_len, 0],
            'TIC2': tica[:min_len, 1],
            'Macrostate': labels[:min_len]
        })
        # Clear sheet name: panel letter + system + "_TICA"
        sheet_name = f"{chr(65+idx)}_{sys}_TICA"
        df.to_excel(writer, sheet_name=sheet_name, index=False)

    # 6E: overlay of all systems (system-colored squares)
    df_e = pd.DataFrame(all_rows, columns=['System', 'TIC1', 'TIC2', 'Macrostate'])
    df_e.to_excel(writer, sheet_name='E_Overlay_AllSystems', index=False)

    # 6F: relative free energy per macrostate (ΔG)
    df_f = pd.DataFrame(deltaG, index=SYSTEMS, columns=STATE_LABELS)
    df_f.index.name = 'System'
    df_f.to_excel(writer, sheet_name='F_DeltaG', index=False)

print(f"Combined Excel with sheets (A_WT_TICA, B_D111A_TICA, C_T110A_D111A_TICA, D_T60A_S63A_TICA, E_Overlay_AllSystems, F_DeltaG) saved to: {excel_path}")

# ---------- Create figure and axes ----------
fig, axes = plt.subplots(2, 3, figsize=(20, 12))
plt.subplots_adjust(wspace=0.25, hspace=0.25, left=0.06, right=0.98, top=0.94, bottom=0.06)

def add_label(ax, text):
    ax.text(PANEL_LABEL_POS[0], PANEL_LABEL_POS[1], text,
            transform=ax.transAxes, fontweight='bold', fontsize=LABEL_FONTSIZE,
            va='bottom', ha='center', color='black')

# A–D: Individual systems by macrostate (circles)
for idx, sys in enumerate(SYSTEMS):
    ax = axes.flatten()[idx]
    tica = tica_data[idx]
    labels = macro_labels[idx]
    min_len = min(len(tica), len(labels))
    tica = tica[:min_len]
    labels = labels[:min_len]
    for s in range(1, 6):
        mask = (labels == s)
        if np.sum(mask) > 0:
            ax.scatter(tica[mask,0], tica[mask,1],
                       c=STATE_COLORS[s-1], marker='o', s=2, alpha=0.6,
                       label=STATE_LABELS[s-1] if idx == 0 else "")
    ax.set_title(sys, fontweight='bold')
    ax.set_xlabel('TIC 1')
    ax.set_ylabel('TIC 2')
    ax.margins(x=0.05, y=0.05)
    add_label(ax, chr(65+idx))

# Add macrostate legend only in panel A (lower right)
ax1 = axes.flatten()[0]
from matplotlib.lines import Line2D
macro_handles = [Line2D([0], [0], marker='o', color='w', markerfacecolor=c, markersize=6, label=l)
                 for c, l in zip(STATE_COLORS, STATE_LABELS)]
leg_a = ax1.legend(handles=macro_handles, loc='lower right',
                   fontsize=LABEL_FONTSIZE, frameon=False, title='Macrostate')
leg_a._legend_box.align = "left"

# E: Overlay by system (squares)
ax5 = axes.flatten()[4]
for sys in SYSTEMS:
    idx = SYSTEMS.index(sys)
    tica = tica_data[idx]
    ax5.scatter(tica[:,0], tica[:,1], c=COLORS[sys], marker='s', s=2, alpha=0.4, label=sys, rasterized=True)
ax5.set_xlabel('TIC 1')
ax5.set_ylabel('TIC 2')
ax5.margins(x=0.05, y=0.05)
leg_e = ax5.legend(loc='lower right', fontsize=LABEL_FONTSIZE, frameon=False, title='System')
leg_e._legend_box.align = "left"
add_label(ax5, 'E')

# F: Relative free energy (ΔG) per macrostate
ax6 = axes.flatten()[5]
x = np.arange(len(STATE_LABELS))
bar_width = 0.15
all_bars_f = []
for si, sys_name in enumerate(SYSTEMS):
    bars = ax6.bar(x + si * bar_width, deltaG[si, :], bar_width,
                   color=COLORS[sys_name], label=sys_name)
    all_bars_f.append(bars)

ax6.set_xticks(x + bar_width * 1.5)
ax6.set_xticklabels(STATE_LABELS)
ax6.set_xlabel('Macrostate')
ax6.set_ylabel('ΔG (kJ/mol)')
ax6.margins(x=0.05)
leg_f = ax6.legend(loc='upper right', fontsize=LABEL_FONTSIZE, title='System', frameon=False)
leg_f._legend_box.align = "left"
add_label(ax6, 'F')

# Define offsets per macrostate (system order: WT, D111A, T110A_D111A, T60A_S63A)
offset_dict = {
    'M1': [-0.05, 0.05, 0.03, 0.03],
    'M2': [-0.05, 0.05, 0.03, 0.03],
    'M3': [0.0, -0.03, -0.03, 0.0],
    'M4': [0.0, 0.0, 0.0, 0.0],
    'M5': [-0.14, -0.05, 0.05, 0.14]
}

for state_idx, state_label in enumerate(STATE_LABELS):
    offsets = offset_dict.get(state_label, [0.0, 0.0, 0.0, 0.0])
    for si in range(len(SYSTEMS)):
        bar = all_bars_f[si][state_idx]
        height = deltaG[si, state_idx]
        if height < 0.01:
            continue
        x_center = bar.get_x() + bar.get_width() / 2.0
        ax6.text(x_center + offsets[si], height + 0.2,
                f'{height:.1f}', ha='center', va='bottom', fontsize=7)

ymax = max(75, np.nanmax(deltaG) * 1.35)
if np.isfinite(ymax) and ymax > 0:
    ax6.set_ylim(0, ymax)
else:
    ax6.set_ylim(0, 75)

for ext in ('png', 'jpg', 'pdf'):
    fig.savefig(OUTPUT_DIR / f'Figure_6.{ext}', dpi=DPI, bbox_inches='tight', facecolor='white')
plt.close()
print(f'Figure 6 plots saved to {OUTPUT_DIR}')