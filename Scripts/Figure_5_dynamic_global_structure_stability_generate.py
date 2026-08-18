#!/usr/bin/env python3
"""
Figure 5 – Combined structural dynamics analysis.
Integrates RMSD, Rg, SASA (A–C), RMSF of FcμR-D1 (D), COM distance time series (E),
and equilibrium COM distance bar plot (F). All Outputs saved to fig7_output/.
"""

import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from config import (SYSTEMS, COLORS, DATA_DIR, OUT_DIR, DPI, MPL_RCPARAMS,
                    TIME_MAX_NS, EQUILIBRATION_NS, PDB_DIR)

plt.rcParams.update(MPL_RCPARAMS)

OUTPUT_DIR = OUT_DIR / "fig5_output"
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

SMOOTH_WINDOW = 15
COM_SMOOTH_WINDOW = 10
LINE_WIDTH = 1.5
ALPHA = 0.9
LABEL_FONTSIZE = 9
PANEL_LABEL_POS = (-0.12, 1.06)   # outside top-left, unchanged

# ---------- Helper functions ----------
def get_data_path(sys, fname):
    return DATA_DIR / sys / fname

def load_raw(sys, fname):
    f = get_data_path(sys, fname)
    if not f.exists():
        return None, None
    data = np.loadtxt(f, comments=['@', '#'])
    t = data[:, 0] / 1000.0
    y = data[:, 1]
    mask = t <= TIME_MAX_NS
    t, y = t[mask], y[mask]
    return t, y

def load_and_smooth(sys, fname, window=SMOOTH_WINDOW):
    t, y = load_raw(sys, fname)
    if t is None:
        return None, None
    y_smooth = pd.Series(y).rolling(window, center=True, min_periods=1).mean()
    return t, y_smooth

def get_chain_residues(pdb_path, chain_id):
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

def build_residue_map(sys_name):
    map_file = DATA_DIR / sys_name / 'residue_map_C.csv'
    if map_file.exists():
        map_file.unlink()
    pdb_path = PDB_DIR / sys_name / 'protein.pdb'
    if not pdb_path.exists():
        pdb_path = DATA_DIR / sys_name / 'protein.pdb'
    if not pdb_path.exists():
        return None, None
    residues = get_chain_residues(pdb_path, 'C')
    if not residues:
        return None, None
    map_file.parent.mkdir(parents=True, exist_ok=True)
    with open(map_file, 'w') as f:
        f.write('index,chain,resname,resid\n')
        for idx, (name, rid, s, e) in enumerate(residues):
            f.write(f'{idx+1},C,{name},{rid}\n')
    return map_file, residues

def get_chain_start_and_n(pdb_path, chain):
    ca_order = []
    with open(pdb_path) as f:
        for line in f:
            if line.startswith(('ATOM', 'HETATM')) and line[12:16].strip() == 'CA':
                ca_order.append(line[21].strip())
    start = ca_order.index(chain)
    n = ca_order.count(chain)
    return start, n

# ---------- Data loading ----------
def get_rmsf_data():
    results = []
    max_len = 0
    final_labels = None
    for sys_name in SYSTEMS:
        rmsf_f = DATA_DIR / sys_name / 'rmsf.xvg'
        if not rmsf_f.exists():
            continue
        map_f, res_info = build_residue_map(sys_name)
        if res_info is None:
            continue
        pdb_f = PDB_DIR / sys_name / 'protein.pdb'
        if not pdb_f.exists():
            pdb_f = DATA_DIR / sys_name / 'protein.pdb'
        if not pdb_f.exists():
            continue
        rmsf_all = np.loadtxt(rmsf_f, comments=['@', '#'])
        start_C, n_C = get_chain_start_and_n(pdb_f, 'C')
        if start_C + n_C > len(rmsf_all):
            y_ca = rmsf_all[:, 1]
        else:
            y_ca = rmsf_all[start_C:start_C + n_C, 1]
        n_ca_pdb = sum(1 for _, _, s, e in res_info for _ in range(s, e))
        if len(y_ca) != n_ca_pdb:
            min_len = min(len(y_ca), n_ca_pdb)
            y_ca = y_ca[:min_len]
            trimmed_res = []
            for name, rid, s, e in res_info:
                if s >= min_len:
                    break
                e = min(e, min_len)
                if e > s:
                    trimmed_res.append((name, rid, s, e))
            res_info = trimmed_res
        y_res = []
        labels = []
        for name, rid, s, e in res_info:
            if s < len(y_ca) and e <= len(y_ca):
                y_res.append(np.mean(y_ca[s:e]))
                labels.append(f"{name}{rid}")
        y_vals = np.array(y_res)
        results.append((sys_name, y_vals, labels))
        if len(y_vals) > max_len:
            max_len = len(y_vals)
            final_labels = labels
    return results, max_len, final_labels

def get_com_data():
    ts_smooth = {}
    ts_raw = {}
    means = {}
    cis = {}
    for sys_name in SYSTEMS:
        f = DATA_DIR / sys_name / "com_distance.xvg"
        if not f.exists():
            continue
        d = np.loadtxt(f, comments=["@", "#"])
        t = d[:, 0] / 1000.0
        y = d[:, 1]
        mask = t <= TIME_MAX_NS
        t, y = t[mask], y[mask]
        ts_raw[sys_name] = (t.copy(), y.copy())
        y_smooth = pd.Series(y).rolling(COM_SMOOTH_WINDOW, center=True, min_periods=1).mean()
        ts_smooth[sys_name] = (t, y_smooth)
        eq_mask = t >= EQUILIBRATION_NS
        y_eq = y[eq_mask]
        if len(y_eq) == 0:
            continue
        mean_val = np.mean(y_eq)
        boot = np.random.choice(y_eq, size=(1000, len(y_eq)), replace=True).mean(axis=1)
        ci_low, ci_high = np.percentile(boot, [2.5, 97.5])
        means[sys_name] = mean_val
        cis[sys_name] = (ci_low, ci_high)
    return ts_smooth, ts_raw, means, cis

rmsf_data, rmsf_max_len, rmsf_labels = get_rmsf_data()
com_smooth, com_raw, com_means, com_cis = get_com_data()

# ---------- Data export ----------
# A: RMSD
dfs = []
for sys in SYSTEMS:
    t, y = load_raw(sys, 'rmsd.xvg')
    if t is not None:
        dfs.append(pd.DataFrame({'Time_ns': t, sys: y}))
if dfs:
    df_wide = dfs[0]
    for df in dfs[1:]:
        df_wide = pd.merge(df_wide, df, on='Time_ns', how='outer')
    cols = ['Time_ns'] + [s for s in SYSTEMS if s in df_wide.columns]
    df_wide = df_wide[cols].sort_values('Time_ns').reset_index(drop=True)
    df_wide.to_csv(OUTPUT_DIR / 'Figure_5A_rmsd_time_series.csv', index=False, float_format='%.6f')

# B: Rg
dfs = []
for sys in SYSTEMS:
    t, y = load_raw(sys, 'rg.xvg')
    if t is not None:
        dfs.append(pd.DataFrame({'Time_ns': t, sys: y}))
if dfs:
    df_wide = dfs[0]
    for df in dfs[1:]:
        df_wide = pd.merge(df_wide, df, on='Time_ns', how='outer')
    cols = ['Time_ns'] + [s for s in SYSTEMS if s in df_wide.columns]
    df_wide = df_wide[cols].sort_values('Time_ns').reset_index(drop=True)
    df_wide.to_csv(OUTPUT_DIR / 'Figure_5B_rg_time_series.csv', index=False, float_format='%.6f')

# C: SASA
dfs = []
for sys in SYSTEMS:
    t, y = load_raw(sys, 'sasa.xvg')
    if t is not None:
        dfs.append(pd.DataFrame({'Time_ns': t, sys: y}))
if dfs:
    df_wide = dfs[0]
    for df in dfs[1:]:
        df_wide = pd.merge(df_wide, df, on='Time_ns', how='outer')
    cols = ['Time_ns'] + [s for s in SYSTEMS if s in df_wide.columns]
    df_wide = df_wide[cols].sort_values('Time_ns').reset_index(drop=True)
    df_wide.to_csv(OUTPUT_DIR / 'Figure_5C_sasa_time_series.csv', index=False, float_format='%.6f')

# D: RMSF
rmsf_rows = []
for sys_name, y_vals, labels in rmsf_data:
    for i, lbl in enumerate(labels):
        rmsf_rows.append([sys_name, lbl, y_vals[i]])
pd.DataFrame(rmsf_rows, columns=['System', 'Residue_Label', 'RMSF (nm)']).to_csv(
    OUTPUT_DIR / 'Figure_5D_RMSF_chainC_all_residues.csv', index=False)

# E: COM time series (raw)
ts_dfs = []
for sys_name, (t, y) in com_raw.items():
    ts_dfs.append(pd.DataFrame({'Time_ns': t, sys_name: y}))
if ts_dfs:
    df_ts = ts_dfs[0]
    for df in ts_dfs[1:]:
        df_ts = pd.merge(df_ts, df, on='Time_ns', how='outer')
    cols = ['Time_ns'] + [s for s in SYSTEMS if s in df_ts.columns]
    df_ts = df_ts[cols].sort_values('Time_ns').reset_index(drop=True)
    df_ts.to_csv(OUTPUT_DIR / 'Figure_5E_COM_time_series.csv', index=False, float_format='%.6f')

# F: COM statistics
stat_rows = []
for sys in SYSTEMS:
    if sys in com_means:
        stat_rows.append([sys, com_means[sys], com_cis[sys][0], com_cis[sys][1]])
pd.DataFrame(stat_rows, columns=['System', 'Mean', 'CI_lower', 'CI_upper']).to_csv(
    OUTPUT_DIR / 'Figure_5F_COM_bootstrap_CI.csv', index=False, float_format='%.6f')

# Excel
try:
    with pd.ExcelWriter(OUTPUT_DIR / 'Figure_5_data.xlsx', engine='openpyxl') as writer:
        for prefix, sheet in [('Figure_5A_rmsd_time_series.csv', 'A_RMSD'),
                              ('Figure_5B_rg_time_series.csv', 'B_Rg'),
                              ('Figure_5C_sasa_time_series.csv', 'C_SASA'),
                              ('Figure_5D_RMSF_chainC_all_residues.csv', 'D_RMSF'),
                              ('Figure_5E_COM_time_series.csv', 'E_COM_TS'),
                              ('Figure_5F_COM_bootstrap_CI.csv', 'F_COM_Stats')]:
            f = OUTPUT_DIR / prefix
            if f.exists():
                pd.read_csv(f).to_excel(writer, sheet_name=sheet, index=False)
    print(f"Excel saved to {OUTPUT_DIR / 'Figure_5_all_data.xlsx'}")
except Exception as e:
    print(f"Excel export skipped: {e}")

# ---------- Plotting ----------
fig, axes = plt.subplots(2, 3, figsize=(20, 12))
plt.subplots_adjust(wspace=0.25, hspace=0.25, left=0.06, right=0.98, top=0.94, bottom=0.06)
ax1, ax2, ax3, ax4, ax5, ax6 = axes.flatten()

# Panel labels (outside top-left, unchanged)
def add_panel_label(ax, text):
    ax.text(PANEL_LABEL_POS[0], PANEL_LABEL_POS[1], text,
            transform=ax.transAxes, fontweight='bold', fontsize=LABEL_FONTSIZE,
            va='bottom', ha='center', color='black')

# A: RMSD
add_panel_label(ax1, 'A')
for sys in SYSTEMS:
    t, y = load_and_smooth(sys, 'rmsd.xvg', window=SMOOTH_WINDOW)
    if t is not None:
        ax1.plot(t, y, color=COLORS[sys], lw=LINE_WIDTH, alpha=ALPHA, label=sys)
ax1.axvline(EQUILIBRATION_NS, color='gray', linestyle='--', alpha=0.4)
ax1.set_xlabel('Time (ns)')
ax1.set_ylabel('RMSD (nm)')
ax1.margins(x=0.05)

# B: Rg
add_panel_label(ax2, 'B')
for sys in SYSTEMS:
    t, y = load_and_smooth(sys, 'rg.xvg', window=SMOOTH_WINDOW)
    if t is not None:
        ax2.plot(t, y, color=COLORS[sys], lw=LINE_WIDTH, alpha=ALPHA, label=sys)
ax2.axvline(EQUILIBRATION_NS, color='gray', linestyle='--', alpha=0.4)
ax2.set_xlabel('Time (ns)')
ax2.set_ylabel('Rg (nm)')
ax2.margins(x=0.05)

# C: SASA
add_panel_label(ax3, 'C')
for sys in SYSTEMS:
    t, y = load_and_smooth(sys, 'sasa.xvg', window=SMOOTH_WINDOW)
    if t is not None:
        ax3.plot(t, y, color=COLORS[sys], lw=LINE_WIDTH, alpha=ALPHA, label=sys)
ax3.axvline(EQUILIBRATION_NS, color='gray', linestyle='--', alpha=0.4)
ax3.set_xlabel('Time (ns)')
ax3.set_ylabel('SASA (nm²)')
ax3.margins(x=0.05)

# D: RMSF
add_panel_label(ax4, 'D')
max_len = 0
for sys_name, y_vals, labels in rmsf_data:
    ax4.plot(np.arange(len(y_vals)), y_vals,
             color=COLORS[sys_name], lw=LINE_WIDTH, alpha=ALPHA, label=sys_name)
    if len(y_vals) > max_len:
        max_len = len(y_vals)
        final_labels = labels
if final_labels:
    wt_map_f, wt_info = build_residue_map('WT')
    if wt_info:
        wt_resids = [r[1] for r in wt_info]
        y_top = ax4.get_ylim()[1] if ax4.get_ylim()[1] else 0.5
        marker_resids = [60, 63, 110, 111]
        marker_labels = {60: 'THR60', 63: 'SER63', 110: 'THR110', 111: 'ASP111'}
        for marker_id in marker_resids:
            if marker_id in wt_resids:
                idx = wt_resids.index(marker_id)
                ax4.axvline(x=idx, color='#333333', linestyle='--', alpha=0.7, linewidth=1)
                if marker_id in [60, 110]:
                    dx, ha = (-0.5, 'right')
                else:
                    dx, ha = (0.5, 'left')
                ax4.text(idx + dx, y_top * 0.99, marker_labels[marker_id],
                         rotation=45, ha=ha, va='top', fontsize=LABEL_FONTSIZE,
                         color='#333333', alpha=0.9,
                         bbox=dict(boxstyle='round,pad=0.1', facecolor='white',
                                   alpha=0.7, edgecolor='none'))
ax4.set_xlabel('FcμR-D1 residue')
ax4.set_ylabel('Cα RMSF (nm)')
ax4.margins(x=0.05)

# E: COM distance time series (smoothed)
add_panel_label(ax5, 'E')
for sys_name, (t, y_smooth) in com_smooth.items():
    ax5.plot(t, y_smooth, color=COLORS[sys_name], lw=LINE_WIDTH, alpha=ALPHA, label=sys_name)
ax5.axvline(EQUILIBRATION_NS, color='gray', linestyle='--', alpha=0.4)
ax5.set_xlabel('Time (ns)')
ax5.set_ylabel('COM Distance (nm)')
ax5.margins(x=0.05)

# F: COM distance bar plot (mean ± 95% CI)
add_panel_label(ax6, 'F')
sys_list = [s for s in SYSTEMS if s in com_means]
vals = [com_means[s] for s in sys_list]
err_low = [com_means[s] - com_cis[s][0] for s in sys_list]
err_high = [com_cis[s][1] - com_means[s] for s in sys_list]
ax6.bar(sys_list, vals, yerr=[err_low, err_high],
        color=[COLORS[s] for s in sys_list], alpha=0.85,
        edgecolor='white', capsize=4, width=0.5)
ax6.set_ylabel('Mean COM Distance (nm)')
ax6.set_xlabel('System')
ax6.margins(x=0.05)

# ---------- Add a single legend for all subplots ----------
handles, labels = ax1.get_legend_handles_labels()
if handles:
    fig.legend(handles, labels, loc='upper right', bbox_to_anchor=(0.98, 1.02),
               ncol=1, fontsize=LABEL_FONTSIZE, frameon=False)

# Save figure
basename = 'Figure_5'
for fmt in ('png', 'jpg', 'pdf'):
    out_path = OUTPUT_DIR / f"{basename}.{fmt}"
    fig.savefig(out_path, dpi=DPI, bbox_inches='tight', facecolor='white')
    print(f"Saved: {out_path}")

plt.close()
print(f"All Outputs saved to {OUTPUT_DIR}")