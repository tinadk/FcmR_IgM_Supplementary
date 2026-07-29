#!/usr/bin/env python3
"""
Fig. 7C: COM distance between FcμR and IgM.
Combined figure: time series (left) and equilibrium bar plot (right).
All parameters from config.py.
"""

import sys
from pathlib import Path
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from config import (SYSTEMS, DATA_DIR, OUT_DIR, TIME_MAX_NS, COLORS, DPI,
                    MPL_RCPARAMS, EQUILIBRATION_NS)

# ------------------- Output directory -------------------
OUTPUT_DIR = OUT_DIR / "fig7c_output"
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

plt.rcParams.update(MPL_RCPARAMS)

# ---------- Data collection ----------
means = {}
cis = {}
time_series = {}   # sys_name -> (t, y_smooth)
raw_ts = {}        # sys_name -> (t_raw, y_raw) for export

for sys_name in SYSTEMS:
    f = DATA_DIR / sys_name / "com_distance.xvg"
    if not f.exists():
        print(f"Warning: {f} not found, skipping {sys_name}")
        continue

    d = np.loadtxt(f, comments=["@", "#"])
    t = d[:, 0] / 1000.0          # ps → ns
    y = d[:, 1]

    mask = t <= TIME_MAX_NS
    t, y = t[mask], y[mask]

    # Store raw data for CSV export
    raw_ts[sys_name] = (t.copy(), y.copy())

    eq_mask = t >= EQUILIBRATION_NS
    y_eq = y[eq_mask]

    # Smoothed curve for plotting
    y_smooth = pd.Series(y).rolling(50, center=True, min_periods=1).mean()
    time_series[sys_name] = (t, y_smooth)

    mean_val = np.mean(y_eq)
    boot = np.random.choice(y_eq, size=(1000, len(y_eq)), replace=True).mean(axis=1)
    ci_low, ci_high = np.percentile(boot, [2.5, 97.5])

    means[sys_name] = mean_val
    cis[sys_name] = (ci_low, ci_high)

# ---------- Combined figure ----------
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 5))

# Left panel: time series
for sys_name, (t, y_smooth) in time_series.items():
    ax1.plot(t, y_smooth, color=COLORS[sys_name], linewidth=1.5, alpha=0.85, label=sys_name)

ax1.axvline(EQUILIBRATION_NS, color='gray', linestyle='--', alpha=0.4)
ax1.set_xlabel("Time (ns)")
ax1.set_ylabel("COM Distance (nm)")
ax1.set_title("COM distance between Fc\u03bcR and IgM")
ax1.legend(fontsize=6.6, loc='upper right', framealpha=0.6, edgecolor='none')
ax1.grid(False)
for spine in ax1.spines.values():
    spine.set_visible(True)

# Right panel: bar plot with 95% CI
labels = [s for s in SYSTEMS if s in means]
vals = [means[s] for s in labels]
err_low = [means[s] - cis[s][0] for s in labels]
err_high = [cis[s][1] - means[s] for s in labels]

ax2.bar(
    labels,
    vals,
    yerr=[err_low, err_high],
    color=[COLORS[s] for s in labels],
    alpha=0.85,
    edgecolor='white',
    capsize=4
)

ax2.set_ylabel("Mean COM Distance (nm)")
ax2.set_title("Equilibrium COM distance (mean \u00b1 95% CI)")
ax2.grid(False)
for spine in ax2.spines.values():
    spine.set_visible(True)

plt.tight_layout()

# Save figures
for fmt in ('png', 'jpg', 'pdf'):
    plt.savefig(OUTPUT_DIR / f"Fig7C_COM_distance_combined.{fmt}",
                dpi=DPI, facecolor='white', edgecolor='none')
plt.close()

# ---------- Export data ----------
# Build wide-format time series DataFrame
ts_dfs = []
for sys_name in SYSTEMS:
    if sys_name in raw_ts:
        t, y = raw_ts[sys_name]
        df_sys = pd.DataFrame({'Time_ns': t, sys_name: y})
        ts_dfs.append(df_sys)

if ts_dfs:
    df_ts_wide = ts_dfs[0]
    for df_sys in ts_dfs[1:]:
        df_ts_wide = pd.merge(df_ts_wide, df_sys, on='Time_ns', how='outer')
    cols = ['Time_ns'] + [s for s in SYSTEMS if s in df_ts_wide.columns]
    df_ts_wide = df_ts_wide[cols].sort_values('Time_ns').reset_index(drop=True)

    ts_csv = OUTPUT_DIR / "Fig7C_COM_time_series.csv"
    df_ts_wide.to_csv(ts_csv, index=False, float_format='%.6f')
    print(f"Time series CSV saved to {ts_csv}")

# Build statistics DataFrame
if means:
    stat_rows = []
    for sys in SYSTEMS:
        if sys in means:
            stat_rows.append({
                'System': sys,
                'Mean': means[sys],
                'CI_lower': cis[sys][0],
                'CI_upper': cis[sys][1]
            })
    df_stats = pd.DataFrame(stat_rows)

    stat_csv = OUTPUT_DIR / "Fig7C_COM_bootstrap_CI.csv"
    df_stats.to_csv(stat_csv, index=False, float_format='%.6f')
    print(f"Statistics CSV saved to {stat_csv}")

# Export Excel workbook with two sheets
excel_path = OUTPUT_DIR / "Fig7C_COM_distance.xlsx"
try:
    with pd.ExcelWriter(excel_path, engine='openpyxl') as writer:
        if ts_dfs:
            df_ts_wide.to_excel(writer, sheet_name='Time_Series', index=False)
        if means:
            df_stats.to_excel(writer, sheet_name='Statistics', index=False)
    print(f"Excel saved to {excel_path}")
except ImportError:
    print("pandas or openpyxl not installed; Excel export skipped.")
except Exception as e:
    print(f"Excel export failed: {e}")

print("Done.")