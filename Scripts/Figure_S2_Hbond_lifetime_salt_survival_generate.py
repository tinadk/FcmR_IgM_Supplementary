#!/usr/bin/env python3
"""
Figure S2: Hydrogen bond and salt bridge survival for WT and three mutants.
Uses shared config.py.  Raw ACF data from XVG files (no smoothing).
Exports:
  - figure (png, jpg, pdf)
  - combined raw survival curves (long-format CSV)
  - separate wide-format CSV for each interaction type
  - combined wide-format Excel with sheets:
      A_Hbond_survival, B_Saltbridge_survival, C_Lifetimes_tau
All Outputs saved to output/Figure_S2_Hbond_salt_survival/
"""

import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from config import DATA_DIR, OUT_DIR, SYSTEMS, COLORS, DPI, MPL_RCPARAMS

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.integrate import simpson
import warnings

warnings.filterwarnings('ignore')
plt.rcParams.update(MPL_RCPARAMS)

# ------------------- Output directory -------------------
OUTPUT_DIR = OUT_DIR / "Figure_S2_Hbond_salt_survival"
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

# ------------------- Figure-specific settings -------------------
LINE_WIDTH = 1.5
ALPHA = 1.0
FONT_SIZE = 9
LEGEND_FONT_SIZE = 9
TICK_FONT_SIZE = 9

YLABEL      = 'Survival probability'
SUPTITLE    = ''

OUTPUT_BASENAME = 'Figure_S2_Hbond_salt_survival'
OUTPUT_CSV_ALL   = OUTPUT_DIR / 'Figure_S2_survival_long_format.csv'

OUTPUT_HBOND_CSV   = OUTPUT_DIR / 'Figure_S2_Hbond_survival_wide.csv'
OUTPUT_SALT_CSV    = OUTPUT_DIR / 'Figure_S2_Saltbridge_survival_wide.csv'
OUTPUT_COMBINED_EXCEL = OUTPUT_DIR / 'Figure_S2_saltbridge_survival_hbond_lifetime_data.xlsx'

# ------------------- Data loading -------------------
def load_hbond_acf(sys):
    f = DATA_DIR / sys / f"{sys}_hbond_acf.xvg"
    if not f.exists():
        return None
    try:
        data = np.loadtxt(f, comments=['@', '#'])
    except Exception as e:
        print(f"Error loading {f}: {e}")
        return None
    if data.shape[1] < 2:
        return None
    time_ns = data[:, 0] / 1000.0
    survival = data[:, 1]
    return time_ns, survival

def load_saltbridge_xvg(sys):
    f = DATA_DIR / sys / f"{sys}_saltbridge_survival.xvg"
    if not f.exists():
        return None
    try:
        data = np.loadtxt(f, comments=['@', '#'])
    except Exception as e:
        print(f"Error loading {f}: {e}")
        return None
    if data.shape[1] < 2:
        return None
    time_ns = data[:, 0] / 1000.0
    survival = data[:, 1]
    return time_ns, survival

def compute_tau(time_ns, survival):
    if len(time_ns) > 1 and np.max(survival) > 0.05:
        return simpson(survival, time_ns)
    return 0.0

# ------------------- Wide-format DataFrame builder -------------------
def build_wide_df(loader_func):
    series_dict = {}
    for sys in SYSTEMS:
        data = loader_func(sys)
        if data is not None:
            t, s = data
            series_dict[sys] = pd.Series(s, index=t)
    if not series_dict:
        return None
    df = pd.DataFrame(series_dict)
    df.index.name = 'Time_ns'
    df = df.reset_index()
    cols = ['Time_ns'] + [s for s in SYSTEMS if s in df.columns]
    return df[cols]

# ------------------- Plotting -------------------
def plot_curves(ax, loader_func, ylabel):
    taus = {}
    for sys in SYSTEMS:
        data = loader_func(sys)
        if data is None:
            continue
        t, s = data
        tau = compute_tau(t, s)
        taus[sys] = tau
        ax.plot(t, s,
                color=COLORS[sys],
                linestyle='-',
                linewidth=LINE_WIDTH,
                alpha=ALPHA,
                label=f"{sys}  tau={tau:.2f} ns")
    ax.set_xlabel('Lag time (ns)', fontsize=FONT_SIZE)
    ax.set_ylabel(ylabel, fontsize=FONT_SIZE)
    ax.tick_params(axis='both', labelsize=TICK_FONT_SIZE)
    ax.grid(False)
    ax.legend(loc='upper right', fontsize=LEGEND_FONT_SIZE)
    ax.set_ylim(-0.05, 1.05)
    ax.margins(x=0.05)
    return taus

# ------------------- Long-format CSV -------------------
def export_raw_curves():
    with open(OUTPUT_CSV_ALL, 'w') as f:
        f.write('Interaction,System,Time_ns,Survival\n')
        for sys in SYSTEMS:
            data = load_hbond_acf(sys)
            if data is not None:
                t, s = data
                for ti, si in zip(t, s):
                    f.write(f'H-bond,{sys},{ti:.6f},{si:.8f}\n')
        for sys in SYSTEMS:
            data = load_saltbridge_xvg(sys)
            if data is not None:
                t, s = data
                for ti, si in zip(t, s):
                    f.write(f'Salt bridge,{sys},{ti:.6f},{si:.8f}\n')
    print(f"Long-format survival curves saved to {OUTPUT_CSV_ALL}")

# ------------------- Wide-format exports -------------------
def export_wide_tables():
    # Build DataFrames
    df_hb = build_wide_df(load_hbond_acf)
    df_sb = build_wide_df(load_saltbridge_xvg)

    # Export CSVs
    if df_hb is not None:
        df_hb.to_csv(OUTPUT_HBOND_CSV, index=False, float_format='%.8f')
        print(f"H-bond wide CSV saved to {OUTPUT_HBOND_CSV}")
    if df_sb is not None:
        df_sb.to_csv(OUTPUT_SALT_CSV, index=False, float_format='%.8f')
        print(f"Salt bridge wide CSV saved to {OUTPUT_SALT_CSV}")

# ------------------- Main -------------------
def main():
    hbond_exists = any(load_hbond_acf(s) is not None for s in SYSTEMS)
    salt_exists  = any(load_saltbridge_xvg(s) is not None for s in SYSTEMS)

    if not hbond_exists and not salt_exists:
        print(f"No data found in {DATA_DIR}.")
        return

    n_plots = int(hbond_exists) + int(salt_exists)
    fig, axes = plt.subplots(1, n_plots, figsize=(12, 5))
    plt.subplots_adjust(wspace=0.25, left=0.06, right=0.98, top=0.90, bottom=0.1)
    if n_plots == 1:
        axes = [axes]

    labels = ['A', 'B']
    for idx, ax in enumerate(axes):
        ax.text(-0.12, 1.06, labels[idx], transform=ax.transAxes,
                fontweight='bold', fontsize=9, va='bottom', ha='center', color='black')

    all_taus = {}
    idx = 0
    df_hb = None
    df_sb = None

    if hbond_exists:
        taus_hb = plot_curves(axes[idx], load_hbond_acf, YLABEL)
        all_taus['H-bond'] = taus_hb
        df_hb = build_wide_df(load_hbond_acf)
        idx += 1
    if salt_exists:
        taus_sb = plot_curves(axes[idx], load_saltbridge_xvg, YLABEL)
        all_taus['Salt bridge'] = taus_sb
        df_sb = build_wide_df(load_saltbridge_xvg)

    # Save figure
    for fmt in ('png', 'jpg', 'pdf'):
        out_path = OUTPUT_DIR / f"{OUTPUT_BASENAME}.{fmt}"
        plt.savefig(out_path, dpi=DPI, facecolor='white')
        print(f"Saved: {out_path}")
    plt.close(fig)

    # Export CSVs
    export_raw_curves()
    export_wide_tables()

    # ----- Combined Excel with all three sheets -----
    with pd.ExcelWriter(OUTPUT_COMBINED_EXCEL, engine='openpyxl') as writer:
        if df_hb is not None:
            df_hb.to_excel(writer, sheet_name='A_Hbond_survival', index=False)
        if df_sb is not None:
            df_sb.to_excel(writer, sheet_name='B_Saltbridge_survival', index=False)
        if all_taus:
            df_tau = pd.DataFrame(all_taus, index=SYSTEMS).T
            df_tau.index.name = 'Interaction'
            # Reset index to make 'Interaction' a column for clarity
            df_tau.reset_index(inplace=True)
            df_tau.to_excel(writer, sheet_name='C_Lifetimes_tau', index=False)

    print(f"Combined Excel with sheets (A_Hbond_survival, B_Saltbridge_survival, C_Lifetimes_tau) saved to {OUTPUT_COMBINED_EXCEL}")

if __name__ == "__main__":
    main()