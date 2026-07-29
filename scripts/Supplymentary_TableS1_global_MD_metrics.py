#!/usr/bin/env python3
"""
Table S1 – Global MD metrics with bootstrap 95% CI.
Uses core.compute and core.stats for statistics.
Reads data: data/<system>/{rmsd, rg, sasa, hbond}.xvg
Outputs: output/tables/Table S1. global MD metrics.xlsx
"""

import sys
from pathlib import Path
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from config import SYSTEMS, DATA_DIR, OUT_DIR, TIME_MAX_NS
from core.compute import trim_equilibration
from core.stats import bootstrap_mean_ci

import numpy as np
import pandas as pd

# ---------- Output directory ----------
OUTPUT_DIR = OUT_DIR / 'tables'
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

# ---------- File definitions ----------
metrics = {
    'RMSD (nm)': 'rmsd.xvg',
    'Rg (nm)': 'rg.xvg',
    'SASA (nm²)': 'sasa.xvg',
    'Hbonds': 'hbond.xvg'          # number of hydrogen bonds over time
}

rows = []

for sys in SYSTEMS:
    row = [sys]
    for metric_name, fname in metrics.items():
        fpath = DATA_DIR / sys / fname
        if not fpath.exists():
            print(f'  WARNING: {fpath} not found, filling NaN')
            row.extend([np.nan, np.nan, np.nan])
            continue

        # Load, truncate to TIME_MAX_NS, remove equilibration
        data = np.loadtxt(fpath, comments=['@', '#'])
        t = data[:, 0] / 1000.0          # ps → ns
        y = data[:, 1]
        mask = t <= TIME_MAX_NS
        t, y = t[mask], y[mask]
        t_eq, y_eq = trim_equilibration(t, y)

        # Bootstrap statistics
        mean_val, ci_low, ci_high = bootstrap_mean_ci(y_eq)

        row.extend([round(mean_val, 3), round(ci_low, 3), round(ci_high, 3)])
    rows.append(row)

# ---------- Build DataFrame ----------
columns = ['System']
for m in metrics.keys():
    columns += [f'{m} Mean', f'{m} CI_low', f'{m} CI_high']

df = pd.DataFrame(rows, columns=columns)

# ---------- Save Excel ----------
xlsx_path = OUTPUT_DIR / 'Table S1. global MD metrics.xlsx'
df.to_excel(xlsx_path, index=False)
print(f'Table S1 saved to {xlsx_path}')