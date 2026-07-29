#!/usr/bin/env python3
"""
Table S2 – MSM validation key parameters.
Uses only the last two CK error columns (Set2 & Set3).
Outputs: output/tables/Table S2. MSM validation.xlsx
"""

import sys
from pathlib import Path
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from config import OUT_DIR
import numpy as np
import pandas as pd

MSM_DIR = ROOT / 'MSM_final' / 'archive' / 'complete' / 'data' / 'msm'
OUTPUT_DIR = OUT_DIR / 'tables'
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

# Load data
its_frames = np.load(MSM_DIR / 'its_k150.npy')
lags_frames = np.load(MSM_DIR / 'lags_k150.npy')
ck_raw = np.load(MSM_DIR / 'ck_errors.npy')
eig = np.real(np.load(MSM_DIR / 'eigenvalues.npy'))

LAG_NS = 1.0
lags_ns = lags_frames * 0.1
its_ns = its_frames * 0.1
lag_idx = np.argmin(np.abs(lags_ns - LAG_NS))

# Implied timescales (first 5)
its_vals = its_ns[lag_idx, :5]

# CK error: ignore first column (counter), use Set2 & Set3
ck_errors = ck_raw[:, 1:] if ck_raw.shape[1] >= 3 else ck_raw
ck_avg = np.mean(ck_errors[lag_idx, :])

# Eigenvalue gap after 5th eigenvalue
spec = eig[eig < 0.9999]
eig_gap = spec[4] - spec[5] if len(spec) > 5 else np.nan

# Build and save table
data = {
    'Parameter': [
        'Lag time (ns)',
        'Number of macrostates',
        'ITS1 (ns)', 'ITS2 (ns)', 'ITS3 (ns)', 'ITS4 (ns)', 'ITS5 (ns)',
        'Average CK error',
        'Eigenvalue gap (λ5‑λ6)'
    ],
    'Value': [
        LAG_NS,
        5,
        round(its_vals[0], 2), round(its_vals[1], 2), round(its_vals[2], 2),
        round(its_vals[3], 2), round(its_vals[4], 2),
        round(ck_avg, 4),
        round(eig_gap, 4)
    ]
}

df = pd.DataFrame(data)
xlsx_path = OUTPUT_DIR / 'Table S2. MSM validation.xlsx'
df.to_excel(xlsx_path, index=False)
print(f'Table S2 saved to {xlsx_path}')