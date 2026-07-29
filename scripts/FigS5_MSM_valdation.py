#!/usr/bin/env python3
"""
Fig. S5 – MSM validation (implied timescales, CK test, eigenvalue spectrum,
          timescale separation). Single combined 2×2 panel.
Reads data and parameters from config.py.  No titles.
"""

import sys
from pathlib import Path
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from config import OUT_DIR, DPI, MPL_RCPARAMS

import numpy as np
import matplotlib.pyplot as plt

plt.rcParams.update(MPL_RCPARAMS)

# ---------- Paths ----------
MSM_DIR = ROOT / 'MSM_final' / 'archive' / 'complete' / 'data' / 'msm'
OUTPUT_DIR = OUT_DIR / 'Supplementary_figS5_MSM_validation'
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

# ---------- Load data ----------
its = np.load(MSM_DIR / 'its_k150.npy')
lags = np.load(MSM_DIR / 'lags_k150.npy') * 0.1          # frames → ns
its_ns = its * 0.1

ck = np.load(MSM_DIR / 'ck_errors.npy')
lags_ck = (np.arange(1, ck.shape[0] + 1) * 0.1).astype(float)

eig = np.real(np.load(MSM_DIR / 'eigenvalues.npy'))
spec = eig[eig < 0.9999]                                 # discard trivial eigenvalue ~1
t_scale = -1.0 / np.log(spec[:20]) * 0.1                 # first 20 timescales (ns)

# ---------- CSV export ----------
# Implied timescales
its_header = 'lag_ns,' + ','.join(f'ITS{i+1}' for i in range(its_ns.shape[1]))
its_csv = np.column_stack([lags, its_ns])
np.savetxt(OUTPUT_DIR / 'Supplementary_FigS5_ITS.csv', its_csv, delimiter=',',
           header=its_header, comments='', fmt='%.6f')

# CK errors
ck_header = 'lag_ns,' + ','.join(f'CK_error_set{i+1}' for i in range(ck.shape[1]))
ck_csv = np.column_stack([lags_ck, ck])
np.savetxt(OUTPUT_DIR / 'Supplementary_FigS5_CK.csv', ck_csv, delimiter=',',
           header=ck_header, comments='', fmt='%.6f')

# Eigenvalues and timescales (first 20 timescales, NaN padded)
timescale_padded = np.pad(t_scale, (0, len(spec) - len(t_scale)),
                          constant_values=np.nan)
eig_csv = np.column_stack([np.arange(1, len(spec) + 1), spec, timescale_padded])
np.savetxt(OUTPUT_DIR / 'Supplementary_FigS5_eigenvalues.csv', eig_csv, delimiter=',',
           header='index,eigenvalue,timescale_ns', comments='', fmt='%.6f')

# ---------- Combined 2×2 figure ----------
fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(14, 10))

# Top‑left: Implied timescales
for i in range(min(5, its_ns.shape[1])):
    ax1.plot(lags, its_ns[:, i], marker='o', label=f'ITS {i+1}')
ax1.set_xlabel('Lag time (ns)')
ax1.set_ylabel('Timescale (ns)')
ax1.legend(fontsize=8)

# Top‑right: Chapman–Kolmogorov test
for i in range(ck.shape[1]):
    ax2.plot(lags_ck, ck[:, i], marker='o', label=f'Set {i+1}')
ax2.set_xlabel('Lag time (ns)')
ax2.set_ylabel('CK error')
ax2.legend(fontsize=8)

# Bottom‑left: Eigenvalue spectrum
ax3.plot(range(1, len(spec) + 1), spec, 'o-', markersize=4)
ax3.set_xlabel('Eigenvalue index')
ax3.set_ylabel('Eigenvalue (real part)')

# Bottom‑right: Timescale separation
ax4.plot(range(1, len(t_scale) + 1), t_scale, 's-', color='crimson')
ax4.set_xlabel('Eigenvalue index')
ax4.set_ylabel('Implied timescale (ns)')

# No suptitle, no individual subplot titles
plt.tight_layout()

# Save figure
for ext in ('png', 'jpg', 'pdf'):
    plt.savefig(OUTPUT_DIR / f'FigS5_MSM_validation.{ext}',
                dpi=DPI, bbox_inches='tight', facecolor='white')
plt.close()

print(f'Supplementary_Fig. S5 and CSV files saved to {OUTPUT_DIR}')