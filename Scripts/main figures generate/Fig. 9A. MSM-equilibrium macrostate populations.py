#!/usr/bin/env python3
"""
Fig. 9A – Macrostate populations bar chart.
Data from MSM PCCA.  Uses config.py for systems, colors, output, DPI.
No top title, x‑axis labeled 'Macrostate', y‑axis 'Population (%)'.
"""

import sys
from pathlib import Path
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from config import SYSTEMS, COLORS, OUT_DIR, DPI, MPL_RCPARAMS

import numpy as np
import matplotlib.pyplot as plt
import matplotlib
matplotlib.use('Agg')
import warnings

warnings.filterwarnings('ignore')
plt.rcParams.update(MPL_RCPARAMS)

# ---------- Data paths ----------
MSM_DIR = ROOT / 'MSM_final' / 'archive' / 'complete' / 'data' / 'msm'
TICA_DIR = ROOT / 'MSM_final' / 'archive' / 'complete' / 'data' / 'tica'

OUTPUT_DIR = OUT_DIR / 'fig9a_macrostate_populations'
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

# ---------- Load data ----------
pcca_membership = np.load(MSM_DIR / 'pcca_membership.npy')
dtraj_all = np.load(MSM_DIR / 'dtraj_k150.npy')
system_labels = np.load(TICA_DIR / 'system_labels.npy')
macro_assign = pcca_membership[dtraj_all]
n_macro = pcca_membership.max() + 1
macro_labels = [f'M{i+1}' for i in range(n_macro)]

# ---------- Compute per‑system populations (%) ----------
sys_names = list(SYSTEMS)
sys_counts = []
for si in range(len(sys_names)):
    mask = system_labels == si
    counts = np.bincount(macro_assign[mask], minlength=n_macro)
    sys_counts.append(counts / counts.sum() * 100)

# ---------- Plot ----------
fig, ax = plt.subplots(figsize=(8, 6))
bar_width = 0.15
x = np.arange(n_macro)
all_bars = []

color_list = [COLORS[s] for s in sys_names]

for si in range(len(sys_names)):
    bars = ax.bar(x + si * bar_width, sys_counts[si], bar_width,
                  color=color_list[si], label=sys_names[si])
    all_bars.append(bars)

# Value labels (same logic as original)
for state_idx in range(n_macro):
    heights = [sys_counts[si][state_idx] for si in range(len(sys_names))]
    max_idx = np.argmax(heights)
    for si in range(len(sys_names)):
        bar = all_bars[si][state_idx]
        height = heights[si]
        if height == 0:
            continue
        x_center = bar.get_x() + bar.get_width() / 2.0
        shift = 0.0
        if not (state_idx == 0 and si == 1):
            if si < max_idx:
                shift = -0.08
            elif si > max_idx:
                shift = 0.08
        extra_offset = 5.0 if (state_idx == 0 and si == 1) else 0.0
        ax.text(x_center + shift, height + 0.3 + extra_offset,
                f'{height:.1f}%', ha='center', va='bottom', fontsize=6)

max_y = max([max(sc) for sc in sys_counts]) * 1.25
ax.set_ylim(0, max_y)
ax.set_xticks(x + bar_width * 1.5)
ax.set_xticklabels(macro_labels)
ax.set_xlabel('Macrostate')                # added x‑axis label
ax.set_ylabel('Population (%)')            # changed from 'Proportion (%)'
# top title removed
ax.legend(loc='upper right', fontsize=9)

# ---------- Save ----------
for ext in ('png', 'jpg', 'pdf'):
    plt.savefig(OUTPUT_DIR / f'Fig9A_Macrostate_populations.{ext}',
                dpi=DPI, bbox_inches='tight', facecolor='white')
plt.close()

# ---------- Export data table ----------
import pandas as pd
data = {}
for si, sys_name in enumerate(sys_names):
    mask = system_labels == si
    counts = np.bincount(macro_assign[mask], minlength=n_macro)
    proportions = counts / counts.sum() * 100
    data[sys_name] = proportions

df = pd.DataFrame(data, index=macro_labels)
df = df[list(sys_names)]
df.to_csv(OUTPUT_DIR / 'Supplementary Table S3 Macrostate Populations.csv')
df.to_excel(OUTPUT_DIR / 'Supplementary Table S3 Macrostate Populations.xlsx')

print(f'Fig. 9A and table saved to {OUTPUT_DIR}')