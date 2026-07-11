import sys
from pathlib import Path
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from config import SYSTEMS, DATA_DIR, OUT_DIR, TIME_MAX_NS
from definitions import EQUILIBRATION_NS

plt.rcParams.update({'font.family': 'Arial', 'font.size': 9})

colors = {
    'WT': '#2C3E50',
    'D111A': '#E74C3C',
    'T60A_S63A': '#2980B9',
    'T110A_D111A': '#27AE60'
}

means = {}
cis = {}

# ---------- Figure 1: Time series ----------
fig1, ax1 = plt.subplots(figsize=(8, 5))

for s, c in colors.items():
    f = DATA_DIR / s / "com_distance.xvg"
    if not f.exists():
        continue

    d = np.loadtxt(f, comments=["@", "#"])
    t = d[:, 0] / 1000.0
    y = d[:, 1]

    mask = t <= TIME_MAX_NS
    t, y = t[mask], y[mask]

    eq_mask = t >= EQUILIBRATION_NS
    y_eq = y[eq_mask]

    y_smooth = pd.Series(y).rolling(50, center=True, min_periods=1).mean()
    ax1.plot(t, y_smooth, color=c, linewidth=1.5, alpha=0.85, label=s)

    mean = np.mean(y_eq)
    boot = np.random.choice(y_eq, size=(1000, len(y_eq)), replace=True).mean(axis=1)
    ci_low, ci_high = np.percentile(boot, [2.5, 97.5])

    means[s] = mean
    cis[s] = (ci_low, ci_high)

ax1.axvline(EQUILIBRATION_NS, color='gray', linestyle='--', alpha=0.4)
ax1.set(xlabel="Time (ns)", ylabel="COM Distance (nm)", title="COM Distance (FcμR–IgM)")
ax1.legend(loc='center left', bbox_to_anchor=(0.15, 0.9), frameon=False)
ax1.grid(visible=False)
for spine in ax1.spines.values():
    spine.set_visible(True)

plt.tight_layout()
plt.savefig(OUT_DIR / "S12_COM_distance.png", dpi=1600, facecolor='white', edgecolor='none')
plt.savefig(OUT_DIR / "S12_COM_distance.jpg", dpi=1600, facecolor='white', edgecolor='none')
plt.savefig(OUT_DIR / "S12_COM_distance.pdf", dpi=1600, facecolor='white', edgecolor='none')
plt.close()

# ---------- Figure 2: Bar plot with CI ----------
fig2, ax2 = plt.subplots(figsize=(8, 5))

labels = list(means.keys())
vals = [means[s] for s in labels]
err_low = [means[s] - cis[s][0] for s in labels]
err_high = [cis[s][1] - means[s] for s in labels]

ax2.bar(
    labels,
    vals,
    yerr=[err_low, err_high],
    color=[colors[s] for s in labels],
    alpha=0.85,
    edgecolor='white',
    capsize=4
)

ax2.set_ylabel("Mean COM Distance (nm)")
ax2.set_title("Equilibrium COM Distance (mean ± 95% CI)")
ax2.grid(False)
for spine in ax2.spines.values():
    spine.set_visible(True)

plt.tight_layout()
plt.savefig(OUT_DIR / "S12_Equilibrium_COM_distance_mean.png", dpi=1600, facecolor='white', edgecolor='none')
plt.savefig(OUT_DIR / "S12_Equilibrium_COM_distance_mean.jpg", dpi=1600, facecolor='white', edgecolor='none')
plt.savefig(OUT_DIR / "S12_Equilibrium_COM_distance_mean.pdf", dpi=1600, facecolor='white', edgecolor='none')
plt.close()

print("Two separate output saved: S12_COM_distance and S12_Equilibrium_COM_distance_mean")