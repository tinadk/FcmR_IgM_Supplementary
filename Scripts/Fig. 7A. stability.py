import sys; from pathlib import Path
ROOT = Path(__file__).resolve().parents[1]; sys.path.insert(0, str(ROOT))
import numpy as np; import matplotlib.pyplot as plt; import pandas as pd
from config import SYSTEMS, COLORS, DATA_DIR, OUT_DIR, DPI, TIME_MAX_NS
from core.compute import trim_equilibration; from core.stats import bootstrap_mean_ci
from definitions import EQUILIBRATION_NS

# =================== Figure: 1 row × 3 panels ===================
fig, axes = plt.subplots(1, 3, figsize=(24, 6))

# ---- RMSD ----
ax = axes[0]
for sys, c in COLORS.items():
    f = DATA_DIR / sys / "rmsd.xvg"
    if not f.exists(): continue
    d = np.loadtxt(f, comments=["@","#"]); t = d[:,0]/1000.0; y = d[:,1]
    m = t <= TIME_MAX_NS; t, y = t[m], y[m]
    t_eq, y_eq = trim_equilibration(t, y)
    y_smooth = pd.Series(y).rolling(50, center=True, min_periods=1).mean()
    ax.plot(t, y_smooth, color=c, linewidth=1.8, alpha=0.9, label=sys)
ax.axvline(EQUILIBRATION_NS, color='gray', linestyle='--', alpha=0.5, linewidth=1)
ax.set(xlabel="Time (ns)", ylabel="RMSD (nm)", title="Backbone RMSD")
ax.legend(loc='lower right', fontsize=10)
ax.set_ylim(0, None)

# ---- Rg ----
ax = axes[1]
for sys, c in COLORS.items():
    f = DATA_DIR / sys / "rg.xvg"
    if not f.exists(): continue
    d = np.loadtxt(f, comments=["@","#"]); t = d[:,0]/1000.0; y = d[:,1]
    m = t <= TIME_MAX_NS; t, y = t[m], y[m]
    t_eq, y_eq = trim_equilibration(t, y)
    y_smooth = pd.Series(y).rolling(50, center=True, min_periods=1).mean()
    ax.plot(t, y_smooth, color=c, linewidth=1.8, alpha=0.9, label=sys)
ax.axvline(EQUILIBRATION_NS, color='gray', linestyle='--', alpha=0.5, linewidth=1)
ax.set(xlabel="Time (ns)", ylabel="Rg (nm)", title="Radius of Gyration")
ax.legend(loc='lower right', fontsize=10)
ax.set_ylim(0, None)

# ---- SASA ----
ax = axes[2]
for sys, c in COLORS.items():
    f = DATA_DIR / sys / "sasa.xvg"
    if not f.exists(): continue
    d = np.loadtxt(f, comments=["@","#"]); t = d[:,0]/1000.0; y = d[:,1]
    m = t <= TIME_MAX_NS; t, y = t[m], y[m]
    t_eq, y_eq = trim_equilibration(t, y)
    y_smooth = pd.Series(y).rolling(50, center=True, min_periods=1).mean()
    ax.plot(t, y_smooth, color=c, linewidth=1.8, alpha=0.9, label=sys)
ax.axvline(EQUILIBRATION_NS, color='gray', linestyle='--', alpha=0.5, linewidth=1)
ax.set(xlabel="Time (ns)", ylabel="SASA (nm²)", title="Solvent-Accessible Surface Area")
ax.legend(loc='lower right', fontsize=10)
ax.set_ylim(110, None)
ymin, ymax = ax.get_ylim()
ax.set_ylim(ymin, ymax * 1.05)

plt.tight_layout()
plt.savefig(OUT_DIR/"S30_rmsd_rg_sasa_combine.png", dpi=DPI, bbox_inches='tight', facecolor='white')
plt.savefig(OUT_DIR/"S30_rmsd_rg_sasa_combine.jpg", dpi=DPI, bbox_inches='tight', facecolor='white')
plt.savefig(OUT_DIR/"S30_rmsd_rg_sasa_combine.pdf", dpi=DPI, bbox_inches='tight', facecolor='white')
plt.close()
print("Combined figure saved.")