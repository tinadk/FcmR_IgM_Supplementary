#!/usr/bin/env python3
"""
Fig. 9 MSM analysis - main manuscript figures.
A) Macrostate populations bar chart
B) TPT net flux networks (WT vs T60A_S63A)
C) MFPT matrices (WT vs T60A_S63A)
"""

import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
import warnings

warnings.filterwarnings('ignore')

# ================= CONFIGURATION =================
BASE_MSM = Path('MSM_rebuild_final_part3/LEGACY_ARCHIVE/ARCHIVE_COMPLETE')
MSM_DIR = BASE_MSM / 'data/msm'
TICA_DIR = BASE_MSM / 'data/tica'
OUT_DIR = Path('use_MSM_publication_figures')
OUT_DIR.mkdir(exist_ok=True)

DPI = 2400
LAG_TIME = 10          # frames per step
FRAME_TO_NS = 0.1      # 100 ps per frame -> 1 frame = 0.1 ns
COLORS = {
    'WT': '#2C3E50',
    'D111A': '#E74C3C',
    'T110A_D111A': '#27AE60',
    'T60A_S63A': '#2980B9'
}

# ================= LOAD COMMON DATA =================
print('Loading data...')
pcca_membership = np.load(MSM_DIR / 'pcca_membership.npy')
macro_pops_global = np.load(MSM_DIR / 'pcca_macro_populations.npy')
dtraj_all = np.load(MSM_DIR / 'dtraj_k150.npy')
system_labels = np.load(TICA_DIR / 'system_labels.npy')
macro_assign = pcca_membership[dtraj_all]
n_macro = pcca_membership.max() + 1
macro_labels = [f'M{i+1}' for i in range(n_macro)]

sys_names = ['WT', 'D111A', 'T110A_D111A', 'T60A_S63A']

# ================= A. MACROSTATE POPULATIONS =================
print('A. Macrostate population bar chart...')
sys_counts = []
for si in range(4):
    mask = system_labels == si
    counts = np.bincount(macro_assign[mask], minlength=n_macro)
    sys_counts.append(counts / counts.sum() * 100)

fig, ax = plt.subplots(figsize=(8, 6))
bar_width = 0.15
x = np.arange(n_macro)
all_bars = []

for si in range(4):
    bars = ax.bar(x + si * bar_width, sys_counts[si], bar_width,
                  color=COLORS[sys_names[si]], label=sys_names[si])
    all_bars.append(bars)

for state_idx in range(n_macro):
    heights = [sys_counts[si][state_idx] for si in range(4)]
    max_idx = np.argmax(heights)
    for si in range(4):
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
ax.set_xticklabels([f'M{i+1}' for i in range(n_macro)])
ax.set_ylabel('Proportion (%)')
ax.set_title('Macrostate populations per system')
ax.legend(loc='upper right', fontsize=9)
save_fig(fig, 'Fig9A_Macrostate_populations')

# ================= HELPER: PER-SYSTEM MACRO TRANSITION MATRIX =================


# ================= B. TPT NET FLUX NETWORKS (WT vs T60A_S63A) =================
print('B. TPT net flux networks for WT and T60A_S63A...')
try:
    import networkx as nx
    has_nx = True
except ImportError:
    has_nx = False
    print("networkx not installed, skipping TPT networks.")

if has_nx:
    target_systems = {'WT': 0, 'T60A_S63A': 3}
    fig, axes = plt.subplots(1, 2, figsize=(16, 7))
    for ax, (name, sys_id) in zip(axes, target_systems.items()):
        P = build_macro_P(sys_id)
        pi = stationary_from_P(P)
        G = nx.DiGraph()
        for i in range(n_macro):
            if pi[i] > 1e-6:
                G.add_node(i, weight=pi[i])
        for i in range(n_macro):
            for j in range(n_macro):
                if i != j and not np.isnan(P[i, j]) and pi[i] > 1e-6 and pi[j] > 1e-6:
                    flux = pi[i] * P[i, j]
                    if flux > 1e-6:
                        G.add_edge(i, j, flux=flux)

        pos = nx.circular_layout(G)
        node_size = [G.nodes[i]['weight'] * 5000 for i in G.nodes()]
        nx.draw_networkx_nodes(G, pos, node_size=node_size,
                               node_color='lightblue', ax=ax)
        fluxes = [G[u][v]['flux'] for u, v in G.edges()]
        if fluxes:
            max_flux = max(fluxes)
            nx.draw_networkx_edges(
                G, pos, edgelist=G.edges(),
                width=[3 * f / max_flux for f in fluxes],
                edge_color='gray', alpha=0.6, ax=ax)
        labels = {i: f'M{i+1}' for i in G.nodes()}
        nx.draw_networkx_labels(G, pos, labels=labels, font_size=10, ax=ax)
        ax.set_title(name)
        ax.axis('off')
    plt.suptitle('TPT Net Flux Networks', fontweight='bold')
    plt.tight_layout()
    save_fig(fig, 'Fig9B_TPT_flux_networks')

# ================= C. MFPT MATRICES (WT vs T60A_S63A) =================
print('C. MFPT matrices for WT and T60A_S63A...')

mfpt_wt = mfpt_for_system(0)
mfpt_t60a = mfpt_for_system(3)

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))
plot_single_mfpt(ax1, mfpt_wt, 'WT')
plot_single_mfpt(ax2, mfpt_t60a, 'T60A_S63A')
plt.suptitle('Mean First Passage Time matrices', fontweight='bold')
plt.tight_layout()
save_fig(fig, 'Fig9C_MFPT_matrices')

print('All main Fig. 9 panels saved in', OUT_DIR)