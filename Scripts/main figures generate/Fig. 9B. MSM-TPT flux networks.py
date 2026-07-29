#!/usr/bin/env python3
"""
Fig. 9B – TPT net flux networks for all four systems.
Node area ∝ population fraction (consistent with Fig. 9A).
Arrow thickness ∝ net flux. Edges with flux ≤ 1% of maximum are omitted.
Exports: 2×2 figure, node population matrix, edge fluxes CSV, combined XLSX.
Uses config.py for settings.
"""

import sys
from pathlib import Path
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from config import SYSTEMS, OUT_DIR, DPI, MPL_RCPARAMS

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import warnings

warnings.filterwarnings('ignore')
plt.rcParams.update(MPL_RCPARAMS)

# ---------- Paths ----------
MSM_DIR = ROOT / 'MSM_final' / 'archive' / 'complete' / 'data' / 'msm'
TICA_DIR = ROOT / 'MSM_final' / 'archive' / 'complete' / 'data' / 'tica'

OUTPUT_DIR = OUT_DIR / 'fig9b_tpt_flux'
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

# ---------- Load common data ----------
pcca_membership = np.load(MSM_DIR / 'pcca_membership.npy')
dtraj_all = np.load(MSM_DIR / 'dtraj_k150.npy')
system_labels = np.load(TICA_DIR / 'system_labels.npy')
macro_assign = pcca_membership[dtraj_all]
n_macro = pcca_membership.max() + 1
macro_labels = [f'M{i+1}' for i in range(n_macro)]

sys_to_id = {name: i for i, name in enumerate(SYSTEMS)}

# ---------- Compute population fractions (exactly Fig. 9A data) ----------
pop_fractions = {}
for si, sys_name in enumerate(SYSTEMS):
    mask = system_labels == si
    counts = np.bincount(macro_assign[mask], minlength=n_macro)
    pop_fractions[sys_name] = counts / counts.sum()

# ---------- Helper: build macrostate transition matrix (adjacent frames) ----------
def build_macro_P(sid):
    mask = system_labels == sid
    traj = macro_assign[mask]
    C = np.zeros((n_macro, n_macro))
    for t in range(len(traj) - 1):
        C[traj[t], traj[t+1]] += 1
    row_sum = C.sum(axis=1)
    P = np.full((n_macro, n_macro), np.nan)
    valid = row_sum > 0
    P[valid] = C[valid] / row_sum[valid, None]
    return P

# ---------- Plot & collect data ----------
try:
    import networkx as nx
    has_nx = True
except ImportError:
    has_nx = False
    print("networkx not installed, skipping.")

if has_nx:
    fig, axes = plt.subplots(2, 2, figsize=(18, 14))
    axes = axes.flatten()

    all_edge_data = {}

    for ax, sys_name in zip(axes, SYSTEMS):
        pi = pop_fractions[sys_name]
        P = build_macro_P(sys_to_id[sys_name])

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

        edge_rows = [(f'M{u+1}', f'M{v+1}', f) for u, v, f in G.edges(data='flux')]
        edge_df = pd.DataFrame(edge_rows, columns=['Source', 'Target', 'Net_Flux'])
        all_edge_data[sys_name] = edge_df

        pos = nx.circular_layout(G)
        node_size = [max(200, G.nodes[i]['weight'] * 5000) for i in G.nodes()]
        nx.draw_networkx_nodes(G, pos, node_size=node_size,
                               node_color='lightblue', ax=ax)

        fluxes = [G[u][v]['flux'] for u, v in G.edges()]
        if fluxes:
            max_flux = max(fluxes)
            edges_to_draw = [(u, v) for u, v in G.edges() if G[u][v]['flux'] / max_flux > 0.01]
            widths = [max(0.5, 3 * G[u][v]['flux'] / max_flux) for u, v in edges_to_draw]
            nx.draw_networkx_edges(G, pos, edgelist=edges_to_draw,
                                   width=widths, edge_color='gray', alpha=0.6, ax=ax)

        labels = {i: f'M{i+1}' for i in G.nodes()}
        nx.draw_networkx_labels(G, pos, labels=labels, font_size=10, ax=ax)
        ax.set_title(sys_name)
        ax.axis('off')

    plt.tight_layout()

    for ext in ('png', 'jpg', 'pdf'):
        plt.savefig(OUTPUT_DIR / f'Fig9B_TPT_all_systems.{ext}',
                    dpi=DPI, bbox_inches='tight', facecolor='white')
    plt.close()

    # Export data
    node_matrix = pd.DataFrame(pop_fractions, index=macro_labels)
    node_matrix.index.name = 'Macrostate'
    node_matrix = node_matrix[list(SYSTEMS)]
    node_matrix_pct = node_matrix * 100

    node_matrix.to_csv(OUTPUT_DIR / 'Fig9B_TPT_node_population_matrix.csv', float_format='%.10f')
    node_matrix_pct.to_csv(OUTPUT_DIR / 'Fig9B_TPT_node_population_percent.csv', float_format='%.4f')

    node_long = node_matrix.reset_index().melt(id_vars='Macrostate', var_name='System', value_name='Population_Fraction')
    node_long.to_csv(OUTPUT_DIR / 'Fig9B_TPT_node_population_long.csv', index=False)

    edge_combined = pd.concat([df.assign(System=sys) for sys, df in all_edge_data.items()], ignore_index=True)[['System','Source','Target','Net_Flux']]
    edge_combined.to_csv(OUTPUT_DIR / 'Fig9B_TPT_edge_fluxes.csv', index=False)

    with pd.ExcelWriter(OUTPUT_DIR / 'Fig9B_TPT_data.xlsx', engine='openpyxl') as writer:
        node_matrix.to_excel(writer, sheet_name='Node_Populations')
        node_matrix_pct.to_excel(writer, sheet_name='Node_Populations_Percent')
        node_long.to_excel(writer, sheet_name='Node_Populations_Long', index=False)
        edge_combined.to_excel(writer, sheet_name='Edge_Fluxes', index=False)

    print(f'All Fig. 9B outputs saved to {OUTPUT_DIR}')
else:
    print("networkx is required.")