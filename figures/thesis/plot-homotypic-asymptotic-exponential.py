import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt

from formatting import (set_mpl_customisation, cbar_add_central_tick)
from plot_survival_difference import (plot_survival_frequency,
        plot_homotypic_viability_curve, plot_point)
from sys import argv

# Process arguments
if len(argv) != 4:
    print('Usage: {} DATA_WM_CSV DATA_VERTEX_CSV IMAGE_PNG'.format(argv[0]))
    exit(1)

data_wm_csv = argv[1]
data_vertex_csv = argv[2]
image_png = argv[3]

# Customise matplotlib
set_mpl_customisation()

# Read datafile
df_wm = pd.read_csv(data_wm_csv)
df_vertex = pd.read_csv(data_vertex_csv)

df_vertex['num_divisions'] = df_vertex['final_cell_count'] - df_vertex['initial_cell_count'] + df_vertex['num_apoptosis'] + df_vertex['num_extrusion']

df_vertex['n0'] = df_vertex['initial_cell_count']

df_vertex['num_deaths'] = df_vertex['num_apoptosis'] + df_vertex['num_extrusion']

df_vertex['num_divisions_B'] = df_vertex['num_divisions']
df_vertex['num_deaths_B'] = df_vertex['num_deaths']

# Compute survival frequencies
df_wm['theta_B'] = df_wm['num_divisions_B'] / (df_wm['num_divisions_B'] + df_wm['num_deaths_B'])
df_vertex['theta_B'] = df_vertex['num_divisions_B'] / (df_vertex['num_divisions_B'] + df_vertex['num_deaths_B'])

# Get unique eta_A and beta_A values
assert len(df_wm['eta_A'].unique()) == 1
assert len(df_wm['beta_A'].unique()) == 1
assert len(df_vertex['eta_A'].unique()) == 1
assert len(df_vertex['beta_A'].unique()) == 1

# Create subfigures
fig, axes = plt.subplots(1, 2)

ax_wm = axes[0]
ax_vertex = axes[1]

im = plot_survival_frequency(ax_wm, df_wm)
im = plot_survival_frequency(ax_vertex, df_vertex)

ax_wm.set_xlabel('$\\beta$')
ax_wm.set_ylabel('$\\eta$')

ax_vertex.set_xlabel('$\\beta$')
ax_vertex.set_ylabel('$\\eta$')

plot_homotypic_viability_curve(ax_wm, df_wm)
plot_homotypic_viability_curve(ax_vertex, df_vertex)

# Mark simulations
plot_point(ax_wm, df_wm, 0.5, 0.2, 'gd')
plot_point(ax_wm, df_wm, 0.5, 0.05, 'go')

# Add colourbar
fig.subplots_adjust(right=0.8)
cbar_ax = fig.add_axes([0.85, 0.15, 0.05, 0.7])
cbar = fig.colorbar(im, cax=cbar_ax)
cbar_add_central_tick(cbar)

cbar_ax.set_ylabel(r'$\overline{\lambda}: \textrm{Estimated survival frequency}$')

# Formatting
fig.set_figheight(4)
fig.set_figwidth(8)

# Well-mixed model formatting
ax_wm.set_title('$\\textrm{Well-mixed}$')

eta_B_num_max_wm = len(df_wm['eta_B'].unique())

#eta_y = np.array([0.02, 0.10, 0.20, 0.30, 0.40, 0.50])
eta_y = np.array([0.01, 0.05, 0.10, 0.15, 0.20, 0.25])
yticks = eta_B_num_max_wm - eta_y / 0.01
yticklabels = [r'${:.2f}$'.format(elem) for elem in eta_y]

ax_wm.set_yticks(yticks)
ax_wm.set_yticklabels(yticklabels)

beta_x = np.array([0.05, 0.3, 0.5, 0.7, 0.95])
xticks = (beta_x - 0.05) / 0.05
xticklabels = [r'${:.2f}$'.format(elem) for elem in beta_x]

ax_wm.set_xticks(xticks)
ax_wm.set_xticklabels(xticklabels)

# Vertex model formatting
ax_vertex.set_title('$\\textrm{Vertex}$')

eta_B_num_max_vertex = len(df_vertex['eta_B'].unique())

#eta_y = np.array([0.04, 0.10, 0.20, 0.30, 0.40, 0.48])
eta_y = np.array([0.02, 0.05, 0.10, 0.15, 0.20, 0.24])
yticks = eta_B_num_max_vertex - eta_y / 0.02
yticklabels = [r'${:.2f}$'.format(elem) for elem in eta_y]

ax_vertex.set_yticks(yticks)
ax_vertex.set_yticklabels(yticklabels)

beta_x = np.array([0.1, 0.3, 0.5, 0.7, 0.9])
xticks = (beta_x - 0.10) / 0.10
xticklabels = [r'${:.2f}$'.format(elem) for elem in beta_x]

ax_vertex.set_xticks(xticks)
ax_vertex.set_xticklabels(xticklabels)

# Save figure
#fig.tight_layout()
fig.savefig(image_png, bbox_inches='tight')
