import pandas as pd
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt

from plot_survival_difference import (plot_survival_frequency,
        plot_homotypic_viability_curve, plot_lower_curve, plot_upper_curve,
        plot_rho_curve, plot_point)
from formatting import (cbar_add_central_tick, set_mpl_customisation)
from sys import argv

# Process arguments
if len(argv) != 6:
    print('Usage: {} DATA_WM_I_CSV DATA_WM_II_CSV '
    'DATA_VERTEX_I_CSV DATA_VERTEX_II_CSV IMAGE_PNG'.format(argv[0]))
    exit(1)

data_wm_I_csv = argv[1]
data_wm_II_csv = argv[2]
data_vertex_I_csv = argv[3]
data_vertex_II_csv = argv[4]
image_png = argv[5]

# Customise matplotlib
set_mpl_customisation()

# Read datafile
dfwmI = pd.read_csv(data_wm_I_csv)
dfwmII = pd.read_csv(data_wm_II_csv)

dfwms = [dfwmI, dfwmII]

dfvertexI = pd.read_csv(data_vertex_I_csv)
dfvertexII = pd.read_csv(data_vertex_II_csv)

dfvertexs = [dfvertexI, dfvertexII]

# Vertex data preprocessing
for df in dfvertexs:

    df['num_divisions'] = df['final_cell_count'] - df['initial_cell_count'] + df['num_apoptosis'] + df['num_extrusion']

    df['n0'] = df['initial_cell_count']

    df['num_deaths'] = df['num_apoptosis'] + df['num_extrusion']
    #df['num_deaths'] = df['num_apoptosis']

    df['n0_B'] = df['n0']
    df['num_divisions_B'] = df['num_divisions']
    df['num_deaths_B'] = df['num_deaths']

# Compute survival frequencies
for df in dfwms + dfvertexs:
    df['theta_B'] = df['num_divisions_B'] / (df['num_divisions_B'] + df['num_deaths_B'])

# Check unique eta_A and beta_A values
for df in dfwms + dfvertexs:
    assert len(df['eta_A'].unique()) == 1
    assert len(df['beta_A'].unique()) == 1

# Get unique rho values
rho_wms = []
beta_B_num_wms = []
eta_B_num_wms = []
for df in dfwms:
    rho_B = df['rho_B'].unique()

    assert len(rho_B) == 1

    rho_wms.append(rho_B[0])

    beta_B_num_wms.append(len(df['beta_B'].unique()))
    eta_B_num_wms.append(len(df['eta_B'].unique()))

beta_B_num_max_wm = np.max(beta_B_num_wms)
eta_B_num_max_wm = np.max(eta_B_num_wms)
assert np.all(eta_B_num_wms == eta_B_num_max_wm)

long_beta_B_array_wm = dfwms[0]['beta_B'].unique()

rho_vertexs = []
beta_B_num_vertexs = []
eta_B_num_vertexs = []
for df in dfvertexs:
    rho_B = df['rho_B'].unique()

    assert len(rho_B) == 1

    rho_vertexs.append(rho_B[0])

    beta_B_num_vertexs.append(len(df['beta_B'].unique()))
    eta_B_num_vertexs.append(len(df['eta_B'].unique()))

beta_B_num_max_vertex = np.max(beta_B_num_vertexs)
eta_B_num_max_vertex = np.max(eta_B_num_vertexs)
assert np.all(eta_B_num_vertexs == eta_B_num_max_vertex)

long_beta_B_array_vertex = dfvertexs[0]['beta_B'].unique()

# Create subfigures
fig, axes = plt.subplots(2, 2, sharex='col', sharey='col')

axes = np.transpose(axes)

#from IPython import embed; embed()

axes_wm = list(axes[0])
axes_vertex = list(axes[1])

# Loop over well-mixed cross sections
for df, ax, rho, beta_B_num, eta_B_num, cross_section, cell_based in zip(
        dfwms + dfvertexs,
        axes_wm + axes_vertex,
        rho_wms + rho_vertexs,
        beta_B_num_wms + beta_B_num_vertexs,
        eta_B_num_wms + eta_B_num_vertexs,
        ['I', 'II'] * 2,
        ['well-mixed'] * 2 + ['vertex'] * 2
        ):

    # Set title
    ax.set_title('$\\textrm{{Cross Section {}: {}}}$'.format(cross_section, cell_based))

    # Define extent
    if cell_based == 'well-mixed':
        beta_B_num_max = beta_B_num_max_wm
        eta_B_num_max = eta_B_num_max_wm
        long_beta_B_array = long_beta_B_array_wm
    elif cell_based == 'vertex':
        beta_B_num_max = beta_B_num_max_vertex
        eta_B_num_max = eta_B_num_max_vertex
        long_beta_B_array = long_beta_B_array_vertex

    extent = (
            (beta_B_num_max - beta_B_num) - 0.5,
            beta_B_num_max - 0.5,
            eta_B_num_max - 0.5,
            -0.5)

    # Plot survival frequency
    im = plot_survival_frequency(ax, df, extent=extent)

    # Plot curves
    plot_homotypic_viability_curve(ax, df, coef=1, beta_B_array=long_beta_B_array)

    plot_lower_curve(ax, df, rho, beta_B_array=long_beta_B_array)

    plot_upper_curve(ax, df, rho, beta_B_array=long_beta_B_array)

    plot_rho_curve(ax, df, rho, beta_B_array=long_beta_B_array)

    # Plot D and E simulations
    if cell_based == 'well-mixed':
        if cross_section == 'I':
            plot_point(ax, df, 0.5, 0.16, 'gd', long_beta_B_array)
        elif cross_section == 'II':
            plot_point(ax, df, 0.5, 0.12, 'g^', long_beta_B_array)
            plot_point(ax, df, 0.5, 0.08, 'go', long_beta_B_array)

    # Formatting
    ax.set_xlabel(r'$\beta$')
    ax.set_ylabel(r'$\eta$')

# Add colourbar
fig.subplots_adjust(right=0.8)
cbar_ax = fig.add_axes([0.85, 0.15, 0.05, 0.7])
cbar = fig.colorbar(im, cax=cbar_ax)
cbar_add_central_tick(cbar)

cbar_ax.set_ylabel(r'$\overline{\lambda}: \textrm{Estimated survival frequency}$')

# Formatting
fig.set_figheight(8)
fig.set_figwidth(8)

# Well-mixed model formatting
ax_wm = axes_wm[0]

eta_y = np.array([0.02, 0.10, 0.20, 0.30, 0.40, 0.50])
#eta_y = np.array([0.02, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45, 0.50])
yticks = eta_B_num_max_wm - eta_y / 0.02
yticklabels = [r'${:.2f}$'.format(elem) for elem in eta_y]

ax_wm.set_yticks(yticks)
ax_wm.set_yticklabels(yticklabels)

beta_x = np.array([0.1, 0.3, 0.5, 0.7, 0.95])
xticks = (beta_x - 0.1) / 0.05
xticklabels = [r'${:.2f}$'.format(elem) for elem in beta_x]

ax_wm.set_xticks(xticks)
ax_wm.set_xticklabels(xticklabels)

# Vertex model formatting
ax_vertex = axes_vertex[0]

eta_y = np.array([0.04, 0.10, 0.20, 0.30, 0.40, 0.48])
#eta_y = np.array([0.04, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45, 0.48])
yticks = eta_B_num_max_vertex - eta_y / 0.04
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
