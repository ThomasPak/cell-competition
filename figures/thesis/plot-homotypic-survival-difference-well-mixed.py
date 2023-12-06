import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt

from formatting import (set_mpl_customisation, cbar_add_central_tick)
from plot_survival_difference import (plot_homotypic_survival_difference,
        plot_neutral_coexistence_point, plot_neutral_competition_curve,
        plot_neutral_simulations)
from sys import argv

# Process arguments
if len(argv) != 5:
    print('Usage: {} COMPETITION_DATA_CSV CONTROL_DATA_CSV CELL_TYPE IMAGE_PNG'.format(argv[0]))
    exit(1)

competition_data_csv = argv[1]
control_data_csv = argv[2]
cell_type = argv[3]
image_png = argv[4]

assert cell_type == 'A' or cell_type == 'B'

# Customise matplotlib
set_mpl_customisation()

# Read datafile
df_competition = pd.read_csv(competition_data_csv)
df_control = pd.read_csv(control_data_csv)

# Compute survival frequency
df_competition['theta_A'] = df_competition['num_divisions_A'] / (df_competition['num_divisions_A'] + df_competition['num_deaths_A'])
df_competition['theta_B'] = df_competition['num_divisions_B'] / (df_competition['num_divisions_B'] + df_competition['num_deaths_B'])

df_control['theta_A'] = df_control['num_divisions_A'] / (df_control['num_divisions_A'] + df_control['num_deaths_A'])
df_control['theta_B'] = df_control['num_divisions_B'] / (df_control['num_divisions_B'] + df_control['num_deaths_B'])

# Get unique eta_A and beta_A values
unique_eta_As = df_competition['eta_A'].unique()
unique_beta_As = df_competition['beta_A'].unique()

# Create subfigures
fig, axes = plt.subplots(len(unique_eta_As), len(unique_beta_As),
        sharex=True, sharey=True)

# Invert vertical ordering of axes
axes = axes[::-1, :]

# Initialise ims
ims = []
max_abs_val = -np.inf

# Loop over eta_A and beta_A
for i, eta_A in enumerate(unique_eta_As):
    for j, beta_A in enumerate(unique_beta_As):

        # Select current axes
        ax = axes[i,j]

        # Select eta_A and beta_A
        cur_df = df_competition[(df_competition['eta_A'] == eta_A) &
                (df_competition['beta_A'] == beta_A)]

        if cell_type == 'A':

            cur_df_control = df_control[np.isclose(df_control['eta_B'], eta_A) &
                    np.isclose(df_control['beta_B'], beta_A)]

        else:

            cur_df_control = df_control

        # Set title
        ax.set_title(r'$\eta_B = {}, \beta_B = {}$'.format(eta_A, beta_A))

        im, cur_max_abs_val = plot_homotypic_survival_difference(ax, cur_df,
                cur_df_control, beta_A, eta_A, cell_type)

        ims.append(im)
        max_abs_val = max(max_abs_val, cur_max_abs_val)

        plot_neutral_competition_curve(ax, cur_df, beta_A, eta_A)

        plot_neutral_coexistence_point(ax, cur_df, beta_A, eta_A)

        #plot_neutral_simulations(ax, cur_df, cur_df_control, cell_type, 0.001)

# Set clim
for im in ims:
    im.set_clim(-max_abs_val, max_abs_val)

# Add colourbar
fig.subplots_adjust(right=0.8)
cbar_ax = fig.add_axes([0.85, 0.15, 0.05, 0.7])
cbar = fig.colorbar(im, cax=cbar_ax)
cbar_add_central_tick(cbar)

cbar_ax.set_ylabel("$\\overline{{\\Delta^{{=}}_{{{0}|{1}}}}}: \\textrm{{Estimated homotypic survival difference}}$".format('A' if cell_type == 'B' else 'B', cell_type))

# Formatting
fig.set_figheight(10)
fig.set_figwidth(8)

ax = axes[0,0]

yticks_old = ax.get_yticks()
yticklabels_old = ax.get_yticklabels()

ax.set_yticks(yticks_old[::8])
ax.set_yticklabels(yticklabels_old[::8])

xticks_old = ax.get_xticks()
xticklabels_old = ax.get_xticklabels()

ax.set_xticks(xticks_old[::6])
ax.set_xticklabels(xticklabels_old[::6])

eta_B_num_max_vertex = len(df_competition['eta_B'].unique())

eta_y = np.array([0.01, 0.05, 0.10, 0.15, 0.20, 0.25])
yticks = eta_B_num_max_vertex - eta_y / 0.01
yticklabels = [r'${:.2f}$'.format(elem) for elem in eta_y]

ax.set_yticks(yticks)
ax.set_yticklabels(yticklabels)

beta_x = np.array([0.05, 0.3, 0.5, 0.7, 0.95])
xticks = (beta_x - 0.05) / 0.05
xticklabels = [r'${:.2f}$'.format(elem) for elem in beta_x]

ax.set_xticks(xticks)
ax.set_xticklabels(xticklabels)

# Sketch regions for cell type B
if cell_type == 'B':
    ax = axes[1, 1]

    point1 = [0.75, 0.16]
    point2 = [0.95 + 0.05/2, 0.21]
    point3 = [0.95 + 0.05/2, 0.03]
    point4 = [0.5, 0.16]

    x_to_i = lambda x: x / 0.05 - 1
    y_to_j = lambda y: 25 - y / 0.01

    xy_to_ij = lambda point: [x_to_i(point[0]), y_to_j(point[1])]

    point1 = xy_to_ij(point1)
    point2 = xy_to_ij(point2)
    point3 = xy_to_ij(point3)
    point4 = xy_to_ij(point4)

    color = 'grey'
    alpha = 1

    ax.plot([point1[0], point2[0]], [point1[1], point2[1]], color=color, alpha=alpha)
    ax.plot([point1[0], point3[0]], [point1[1], point3[1]], color=color, alpha=alpha)
    ax.plot([point1[0], point4[0]], [point1[1], point4[1]], color=color, alpha=alpha)

# Save figure
fig.savefig(image_png, bbox_inches='tight')
