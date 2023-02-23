import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt

from formatting import (set_mpl_customisation, cbar_add_central_tick)
from plot_survival_difference import (plot_homotypic_survival_difference,
        plot_neutral_coexistence_point, plot_neutral_competition_curve)
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

# Switcharoo
cell_type = 'B' if cell_type == 'A' else 'A'

# Customise matplotlib
set_mpl_customisation('combination')

# Read datafile
df_competition = pd.read_csv(competition_data_csv)
df_control = pd.read_csv(control_data_csv)

# Vertex data preprocessing
df_competition['num_divisions_A'] = df_competition['final_cell_count_A'] - df_competition['initial_cell_count_A'] + df_competition['num_apoptosis_A'] + df_competition['num_extrusion_A']
df_competition['num_divisions_B'] = df_competition['final_cell_count_B'] - df_competition['initial_cell_count_B'] + df_competition['num_apoptosis_B'] + df_competition['num_extrusion_B']

df_competition['n0_A'] = df_competition['initial_cell_count_A']
df_competition['n0_B'] = df_competition['initial_cell_count_B']

df_competition['num_deaths_A'] = df_competition['num_apoptosis_A'] + df_competition['num_extrusion_A']
df_competition['num_deaths_B'] = df_competition['num_apoptosis_B'] + df_competition['num_extrusion_B']

df_control['final_cell_count_A'] = df_control['final_cell_count']
df_control['initial_cell_count_A'] = df_control['initial_cell_count']
df_control['num_apoptosis_A'] = df_control['num_apoptosis']
df_control['num_extrusion_A'] = df_control['num_extrusion']

df_control['final_cell_count_B'] = df_control['final_cell_count']
df_control['initial_cell_count_B'] = df_control['initial_cell_count']
df_control['num_apoptosis_B'] = df_control['num_apoptosis']
df_control['num_extrusion_B'] = df_control['num_extrusion']


df_control['num_divisions_A'] = df_control['final_cell_count_A'] - df_control['initial_cell_count_A'] + df_control['num_apoptosis_A'] + df_control['num_extrusion_A']
df_control['num_divisions_B'] = df_control['final_cell_count_B'] - df_control['initial_cell_count_B'] + df_control['num_apoptosis_B'] + df_control['num_extrusion_B']

df_control['n0_A'] = df_control['initial_cell_count_A']
df_control['n0_B'] = df_control['initial_cell_count_B']

df_control['num_deaths_A'] = df_control['num_apoptosis_A'] + df_control['num_extrusion_A']
df_control['num_deaths_B'] = df_control['num_apoptosis_B'] + df_control['num_extrusion_B']

# Compute survival frequency
df_competition['theta_A'] = df_competition['num_divisions_A'] / (df_competition['num_divisions_A'] + df_competition['num_deaths_A'])
df_competition['theta_B'] = df_competition['num_divisions_B'] / (df_competition['num_divisions_B'] + df_competition['num_deaths_B'])

df_control['theta_A'] = df_control['num_divisions_A'] / (df_control['num_divisions_A'] + df_control['num_deaths_A'])
df_control['theta_B'] = df_control['num_divisions_B'] / (df_control['num_divisions_B'] + df_control['num_deaths_B'])

# Get unique eta_A and beta_A values
unique_eta_As = df_competition['eta_A'].unique()
unique_beta_As = df_competition['beta_A'].unique()
assert len(unique_eta_As) == 1
assert len(unique_beta_As) == 1
eta_A = unique_eta_As[0]
beta_A = unique_beta_As[0]

# Assert that patterns are random and segregated
unique_patterns = sorted(list(df_competition['pattern'].unique()))
assert unique_patterns[0] == 'random'
assert unique_patterns[1] == 'segregated'

unique_patterns_control = sorted(list(df_control['pattern'].unique()))
assert len(unique_patterns_control) == 1
assert unique_patterns_control[0][:7] == 'control'

# Create subfigures
fig, axes = plt.subplots(1, 2, sharex=True, sharey=True)

# Initialise ims
ims = []
max_abs_val = -np.inf

# Loop over patterns
for i, pattern in enumerate(unique_patterns):

    # Select current axes
    ax = axes[i]

    # Select pattern
    cur_df = df_competition[df_competition['pattern'] == pattern]

    if cell_type == 'A':

        cur_df_control = df_control[np.isclose(df_control['eta_B'], eta_A) &
                np.isclose(df_control['beta_B'], beta_A)]

    else:

        cur_df_control = df_control

    # Set title
    ax.set_title("$\\textrm{{{}}}$".format(pattern[0].upper() + pattern[1:]))

    im, cur_max_abs_val = plot_homotypic_survival_difference(ax, cur_df,
            cur_df_control, beta_A, eta_A, cell_type)

    ims.append(im)
    max_abs_val = max(max_abs_val, cur_max_abs_val)

    plot_neutral_competition_curve(ax, cur_df, beta_A, eta_A)

    plot_neutral_coexistence_point(ax, cur_df, beta_A, eta_A)

# Set clim
for im in ims:
    im.set_clim(-max_abs_val, max_abs_val)

# Add colourbar
#fig.subplots_adjust(right=0.8)
#cbar_ax = fig.add_axes([0.85, 0.15, 0.05, 0.7])
fig.subplots_adjust(right=0.75)
cbar_ax = fig.add_axes([0.78, 0.15, 0.05, 0.7])
cbar = fig.colorbar(im, cax=cbar_ax)
cbar_add_central_tick(cbar)

cbar_ax.set_ylabel("$\\overline{{\\Delta^{{=}}_{{{0}|{1}}}}}: \\textrm{{Estimated homotypic survival difference}}$".format('A' if cell_type == 'B' else 'B', cell_type))

# Formatting
#fig.suptitle(r'$\eta_B = {}, \beta_B = {}$'.format(eta_A, beta_A), fontsize=16)
#fig.set_figheight(4.4)
#fig.set_figwidth(8)
fig.set_figheight(2.5)
fig.set_figwidth(3.544)

ax = axes[0]

eta_B_num_max_vertex = len(df_competition['eta_B'].unique())

eta_y = np.array([0.02, 0.05, 0.10, 0.15, 0.20, 0.24])
yticks = eta_B_num_max_vertex - eta_y / 0.02
yticklabels = [r'${:.2f}$'.format(elem) for elem in eta_y]

ax.set_yticks(yticks)
ax.set_yticklabels(yticklabels)

beta_x = np.array([0.1, 0.3, 0.5, 0.7, 0.9])
xticks = (beta_x - 0.10) / 0.10
xticklabels = [r'${:.2f}$'.format(elem) for elem in beta_x]

ax.set_xticks(xticks)
ax.set_xticklabels(xticklabels)

# Save figure
fig.savefig(image_png)
