import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt

from formatting import (set_mpl_customisation, cbar_add_central_tick)
from plot_survival_difference import (plot_heterotypic_survival_difference,
        plot_neutral_coexistence_point, plot_coexistence_curve)
from sys import argv

# Process arguments
if len(argv) != 3:
    print('Usage: {} DATA_CSV IMAGE_PNG'.format(argv[0]))
    exit(1)

data_csv = argv[1]
image_png = argv[2]

# Customise matplotlib
set_mpl_customisation('combination')

# Read datafile
df = pd.read_csv(data_csv)

# Vertex data preprocessing
df['num_divisions_A'] = df['final_cell_count_A'] - df['initial_cell_count_A'] + df['num_apoptosis_A'] + df['num_extrusion_A']
df['num_divisions_B'] = df['final_cell_count_B'] - df['initial_cell_count_B'] + df['num_apoptosis_B'] + df['num_extrusion_B']

df['n0_A'] = df['initial_cell_count_A']
df['n0_B'] = df['initial_cell_count_B']

df['num_deaths_A'] = df['num_apoptosis_A'] + df['num_extrusion_A']
df['num_deaths_B'] = df['num_apoptosis_B'] + df['num_extrusion_B']

# Compute survival differential
df['theta_A'] = df['num_divisions_A'] / (df['num_divisions_A'] + df['num_deaths_A'])
df['theta_B'] = df['num_divisions_B'] / (df['num_divisions_B'] + df['num_deaths_B'])

df['diff_theta'] = df['theta_B'] - df['theta_A']

# Check unique eta_A and beta_A values
unique_eta_As = df['eta_A'].unique()
unique_beta_As = df['beta_A'].unique()
assert len(unique_eta_As) == 1
assert len(unique_beta_As) == 1
eta_A = unique_eta_As[0]
beta_A = unique_beta_As[0]

# Assert that patterns are random and segregated
unique_patterns = sorted(list(df['pattern'].unique()))
assert unique_patterns[0] == 'random'
assert unique_patterns[1] == 'segregated'

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
    cur_df = df[df['pattern'] == pattern]

    # Set title
    ax.set_title("$\\textrm{{{}}}$".format(pattern[0].upper() + pattern[1:]))

    im, cur_max_abs_val = plot_heterotypic_survival_difference(ax, cur_df, beta_A, eta_A)

    ims.append(im)
    max_abs_val = max(max_abs_val, cur_max_abs_val)

    plot_coexistence_curve(ax, cur_df, beta_A, eta_A)

    plot_neutral_coexistence_point(ax, cur_df, beta_A, eta_A)

# Set clim
for im in ims:
    im.set_clim(-max_abs_val, max_abs_val)

# Add colourbar
fig.subplots_adjust(right=0.75)
cbar_ax = fig.add_axes([0.78, 0.15, 0.05, 0.7])
cbar = fig.colorbar(im, cax=cbar_ax)
cbar_add_central_tick(cbar)

cbar_ax.set_ylabel(r'$\overline{\Delta^{\neq}_{A|B}}: \textrm{Estimated heterotypic survival difference}$')

# Formatting
#fig.suptitle(r'$\eta_B = {}, \beta_B = {}$'.format(eta_A, beta_A), fontsize=7)
#fig.set_figheight(4.4)
#fig.set_figwidth(8)
fig.set_figheight(2.5)
fig.set_figwidth(3.544)

ax = axes[0]

eta_B_num_max_vertex = len(df['eta_B'].unique())

#eta_y = np.array([0.04, 0.10, 0.20, 0.30, 0.40, 0.48])
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
