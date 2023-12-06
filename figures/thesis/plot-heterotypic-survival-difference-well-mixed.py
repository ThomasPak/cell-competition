import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt

from formatting import (set_mpl_customisation, cbar_add_central_tick)
from plot_survival_difference import (plot_heterotypic_survival_difference,
        plot_neutral_coexistence_point, plot_coexistence_curve,
        plot_coexisting_simulations)
from sys import argv

# Process arguments
if len(argv) != 3:
    print('Usage: {} DATA_CSV IMAGE_PNG'.format(argv[0]))
    exit(1)

data_csv = argv[1]
image_png = argv[2]

# Customise matplotlib
set_mpl_customisation()

# Read datafile
df = pd.read_csv(data_csv)

# Compute survival differential
df['theta_A'] = df['num_divisions_A'] / (df['num_divisions_A'] + df['num_deaths_A'])
df['theta_B'] = df['num_divisions_B'] / (df['num_divisions_B'] + df['num_deaths_B'])

df['diff_theta'] = df['theta_B'] - df['theta_A'] 

# Get unique eta_A and beta_A values
unique_eta_As = df['eta_A'].unique()
unique_beta_As = df['beta_A'].unique()

# Create subfigures
fig, axes = plt.subplots(len(unique_eta_As), len(unique_beta_As),
        sharex=True, sharey=True)

# Invert vertical ordering of axes
axes = axes[::-1, :]

# Initialise ims and max_abs_val
ims = []
max_abs_val = -np.inf

# Loop over eta_A and beta_A
for i, eta_A in enumerate(unique_eta_As):
    for j, beta_A in enumerate(unique_beta_As):

        # Select current axes
        ax = axes[i,j]

        # Select eta_A and beta_A
        cur_df = df[(df['eta_A'] == eta_A) & (df['beta_A'] == beta_A)]

        # Set title
        ax.set_title(r'$\eta_B = {}, \beta_B = {}$'.format(eta_A, beta_A))

        im, cur_max_abs_val = plot_heterotypic_survival_difference(ax, cur_df, beta_A, eta_A)

        ims.append(im)
        max_abs_val = max(max_abs_val, cur_max_abs_val)

        plot_coexistence_curve(ax, cur_df, beta_A, eta_A)

        plot_neutral_coexistence_point(ax, cur_df, beta_A, eta_A)

        #plot_coexisting_simulations(ax, cur_df, 0.2)

# Set clim
for im in ims:
    im.set_clim(-max_abs_val, max_abs_val)

# Add colourbar
fig.subplots_adjust(right=0.8)
cbar_ax = fig.add_axes([0.85, 0.15, 0.05, 0.7])
cbar = fig.colorbar(im, cax=cbar_ax)
cbar_add_central_tick(cbar)

cbar_ax.set_ylabel(r'$\overline{\Delta^{\neq}_{A|B}}: \textrm{Estimated heterotypic survival difference}$')

# Formatting
fig.set_figheight(10)
fig.set_figwidth(8)

ax = axes[0,0]

eta_B_num_max_wm = len(df['eta_B'].unique())

#eta_y = np.array([0.02, 0.10, 0.20, 0.30, 0.40, 0.50])
eta_y = np.array([0.01, 0.05, 0.10, 0.15, 0.20, 0.25])
yticks = eta_B_num_max_wm - eta_y / 0.01
yticklabels = [r'${:.2f}$'.format(elem) for elem in eta_y]

ax.set_yticks(yticks)
ax.set_yticklabels(yticklabels)

beta_x = np.array([0.05, 0.3, 0.5, 0.7, 0.95])
xticks = (beta_x - 0.05) / 0.05
xticklabels = [r'${:.2f}$'.format(elem) for elem in beta_x]

ax.set_xticks(xticks)
ax.set_xticklabels(xticklabels)

# Save figure
#fig.tight_layout()
fig.savefig(image_png, bbox_inches='tight')
