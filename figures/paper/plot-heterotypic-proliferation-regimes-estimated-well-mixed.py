import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt

from formatting import (set_mpl_customisation, cbar_add_central_tick)
from plot_survival_difference import (plot_survival_frequency,
        plot_neutral_coexistence_point, plot_coexistence_curve,
        plot_neutral_competition_curve, plot_B_loser_viability_curve,
        plot_A_loser_viability_curve, plot_B_winner_viability_curve)
from sys import argv

# Process arguments
if len(argv) != 5:
    print('Usage: {} DATA_I_CSV DATA_II_CSV DATA_III_CSV IMAGE_PNG'.format(argv[0]))
    exit(1)

data_I_csv = argv[1]
data_II_csv = argv[2]
data_III_csv = argv[3]
image_png = argv[4]

# Customise matplotlib
set_mpl_customisation('combination')

# Read datafile
dfI = pd.read_csv(data_I_csv)
dfII = pd.read_csv(data_II_csv)
dfIII = pd.read_csv(data_III_csv)

dfs = [dfI, dfII, dfIII]

# Compute survival frequencies
for df in dfs:
    df['theta_A'] = df['num_divisions_A'] / (df['num_divisions_A'] + df['num_deaths_A'])
    df['theta_B'] = df['num_divisions_B'] / (df['num_divisions_B'] + df['num_deaths_B'])

# Get unique eta_A and beta_A values
eta_As = []
beta_As = []
for df in dfs:
    eta_A = df['eta_A'].unique()
    beta_A = df['beta_A'].unique()

    assert len(eta_A) == len(beta_A) == 1

    eta_As.append(eta_A[0])
    beta_As.append(beta_A[0])

# Create subfigures
#fig, axes = plt.subplots(len(dfs), 2,
#        sharex=True, sharey=True)

fig, axes = plt.subplots(2, len(dfs),
        sharex=True, sharey=True)

axes_new = axes.transpose()

axes = axes_new

# Loop over cross sections
for df, axs, eta_A, beta_A, cross_section in zip(dfs, axes, eta_As, beta_As,
        ['I', 'II', 'III']):

    axA, axB = axs

    # Set title
    axA.set_title('$\\textrm{{Cross Section {}: A}}$'.format(cross_section))
    axB.set_title('$\\textrm{{Cross Section {}: B}}$'.format(cross_section))

    # This might be confusing because I switched the labels around
    im = plot_survival_frequency(axA, df, beta_A, eta_A, 'B')
    im = plot_survival_frequency(axB, df, beta_A, eta_A, 'A')

    for ax in axs:
        plot_neutral_competition_curve(ax, df, beta_A, eta_A)
        plot_coexistence_curve(ax, df, beta_A, eta_A)
        plot_neutral_coexistence_point(ax, df, beta_A, eta_A)
        plot_B_winner_viability_curve(ax, df, beta_A, eta_A)

        if cross_section in ['I', 'II']:
            plot_B_loser_viability_curve(ax, df, beta_A, eta_A)

        if cross_section in ['II', 'III']:
            plot_A_loser_viability_curve(ax, df, beta_A, eta_A)

# Add colourbar
#fig.subplots_adjust(right=0.8)
fig.subplots_adjust(right=0.8, top=0.92)
cbar_ax = fig.add_axes([0.85, 0.15, 0.05, 0.7])
cbar = fig.colorbar(im, cax=cbar_ax)
cbar.set_ticks(np.arange(0, 1.1, 0.1))

cbar_ax.set_ylabel(r'$\overline{\xi}: \textrm{Estimated survival frequency}$')

# Formatting
fig.set_figheight(5)
fig.set_figwidth(7.48)

ax = axes[-1,0]

eta_B_num_max_wm = len(dfI['eta_B'].unique())

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
fig.savefig(image_png)
