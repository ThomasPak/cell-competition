import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import pandas as pd

from formatting import set_mpl_customisation
from sys import argv

# Process arguments
if len(argv) != 4:
    print('Usage: {} DATA_CSV APOPTOSIS_AT_CHECKPOINT[y/n] IMAGE_PNG'.format(argv[0]))
    exit(1)

data_csv = argv[1]
if argv[2] == 'y':
    apoptosis_at_checkpoint = True
elif argv[2] == 'n':
    apoptosis_at_checkpoint = False
else:
    print('Did not recognise APOPTOSIS_AT_CHECKPOINT[y/n] argument: {}'.format(argv[2]))
    exit(1)

image_png = argv[3]

# Customise matplotlib
set_mpl_customisation()

# Set ghosting
alpha_old = .3
alpha_new = 1

zorder_old = 32
zorder_new = 33

# Import datasets
df = pd.read_csv(data_csv)

#if apoptosis_at_checkpoint:
#    df = df[df['apoptosis_at_checkpoint'] == True]
#else:
#    df = df[df['apoptosis_at_checkpoint'] == False]

etas = df['eta'].unique()
betas = df['beta'].unique()

df['exp_theta'] = df['num_divisions'] / (df['num_divisions'] + df['num_deaths'])

# Theoretical
theta_fun = lambda beta, eta: 1 - np.exp(- eta / (beta * (1 - beta)))

# Create subplots
fig = plt.figure(figsize=(8, 10))
gs = fig.add_gridspec(3, 4)
axes = []
axes.append(fig.add_subplot(gs[0, 0:2]))
axes.append(fig.add_subplot(gs[0, 2:]))
axes.append(fig.add_subplot(gs[1, 0:2]))
axes.append(fig.add_subplot(gs[1, 2:]))
axes.append(fig.add_subplot(gs[2, 1:3]))

# Plot
N = 100
beta_grid = np.linspace(1 / N, 1 - 1/N, N - 2)

alpha_array = [alpha_new]
zorder_array = [zorder_new]

if apoptosis_at_checkpoint == True:
    alpha_array = [alpha_old, alpha_new]
    zorder_array = [zorder_old, zorder_new]

for eta, ax in zip(etas, axes):

    # Plot theoretical
    theta = theta_fun(beta_grid, eta)
    ax.plot(beta_grid, theta, label=r'$\textrm{Predicted}$')

    isLabelled = False

    # Plot data
    exp_theta_data = []
    positions = []
    for beta in betas:

        data_array = [ df[(df["eta"] == eta) &
            (df["beta"] == beta) &
            (df["apoptosis_at_checkpoint"] == False)
            ]["exp_theta"].to_numpy() ]

        if apoptosis_at_checkpoint == True:
            data_array.append(df[(df["eta"] == eta) &
                (df["beta"] == beta) &
                (df["apoptosis_at_checkpoint"] == True)
                ]["exp_theta"].to_numpy())

        for data, zorder, alpha in zip(data_array, zorder_array, alpha_array):

            if np.min(data) != np.max(data):

                boxplot_dict = ax.boxplot(data, positions=[beta],
                        widths=0.05,
                        manage_ticks=False,
                    zorder=zorder, whis=(0, 100))

                for key, value in boxplot_dict.items():
                    for item in value:
                        item.set_alpha(alpha)
            else:
                ax.plot(beta, data[0], 's', color='orange',
                        zorder=zorder, alpha=alpha)

    # Formatting
    ax.set_xlabel(r'$\beta$')
    ax.set_ylabel(r'$\hat{\lambda}: \textrm{Survival frequency}$')

    ax.set_xticks([0, 0.25, 0.5, 0.75, 1])
    ax.set_yticks([0, 0.25, 0.5, 0.75, 1])

    ax.set_ylim([0, 1.05])

    ax.set_title('$\\eta = {}$'.format(eta))
    ax.legend()

plt.tight_layout()
plt.savefig(image_png)
