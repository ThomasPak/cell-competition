import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import pandas as pd

from formatting import set_mpl_customisation
from sys import argv

# Process arguments
if len(argv) != 3:
    print('Usage: {} DATA_CSV IMAGE_PNG'.format(argv[0]))
    exit(1)

data_csv = argv[1]
image_png = argv[2]

# Customise matplotlib
set_mpl_customisation()

# Import datasets
df = pd.read_csv(data_csv)

rhos = df['rho_B'].unique()
etas = df['eta_B'].unique()

df['exp_theta'] = df['num_divisions_B'] / (df['num_divisions_B'] + df['num_deaths_B'])

# Create subplots
fig, axes = plt.subplots(4, 2, figsize=(8, 10))
axes = axes.flatten()

## Create subplots
#fig = plt.figure(figsize=(8, 10))
#gs = fig.add_gridspec(4, 4)
#axes = []
#axes.append(fig.add_subplot(gs[0, 0:2]))
#axes.append(fig.add_subplot(gs[0, 2:]))
#axes.append(fig.add_subplot(gs[1, 0:2]))
#axes.append(fig.add_subplot(gs[1, 2:]))
#axes.append(fig.add_subplot(gs[2, 0:2]))
#axes.append(fig.add_subplot(gs[2, 2:]))
#axes.append(fig.add_subplot(gs[3, 1:-1]))

# Define points
A = [0.4, 0.6]
Bp = [0.2, 0.3]
Bpp = [0.4, 0.45]
Cp = [0.4, 0.2]
Cpp = [0.6, 0.3]
D = [0.1, 0.16]
E = [0.25, 0.12]
E2 = [0.25, 0.08]

points = [A, Bp, Bpp, Cp, Cpp, D, E, E2]
labels = ['A', 'B\'', 'B\'\'', 'C\'', 'C\'\'', 'D', 'E', 'E']

# Plot
theta_fun = lambda beta, eta, rho: np.maximum(0, np.minimum(1,
      (eta - (beta - rho) * (1 - beta)) / (2 * rho * (1 - beta))))

N = 200
beta_grid = np.linspace(1 / N, 1 - 1/N, N - 2)

for ax, point, label in zip(axes, points, labels):

    rho, eta = point

    # Plot theoretical
    theta = theta_fun(beta_grid, eta, rho)
    ax.plot(beta_grid, theta, label=r'$\textrm{Predicted}$')

    isLabelled = False

    # Plot valid beta span
    ax.axvspan(rho, 1, color='green', alpha=0.2)

    # Plot minimiser
    if label in ['B\'', 'B\'\'', 'C\'', 'C\'\'']:
        ax.axvline(1 - np.sqrt(eta), 0, 1, linestyle='dashed',
                   color='black', label=r'$\beta^*$')

    # Plot data
    betas = df[(df["rho_B"] == rho) & (df["eta_B"] == eta)]["beta_B"].unique()
    exp_theta_data = []
    positions = []
    for beta in betas:

        current_data = df[(df["eta_B"] == eta) & (df["rho_B"] == rho) &
            (df["beta_B"] == beta)]["exp_theta"].to_numpy()

        # 2: max cell count reached
        # 3: min cell count reached
        current_status = df[(df["eta_B"] == eta) & (df["rho_B"] == rho) &
            (df["beta_B"] == beta)]["status"].to_numpy()
        extinction_frequency = np.sum(current_status == 3) / current_status.size

        if isLabelled:
            thislabel = None
        else:
            isLabelled = True
            thislabel = r'$\textrm{Extinction frequency}$'
        #ax.plot(beta, extinction_frequency, 'o', color='black', label=thislabel)

        if np.min(current_data) != np.max(current_data):
            exp_theta_data.append(current_data)
            positions.append(beta)
        else:
            ax.plot(beta, current_data[0], 's', color='orange',
                    zorder=32)

    if len(exp_theta_data) > 0:
        ax.boxplot(exp_theta_data, positions=positions,
                widths=0.05,
                manage_ticks=False,
            zorder=32, whis=(0, 100))

    # Formatting
    ax.set_xticks(np.linspace(0, 1, 11))
    ax.set_yticks([0, 0.25, 0.5, 0.75, 1])

    ax.set_xlabel(r'$\beta$')
    ax.set_ylabel(r'$\hat{\lambda}: \textrm{Survival frequency}$')
    ax.set_ylim([-0.05, 1.05])

    ax.set_title('${}: \\rho = {}, \\eta = {}$'.format(label, rho, eta))
    ax.legend()

# Save figure with tight layout
fig.tight_layout()
plt.savefig(image_png)
