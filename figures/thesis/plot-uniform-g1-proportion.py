import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import pandas as pd

from formatting import set_mpl_customisation
from sys import argv

# Process arguments
if len(argv) != 2:
    print('Usage: {} IMAGE_PNG'.format(argv[0]))
    exit(1)

image_png = argv[1]

# Customise matplotlib
set_mpl_customisation()

# Create subplots
fig, axes = plt.subplots(4, 2, figsize=(8, 10))
axes = axes.flatten()

# Create subplots
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

    theta = theta_fun(beta_grid, eta, rho)

    ax.plot(beta_grid, theta)

    ax.axvspan(rho, 1, color='green', alpha=0.2)

    if label in ['B\'', 'B\'\'', 'C\'', 'C\'\'']:
        ax.axvline(1 - np.sqrt(eta), 0, 1, linestyle='dashed',
                   color='black', label=r'$\beta^*$')
        ax.legend()

    ax.set_xticks(np.linspace(0, 1, 11))
    ax.set_yticks([0, 0.25, 0.5, 0.75, 1])

    ax.set_xlabel(r'$\beta$')
    ax.set_ylabel(r'$\lambda$')
    ax.set_ylim([-0.05, 1.05])
    #ax.xlim([0, 1])

    ax.set_title('${}: \\rho = {}, \\eta = {}$'.format(label, rho, eta))


# Save figure with tight layout
fig.tight_layout()
plt.savefig(image_png)
