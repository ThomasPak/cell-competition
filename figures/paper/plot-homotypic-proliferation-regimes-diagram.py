import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt

from formatting import set_mpl_customisation
from sys import argv

# Process arguments
if len(argv) != 2:
    print('Usage: {} IMAGE_PNG'.format(argv[0]))
    exit(1)

image_png = argv[1]

# Customise matplotlib
set_mpl_customisation('combination')

viability_colour = 'black'

extinct_colour = 'tab:grey'
explode_colour = 'tab:red'

def plot_bounds(ax, N=100, eta_max=0.25):

    beta = np.linspace(0, 1, N)

    artists = []
    labels = []

    # Plot viability curve
    p, = ax.plot(beta, np.log(2) * beta * (1 - beta), color=viability_colour)
    artists.append(p)
    labels.append(r'$\lambda = \frac{1}{2}$')
    
    # Fill in regions

    # Extinct
    ax.fill_between(beta, 0, np.log(2) * beta * (1 - beta), alpha=0.5,
            facecolor=extinct_colour)

    # Explode
    ax.fill_between(beta, np.log(2) * beta * (1 - beta), eta_max, alpha=0.5,
            facecolor=explode_colour)

    # Formatting
    ax.set_xlabel(r'$\beta$')
    ax.set_ylabel(r'$\eta$')

    ax.set_xlim([0, 1])
    ax.set_ylim([0, eta_max])

    xticks = [0, 0.5, 1]
    xticklabels = ['$0$', '$0.5$', '$1$']
    ax.set_xticks(xticks)
    ax.set_xticklabels(xticklabels)

    yticks = [0, np.log(2) / 4, eta_max]
    yticklabels = ['$0$', r'$\frac{\ln 2}{4}$', '${}$'.format(eta_max)]

    ax.set_yticks(yticks)
    ax.set_yticklabels(yticklabels)

    ax.legend(artists, labels, fontsize=7)

# Create subplots
#fig, ax = plt.subplots(1, 1, figsize=(2.38, 2.38))
fig, ax = plt.subplots(1, 1, figsize=(1.636, 1.73))

plot_bounds(ax)

fig.tight_layout()
fig.savefig(image_png)
