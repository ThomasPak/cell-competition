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
set_mpl_customisation()

upper_colour = 'tab:grey'
lower_colour = 'tab:cyan'
viability_colour = 'black'
rho_colour = 'black'

lower_style = 'dashed'
upper_style = 'dashed'
rho_style = 'dashed'

extinct_colour = 'tab:grey'
explode_colour = 'tab:red'

extinct_hatch = ''
lower_extinct_hatch = '\\'
explode_hatch = ''
upper_explode_hatch = '/'

def plot_bounds(ax, rho, N=100, eta_max=0.5):

    beta = np.linspace(0, 1, N)

    artists = []
    labels = []

    # Plot rho curve (don't include in legend)
    if rho > 0:
        ax.vlines(rho, 0, eta_max, colors=rho_colour, linestyles=rho_style)

    # Plot viability curve
    p, = ax.plot(beta, beta * (1 - beta), color=viability_colour)
    artists.append(p)
    labels.append(r'$\lambda = \frac{1}{2}$')
    
    # Plot upper curve
    if rho > 0:
        p, = ax.plot(beta, (beta + rho) * (1 - beta), color=upper_colour, linestyle=upper_style)
        artists.append(p)
        labels.append(r'$u(\beta)$')

    # Plot lower curve
    if rho > 0:
        p, = ax.plot(beta, (beta - rho) * (1 - beta), color=lower_colour, linestyle=lower_style)
        artists.append(p)
        labels.append(r'$l(\beta)$')

    # Fill in regions

    # Lower extinct
    ax.fill_between(beta, 0, (beta - rho) * (1 - beta), alpha=0.5,
            facecolor=extinct_colour, hatch=lower_extinct_hatch)

    # Extinct
    if rho > 0:
        ax.fill_between(beta, (beta - rho) * (1 - beta), beta * (1 - beta), alpha=0.5,
                facecolor=extinct_colour, hatch=extinct_hatch)

    # Explode
    if rho > 0:
        ax.fill_between(beta, beta * (1 - beta), (beta + rho) * (1 - beta), alpha=0.5,
                facecolor=explode_colour, hatch=explode_hatch)

    # Upper explode
    ax.fill_between(beta, (beta + rho) * (1 - beta), eta_max, alpha=0.5,
            facecolor=explode_colour, hatch=upper_explode_hatch)

    # Formatting
    ax.set_xlabel(r'$\beta$')
    ax.set_ylabel(r'$\eta$')

    ax.set_xlim([0, 1])
    ax.set_ylim([0, eta_max])

    xticks = [0, 1]
    xticklabels = ['$0$', '$1$']
    if rho > 0:
        xticks.append(rho)
        xticklabels.append(r'$\rho$')
    ax.set_xticks(xticks)
    ax.set_xticklabels(xticklabels)

    yticks = [0, eta_max]
    yticklabels = ['$0$', '${}$'.format(eta_max)]

    if rho > 0:
        yticks.extend([rho, (1 - rho)**2 / 4, (1 + rho)**2 / 4])
        yticklabels.extend(['$\\rho$', '$\\frac{(1 - \\rho)^2}{4}$',
                '$\\frac{(1 + \\rho)^2}{4}$'])
    else:
        yticks.append(1/4)
        yticklabels.append(r'$\frac{1}{4}$')

    ax.set_yticks(yticks)
    ax.set_yticklabels(yticklabels)

    ax.legend(artists, labels, fontsize=10)

# Create subplots
fig, axes = plt.subplots(1, 2, figsize=(8, 4.4))

rhos = [
        0.10,
        0.25,
        ]

# Distinguish between B' and B'', C' and C''
# Not very interesting, so we dump this
#
# fig, axes = plt.subplots(4, 2, figsize=(8, 16))
# axes = np.hstack(axes)
# rhos = [0,
#         0.10,
#         3 - 2 * np.sqrt(2), 
#         0.25,
#         1/3,
#         0.35,
#         (3 - np.sqrt(5))/2,
#         0.5]

for ax, rho in zip(axes, rhos):

    plot_bounds(ax, rho)

    if np.isclose(rho, 3 - 2*np.sqrt(2)):
        rho_label = r'3 - 2 \sqrt{2}'
    elif np.isclose(rho, 1/3):
        rho_label = r'\frac{1}{3}'
    elif np.isclose(rho, (3 - np.sqrt(5))/2):
        rho_label = r'\frac{3 - \sqrt{5}}{2}'
    else:
        rho_label = rho

    if np.isclose(rho, 0.10):
        title = r'$\textrm{{Cross Section I: }} \rho = {}$'.format(rho)
    elif np.isclose(rho, 0.25):
        title = r'$\textrm{{Cross Section II: }} \rho = {}$'.format(rho)
    else:
        title = r'$\rho = {}$'.format(rho_label)

    ax.set_title(title)

fig.tight_layout()
fig.savefig(image_png)
