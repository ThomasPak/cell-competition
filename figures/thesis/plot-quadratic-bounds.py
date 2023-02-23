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

def plot_quadratic_bounds(ax, rho, N=100):

    l = lambda beta: (1 - beta) * (beta - rho)
    u = lambda beta: (1 - beta) * (beta + rho)

    beta_l = np.linspace(rho, 1, N)
    beta_u = np.linspace(-rho, 1, N)

    ax.plot(beta_u, u(beta_u), label='$u(\\beta)$', color=upper_colour)
    ax.plot(beta_l, l(beta_l), label='$l(\\beta)$', color=lower_colour)

    ax.axhline(0, color='black')
    ax.axvline(0, color='black')

    ax.set_xticks([-rho, 0, rho, 1])
    ax.set_xticklabels(['$-\\rho$', '$0$', '$\\rho$', '$1$'])

    ax.set_yticks([0, rho, (1 - rho)**2 / 4, (1 + rho)**2 / 4])
    ax.set_yticklabels(['$0$', '$\\rho$', '$\\frac{(1 - \\rho)^2}{4}$',
        '$\\frac{(1 + \\rho)^2}{4}$'])

    ax.set_xlabel('$\\beta$')
    ax.set_ylabel('$\\eta$')

    ax.set_title('$\\rho = {}$'.format(rho))

fig, axes = plt.subplots(1, 2, figsize=(8, 4))

#eta = .3
rho1 = 0.10
rho2 = 0.25

plot_quadratic_bounds(axes[0], rho1)
plot_quadratic_bounds(axes[1], rho2)

#etaA = (1 + .25)**2 / 4
eta1A = 0.38
eta2A = 0.46
eta1B = 0.28
eta2B = 0.35
eta1D = 0.17
eta2C = 0.20
eta1E = 0.05
eta2E = 0.10

axes[0].hlines(eta1A, 0, 1, label='$A$', color='red', zorder=32)
axes[1].hlines(eta2A, 0, 1, label='$A$', color='red', zorder=32)
axes[0].hlines(eta1B, 0, 1, label='$B$', color='green', zorder=32)
axes[1].hlines(eta2B, 0, 1, label='$B$', color='green', zorder=32)
axes[0].hlines(eta1D, 0, 1, label='$D$', color='gold', zorder=32)
axes[1].hlines(eta2C, 0, 1, label='$C$', color='purple', zorder=32)
axes[0].hlines(eta1E, 0, 1, label='$E$', color='C1', zorder=32)
axes[1].hlines(eta2E, 0, 1, label='$E$', color='C1', zorder=32)

axes[0].legend(loc='best').set_zorder(33)
axes[1].legend(loc='best').set_zorder(33)

axes[0].text((1+rho1)/2, 0.1, '$\\theta = 0$', fontsize=22,
        horizontalalignment='center')
axes[0].text((1-rho1)/2, 0.21, '$0 < \\theta < 1$', fontsize=22,
        horizontalalignment='center')
axes[0].text((1-rho1)/2, 0.335, '$\\theta = 1$', fontsize=22,
        horizontalalignment='center')

axes[1].text((1+rho2)/2, 0.05, '$\\theta = 0$', fontsize=22,
        horizontalalignment='center')
axes[1].text((1-rho2)/2, 0.27, '$0 < \\theta < 1$', fontsize=22,
        horizontalalignment='center')
axes[1].text((1-rho2)/2, 0.41, '$\\theta = 1$', fontsize=22,
        horizontalalignment='center')

fig.tight_layout()
fig.savefig(image_png)
