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

coexistence_colour = 'black'
neutral_colour = 'black'
theta_B_homo_colour = 'tab:blue'
theta_BA_colour = 'tab:grey'
theta_AB_colour = 'tab:cyan'

AB_extinct_colour = 'tab:grey'
A_explode_B_extinct_colour = 'tab:green'
A_extinct_B_explode_colour = 'tab:orange'
AB_explode_colour = 'tab:red'

AB_extinct_hatch = ''
A_explode_B_extinct_hatch = ''
A_extinct_B_explode_hatch = ''
AB_explode_hatch = ''

#A_explode_B_extinct_hatch = '//'
#A_extinct_B_explode_hatch = '\\\\'
#AB_explode_hatch = '//\\\\'

def plot_bounds(ax, eta_A, beta_A, N=100, eta_max=0.25):

    beta_B = np.linspace(0, 1, N)

    artists = []
    labels = []

    # Plot coexistence line
    p, = ax.plot(beta_B, eta_A / beta_A * beta_B, color=coexistence_colour)
    artists.append(p)
    labels.append(r'$\Delta^{\neq}_{A|B} = 0$')

    # Plot neutral competition line
    p = ax.vlines(beta_A, 0, eta_max,
            colors=neutral_colour, linestyles='dashed')
    artists.append(p)
    labels.append(r"$\Delta^{=}_{A|B} = 0$")

    # Plot (beta_A, eta_A)
    ax.plot(beta_A, eta_A, 'o', color='tab:green')

    # Calculate intersection of coexistence line with theta_B = 0.5
    beta_intersect = 1 - eta_A / (beta_A * np.log(2))

    # Plot theta_B' = 0.5, i.e. ln 2 * beta_B * (1 - beta_B)
    label1 = r"$\lambda_A = \frac{1}{2}$"
    if beta_intersect > 0:
        beta_B1 = np.linspace(0, beta_intersect, int(N * beta_intersect))
        beta_B2 = np.linspace(beta_intersect, 1, int(N * (1 - beta_intersect)))
        p, = ax.plot(beta_B1, np.log(2) * beta_B1 * (1 - beta_B1), '-',
                color=theta_B_homo_colour)
        artists.append(p)
        labels.append(label1)

        label2 = None
    else:
        p = None
        beta_B2 = beta_B
        label2 = label1

    p, = ax.plot(beta_B2, np.log(2) * beta_B2 * (1 - beta_B2), ':',
            color=theta_B_homo_colour)
    if label2 != None:
        artists.append(p)
        labels.append(label2)

    # Plot theta_A(B) = 0.5
    if beta_intersect > 0:
        eta_intersect = eta_A / beta_A * (1 - eta_A / (beta_A * np.log(2)))
        p = ax.vlines([beta_intersect, beta_intersect], [eta_intersect, 0],
                [eta_max, eta_intersect], colors=theta_AB_colour,
                linestyles=['solid', 'dotted'])
        artists.append(p)
        labels.append(r'$\xi_{B|A}^{\infty} = \frac{1}{2}$')

    # Plot theta_B(A) = 0.5
    if eta_A > np.log(2) * beta_A * (1 - beta_A):
        p, = ax.plot(beta_B, np.log(2) * (1 - beta_A) * beta_B, theta_BA_colour)
        artists.append(p)
        labels.append(r'$\xi_{A|B}^{\infty} = \frac{1}{2}$')

    # Fill in regions
    # A & B explode
    beta_fill = np.linspace(max(beta_intersect, 0), 1, N)
    ax.fill_between(beta_fill, eta_A / beta_A * beta_fill, eta_max,
            facecolor=AB_explode_colour, edgecolor='black', alpha=0.5,
            hatch=AB_explode_hatch)
    if eta_A > np.log(2) * beta_A * (1 - beta_A):
        ax.fill_between(beta_B,
                np.log(2) * (1 - beta_A) * beta_B,
                eta_A / beta_A * beta_B, facecolor=AB_explode_colour, alpha=0.5,
                edgecolor='black', hatch=AB_explode_hatch)

    # A explodes and B goes extinct
    if eta_A > np.log(2) * beta_A * (1 - beta_A):
        ax.fill_between(beta_B,
                0, np.log(2) * (1 - beta_A) * beta_B,
                facecolor=A_explode_B_extinct_colour, alpha=0.5,
                edgecolor='black', hatch=A_explode_B_extinct_hatch)

    # A goes extinct and B explodes
    beta_fill = np.linspace(0, beta_intersect, N)
    if beta_intersect > 0:
        ax.fill_between(beta_fill,
                np.log(2) * beta_fill * (1 - beta_fill),
                eta_max, facecolor=A_extinct_B_explode_colour, alpha=0.5,
                edgecolor='black', hatch=A_extinct_B_explode_hatch)

    # A & B go extinct
    if beta_intersect > 0:
        ax.fill_between(beta_fill,
                eta_A / beta_A * beta_fill,
                np.log(2) * beta_fill * (1 - beta_fill),
                facecolor=AB_extinct_colour, alpha=0.5,
                edgecolor='black', hatch=AB_extinct_hatch)
    if eta_A < np.log(2) * beta_A * (1 - beta_A):
        ax.fill_between(beta_B,
                0,
                eta_A / beta_A * beta_B,
                facecolor=AB_extinct_colour, alpha=0.5,
                edgecolor='black', hatch=AB_extinct_hatch)

    # Formatting
    ax.set_xlabel(r'$\beta_A$')
    ax.set_ylabel(r'$\eta_A$')

    ax.set_xlim([0, 1])
    ax.set_ylim([0, eta_max])

    ax.set_xticks([0, beta_A, 1])
    ax.set_xticklabels(['$0$', r'$\beta_B$', '$1$' ])

    ax.set_yticks([0, eta_A, eta_max])
    ax.set_yticklabels(['$0$', r'$\eta_B$', '${}$'.format(eta_max)])

    ax.legend(artists, labels, fontsize=7)

# Create subplots
fig, axes = plt.subplots(1, 3, figsize=(7.48, 2.5))
ax1, ax2, ax3 = axes

eta_A1 = 0.2
beta_A1 = 0.2

plot_bounds(ax1, eta_A1, beta_A1)

eta_A2 = 0.2
beta_A2 = 0.8

plot_bounds(ax2, eta_A2, beta_A2)

eta_A3 = 0.1
beta_A3 = 0.4

plot_bounds(ax3, eta_A3, beta_A3)

ax1.set_title(r'$\textrm{{Cross Section I: }} \beta_B = {}, \eta_B = {}$'.format(beta_A1, eta_A1))
ax2.set_title(r'$\textrm{{Cross Section II: }} \beta_B = {}, \eta_B = {}$'.format(beta_A2, eta_A2))
ax3.set_title(r'$\textrm{{Cross Section III: }} \beta_B = {}, \eta_B = {}$'.format(beta_A3, eta_A3))

fig.tight_layout()
fig.savefig(image_png)
