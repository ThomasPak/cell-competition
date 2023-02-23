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

homotypic_survival_probability_legend = True
homotypic_proliferation_regimes_legend = True
competitive_interactions_legend = False
heterotypic_proliferation_regimes_legend = True
heterotypic_competition_regimes_legend = True

orientation = '3-by-2'
#orientation = '2-by-3'

# Set limits
d_max = 2
eta_max = 2

# Homotypic survival probability
etas = [1, 0.5, 0.2, 0.1, 0.05]

# Competitive interactions
d_B = 1
eta_B = 1

# Heterotypic proliferation regimes
d_BI = 0.7
eta_BI = 1

d_BII = 1.7
eta_BII = 1

# Formatting
# Homotypic proliferation regimes
viability_colour = 'black'

extinct_colour = 'tab:grey'
explode_colour = 'tab:red'

# Competitive interactions
text_fs = 11.2

# Heterotypic proliferation regimes
coexistence_colour = 'black'
neutral_colour = 'black'
theta_A_homo_colour = 'tab:blue'
theta_AB_colour = 'tab:grey'
theta_BA_colour = 'tab:cyan'

AB_extinct_colour = 'tab:grey'
A_extinct_B_explode_colour = 'tab:green'
A_explode_B_extinct_colour = 'tab:orange'
AB_explode_colour = 'tab:red'

AB_extinct_hatch = ''
A_extinct_B_explode_hatch = '\\'
A_explode_B_extinct_hatch =  '/'
AB_explode_hatch = '/\\'

# Heterotypic competition regimes

# Create subplots
if orientation == '3-by-2':
    fig, axes = plt.subplots(3, 2, figsize=(8, 10))
    axes = [axes[0, 0], axes[0, 1], axes[1, 0], axes[1, 1], axes[2, 0], axes[2, 1]]
elif orientation == '2-by-3':
    fig, axes = plt.subplots(2, 3, figsize=(10, 6.7))
    axes = [axes[0, 0], axes[0, 1], axes[0, 2], axes[1, 0], axes[1, 1], axes[1, 2]]

## HOMOTYPIC SURVIVAL PROBABILITY ##
ax = axes[0]

theta_fun = lambda d, eta: 1 - np.exp(- eta / d)

N = 100
d_grid = np.linspace(1 / N, d_max, N - 1)

for eta in etas:

    theta = theta_fun(d_grid, eta)

    ax.plot(d_grid, theta, label=r'$\tilde{{\eta}} = {}$'.format(eta))

ax.set_title(r'$\textrm{Homotypic survival probability}$')
ax.set_xlabel(r'$d$')
ax.set_ylabel(r'$\lambda$')

if homotypic_survival_probability_legend:
    ax.legend()

## HOMOTYPIC PROLIFERATION REGIMES ##
def plot_bounds(ax, N=100, d_max=d_max, eta_max=eta_max):

    d = np.linspace(0, d_max, N)

    artists = []
    labels = []

    # Plot viability curve
    p, = ax.plot(d, np.log(2) * d, color=viability_colour)
    artists.append(p)
    labels.append(r'$\lambda = \frac{1}{2}$')

    # Fill in regions
    # Extinct
    ax.fill_between(d, 0, np.log(2) * d, alpha=0.5, facecolor=extinct_colour)

    # Explode
    ax.fill_between(d, np.log(2) * d, eta_max, alpha=0.5,
            facecolor=explode_colour)

    # Formatting
    ax.set_title(r'$\textrm{Homotypic proliferation regimes}$')
    ax.set_xlabel(r'$d$')
    ax.set_ylabel(r'$\tilde{\eta}$')

    ax.set_xlim([0, d_max])
    ax.set_ylim([0, eta_max])

    #xticks = [0, 1, d_max]
    #xticklabels = ['$0$', '$1$', '${}$'.format(d_max)]
    #ax.set_xticks(xticks)
    #ax.set_xticklabels(xticklabels)

    #yticks = [0, np.log(2), eta_max]
    #yticklabels = ['$0$', r'$\ln 2$', '${}$'.format(eta_max)]

    #ax.set_yticks(yticks)
    #ax.set_yticklabels(yticklabels)

    if homotypic_proliferation_regimes_legend:
        ax.legend(artists, labels, fontsize=10)

ax = axes[1]
plot_bounds(ax)

## HETEROTYPIC PROLIFERATION REGIMES ##
def plot_bounds(ax, eta_B, d_B, N=100, d_max=d_max, eta_max=eta_max):

    d_A = np.linspace(0, d_max, N)

    artists = []
    labels = []

    # Plot coexistence line
    p = ax.axhline(eta_B, color='k', linestyle='solid')
    artists.append(p)
    labels.append(r'$\Delta^{\neq}_{A|B} = 0$')

    # Plot neutral competition line
    p = ax.axvline(d_B, color=neutral_colour, linestyle='dashed')
    artists.append(p)
    labels.append(r"$\Delta^{=}_{A|B} = 0$")

    # Plot (d_B, eta_B)
    ax.plot(d_B, eta_B, 'o', color='tab:green')

    # Compute intersection of coexistence line with theta_A = 0.5
    d_intersect = eta_B / np.log(2)

    # Plot theta_A' = 0.5, i.e. ln 2 * d_A
    label1 = r"$\lambda_A = \frac{1}{2}$"
    d_A1 = np.linspace(0, d_intersect, N)
    d_A2 = np.linspace(d_intersect, d_max, N)
    p, = ax.plot(d_A1, np.log(2) * d_A1, ':',
            color=theta_A_homo_colour)

    p, = ax.plot(d_A2, np.log(2) * d_A2, '-',
            color=theta_A_homo_colour)
    artists.append(p)
    labels.append(label1)

    # Plot theta_B(A) = 0.5
    p = ax.vlines([d_intersect, d_intersect], [eta_B, 0],
            [eta_max, eta_B], colors=theta_BA_colour,
            linestyles=['solid', 'dotted'])
    artists.append(p)
    labels.append(r'$\xi_{B|A}^{\infty} = \frac{1}{2}$')

    # Plot theta_A(B) = 0.5
    if eta_B > np.log(2) * d_B:
        p = ax.axhline(np.log(2) * d_B, color=theta_AB_colour)
        artists.append(p)
        labels.append(r'$\xi_{A|B}^{\infty} = \frac{1}{2}$')

    # Fill in regions
    # A & B explode
    d_fill = np.linspace(0, d_intersect, N)
    ax.fill_between(d_fill, eta_B, eta_max,
            facecolor=AB_explode_colour, edgecolor='black', alpha=0.5,
            hatch=AB_explode_hatch)
    if eta_B > np.log(2) * d_B:
        ax.fill_between(d_A,
                np.log(2) * d_B,
                eta_B, facecolor=AB_explode_colour, alpha=0.5,
                edgecolor='black', hatch=AB_explode_hatch)

    # A goes extinct and B explodes
    if eta_B > np.log(2) * d_B:
        ax.fill_between(d_A,
                0, np.log(2) * d_B,
                facecolor=A_extinct_B_explode_colour, alpha=0.5,
                edgecolor='black', hatch=A_extinct_B_explode_hatch)

    # A explodes and B goes extinct
    d_fill = np.linspace(d_intersect, d_max, N)
    ax.fill_between(d_fill,
            np.log(2) * d_fill,
            eta_max, facecolor=A_explode_B_extinct_colour, alpha=0.5,
            edgecolor='black', hatch=A_explode_B_extinct_hatch)

    # A & B go extinct
    ax.fill_between(d_fill,
            eta_B,
            np.log(2) * d_fill,
            facecolor=AB_extinct_colour, alpha=0.5,
            edgecolor='black', hatch=AB_extinct_hatch)
    if eta_B < np.log(2) * d_B:
        ax.fill_between(d_A,
                0,
                eta_B,
                facecolor=AB_extinct_colour, alpha=0.5,
                edgecolor='black', hatch=AB_extinct_hatch)

    # Formatting
    ax.set_xlabel(r'$d_A$')
    ax.set_ylabel(r'$\tilde{\eta}_A$')

    ax.set_xlim([0, d_max])
    ax.set_ylim([0, eta_max])

    ax.set_xticks([0, d_B, d_max])
    ax.set_xticklabels(['$0$', r'$d_B$', '${}$'.format(d_max) ])

    ax.set_yticks([0, eta_B, eta_max])
    ax.set_yticklabels(['$0$', r'$\tilde{\eta}_B$', '${}$'.format(eta_max)])

    if heterotypic_proliferation_regimes_legend:
        ax.legend(artists, labels, fontsize=9)

ax = axes[2]
plot_bounds(ax, eta_BI, d_BI)
ax.set_title(r'$\textrm{Heterotypic proliferation regimes I}$')

ax = axes[3]
plot_bounds(ax, eta_BII, d_BII)
ax.set_title(r'$\textrm{Heterotypic proliferation regimes II}$')

## HOMOTYPIC COMPETITIVE INTERACTIONS ##
ax = axes[4]

# Plot region contours
ax.axhline(eta_B, color='k', label=r'$\Delta^{\neq}_{A|B} = 0$', linestyle='solid')
ax.axvline(d_B, color='k', label=r"$\Delta^{=}_{A|B} = 0$", linestyle='dashed')
ax.plot(d_B, eta_B, 'go', label=r'$\textrm{Neutral coexistence}$')

# Name parameter regimes
ax.text(1.08, 1.7,
        r'''$\textrm{A direct winner}$
        $\textrm{B direct loser}$''', fontsize=text_fs)
ax.text(0.08, 1.7,
#ax.text(0.03, 1.7,
        r'''$\textrm{A indirect winner}$
        $\textrm{B indirect loser}$''', fontsize=text_fs)

ax.text(1.08, 0.1,
#ax.text(1.04, 0.1,
        r'''$\textrm{A indirect loser}$
        $\textrm{B indirect winner}$''', fontsize=text_fs)
ax.text(0.08, 0.1,
        r'''$\textrm{A direct loser}$
        $\textrm{B direct winner}$''', fontsize=text_fs)

ax.text(0.1, 1.12,
        r"$\Delta^{\neq}_{A|B} > 0 $", fontsize=text_fs, color='k')
ax.text(0.1, 0.82,
        r"$\Delta^{\neq}_{A|B} < 0 $", fontsize=text_fs, color='k')

rot = 90
ax.text(0.83, 1.1,
        r"$\Delta^{=}_{A|B} < 0 $", fontsize=text_fs, rotation=rot,
        )
ax.text(1.05, 1.1,
        r"$\Delta^{=}_{A|B} > 0 $", fontsize=text_fs, rotation=rot,
        )

ax.set_xlim([0, d_max])
ax.set_ylim([0, eta_max])

ax.set_xlabel(r'$d_A$')
ax.set_ylabel(r'$\tilde{\eta}_A$')

ax.set_title(r'$\textrm{Competitive interactions}$')

if competitive_interactions_legend:
    ax.legend()

## HETEROTYPIC COMPETITION REGIMES ##
ax = axes[5]

from competition_regimes import plot_bounds_constant_emission

plot_bounds_constant_emission(ax, eta_BI, d_BI, d_max=d_max, eta_max=eta_max,
        show_legend=heterotypic_competition_regimes_legend)

ax.set_title(r'$\textrm{Competition regimes}$')

fig.tight_layout()

fig.savefig(image_png, bbox_inches='tight')
