import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt

coexistence_colour = 'black'
neutral_colour = 'black'
A_winner_viability_curve = 'tab:blue'
A_loser_viability_curve = 'tab:cyan'
B_loser_viability_curve = 'tab:cyan'

cell_competition_colour = 'tab:red'
incomplete_competition_colour = 'tab:orange'
indirection_competition_colour = 'tab:green'

A_winner_alpha = 0.5
B_winner_alpha = 0.5

weak_competition_hatch = 'x'

A_winner_strong_competition_hatch = weak_competition_hatch
A_winner_incomplete_competition_hatch = weak_competition_hatch
A_winner_indirect_competition_hatch = '\\'
B_winner_strong_competition_hatch = weak_competition_hatch
B_winner_incomplete_competition_hatch = weak_competition_hatch
B_winner_indirect_competition_hatch = '\\'


def plot_bounds(ax, eta_B, beta_B, N=100, eta_max=0.25):

    beta_A = np.linspace(0, 1, N)

    artists = []
    labels = []

    # Calculate intersection of coexistence line with theta_A = 0.5
    beta_intersect = 1 - eta_B / (beta_B * np.log(2))

    # Plot homotypic viability curve
    # last for zorder
    p, = ax.plot(beta_A, np.log(2) * beta_A * (1 - beta_A), '-',
            color=A_winner_viability_curve, zorder=32)
    artists.append(p)
    labels.append(r"$\lambda_A = \frac{1}{2}$")

    # Plot coexistence line (partially dotted)
    if beta_intersect > 0:
        beta_A1 = np.linspace(0, beta_intersect, int(N * beta_intersect))
        beta_A2 = np.linspace(beta_intersect, 1, int(N * (1 - beta_intersect)))

        p, = ax.plot(beta_A1, eta_B / beta_B * beta_A1, color=coexistence_colour,
                linestyle='dotted')
        p, = ax.plot(beta_A2, eta_B / beta_B * beta_A2, color=coexistence_colour,
                linestyle='solid', zorder=20)

    else:
        p, = ax.plot(beta_A, eta_B / beta_B * beta_A, color=coexistence_colour, zorder=20)

    artists.append(p)
    labels.append(r'$\Delta^{\neq}_{A|B} = 0$')

    # Plot neutral competition line (partially dotted)
    eta_neutral_intersect = np.log(2) * beta_B * (1 - beta_B)
    p = ax.vlines([beta_B, beta_B], [eta_neutral_intersect, 0],
            [eta_max, eta_neutral_intersect], colors=neutral_colour,
            linestyles=['dashed', 'dotted'])
    artists.append(p)
    labels.append(r"$\Delta^{=}_{L|W} = 0$")

    # Plot (beta_B, eta_B)
    ax.plot(beta_B, eta_B, 'o', color='tab:green', zorder=32)

    # Plot xi_{B|A} = 0.5
    if beta_intersect > 0:
        eta_intersect = eta_B / beta_B * (1 - eta_B / (beta_B * np.log(2)))
        p = ax.vlines([beta_intersect, beta_intersect], [eta_intersect, 0],
                [eta_max, eta_intersect], colors=B_loser_viability_curve,
                linestyles=['solid', 'dotted'])
        #artists.append(p)
        #labels.append(r'$\xi_{B|A}^{\infty} = \frac{1}{2}$')

    # Plot xi_{A|B} = 0.5
    if eta_B > np.log(2) * beta_B * (1 - beta_B):
        beta_A1 = np.linspace(0, beta_B, int(N * beta_B))
        beta_A2 = np.linspace(beta_B, 1, int(N * (1 - beta_B)))

        p, = ax.plot(beta_A1, np.log(2) * (1 - beta_B) * beta_A1,
                A_loser_viability_curve, linestyle='dotted')
        p, = ax.plot(beta_A2, np.log(2) * (1 - beta_B) * beta_A2,
                A_loser_viability_curve, linestyle='solid')

        artists.append(p)
        labels.append(r'$\xi_{L|W}^{\infty} = \frac{1}{2}$')

    # Fill in regimes for A winner
    # Cell competition regime
    beta_fill = np.linspace(0, beta_intersect, N)
    if beta_intersect > 0:
        ax.fill_between(beta_fill,
                np.log(2) * beta_fill * (1 - beta_fill),
                eta_max, facecolor=cell_competition_colour, alpha=A_winner_alpha,
                hatch=A_winner_strong_competition_hatch,
                edgecolor='black',
                linewidth=0
                )

    # Incomplete competition regime
    beta_fill = np.linspace(beta_intersect, beta_B, N)
    ax.fill_between(beta_fill,
            eta_B / beta_B * beta_fill,
            eta_max, facecolor=incomplete_competition_colour, alpha=A_winner_alpha,
            hatch=A_winner_incomplete_competition_hatch,
            edgecolor='black',
            linewidth=0
            )

    # Indirect competition regime
    beta_fill = np.linspace(beta_B, 1, N)
    ax.fill_between(beta_fill,
            eta_B / beta_B * beta_fill,
            eta_max, facecolor=indirection_competition_colour, alpha=A_winner_alpha,
            hatch=A_winner_indirect_competition_hatch,
            edgecolor='black',
            linewidth=0
            )

    # Fill in regimes for B winner
    # Cell competition regime
    beta_fill = np.linspace(beta_B, 1, N)
    ax.fill_between(beta_fill,
            np.log(2) * beta_fill * (1 - beta_fill),
            np.log(2) * (1 - beta_B) * beta_fill,
            facecolor=cell_competition_colour, alpha=B_winner_alpha,
            hatch=B_winner_strong_competition_hatch,
            edgecolor='black',
            linewidth=0
            )

    # Incomplete competition regime
    beta_fill = np.linspace(beta_B, 1, N)
    ax.fill_between(beta_fill,
            np.log(2) * (1 - beta_B) * beta_fill,
            eta_B / beta_B * beta_fill,
            facecolor=incomplete_competition_colour, alpha=B_winner_alpha,
            hatch=B_winner_incomplete_competition_hatch,
            edgecolor='black',
            linewidth=0
            )

    # Indirect competition regime
    beta_fill = np.linspace(beta_intersect, beta_B, N)
    ax.fill_between(beta_fill,
            np.log(2) * beta_fill * (1 - beta_fill),
            eta_B / beta_B * beta_fill,
            facecolor=indirection_competition_colour, alpha=B_winner_alpha,
            hatch=B_winner_indirect_competition_hatch,
            edgecolor='black',
            linewidth=0
            )

    # Formatting
    ax.set_xlabel(r'$\beta_A$')
    ax.set_ylabel(r'$\eta_A$')

    ax.set_xlim([0, 1])
    ax.set_ylim([0, eta_max])

    ax.set_xticks([0, beta_B, 1])
    ax.set_xticklabels(['$0$', r'$\beta_B$', '$1$' ])

    yticks = [0, eta_max]
    yticklabels = ['$0$', '${}$'.format(eta_max)]

    if eta_B < eta_max:
        yticks.append(eta_B)
        yticklabels.append(r'$\eta_B$')

    ax.set_yticks(yticks)
    ax.set_yticklabels(yticklabels)

    l = ax.legend(artists, labels, fontsize=10)
    l.set_zorder(100)

def plot_bounds_constant_emission(ax, eta_B, d_B, N=100, d_max=2, eta_max=2,
        show_legend=False):

    d_A = np.linspace(0, d_max, N)

    artists = []
    labels = []

    # Compute intersection of coexistence line with theta_A = 0.5
    d_intersect = eta_B / np.log(2)

    # Plot homotypic viability curve
    # last for zorder
    p, = ax.plot(d_A, np.log(2) * d_A, '-',
            color=A_winner_viability_curve, zorder=32)
    artists.append(p)
    labels.append(r"$\lambda_A = \frac{1}{2}$")

    # Plot coexistence line (partially dotted)
    d_A1 = np.linspace(0, d_intersect, N)
    d_A2 = np.linspace(d_intersect, d_max, N)

    p = ax.hlines([eta_B, eta_B], [0, d_intersect], [d_intersect, d_max],
            colors=coexistence_colour, linestyles=['solid', 'dotted'])

    artists.append(p)
    labels.append(r'$\Delta^{\neq}_{A|B} = 0$')

    # Plot neutral competition line (partially dotted)
    eta_neutral_intersect = np.log(2) * d_B
    p = ax.vlines([d_B, d_B], [eta_neutral_intersect, 0],
            [eta_max, eta_neutral_intersect], colors=neutral_colour,
            linestyles=['dashed', 'dotted'])
    artists.append(p)
    labels.append(r"$\Delta^{=}_{A|B} = 0$")

    # Plot (d_B, eta_B)
    ax.plot(d_B, eta_B, 'o', color='tab:green', zorder=32)

    # Plot xi_{B|A} = 0.5
    p = ax.vlines([d_intersect, d_intersect], [eta_B, 0],
            [eta_max, eta_B], colors=B_loser_viability_curve,
            linestyles=['solid', 'dotted'])
    #artists.append(p)
    #labels.append(r'$\xi_{B|A}^{\infty} = \frac{1}{2}$')

    # Plot xi_{A|B} = 0.5
    if eta_B > np.log(2) * d_B:
        d_A1 = np.linspace(0, d_B, N)
        d_A2 = np.linspace(d_B, 1, N)

        p = ax.hlines([d_B * np.log(2), d_B * np.log(2)],
                [0, d_B], [d_B, d_max], colors=A_loser_viability_curve,
                linestyles=['solid', 'dotted'])

        artists.append(p)
        labels.append(r'$\xi_{L|W}^{\infty} = \frac{1}{2}$')

    # Fill in regimes for A winner
    # Complete cell competition regime
    d_fill = np.linspace(d_intersect, d_max, N)
    ax.fill_between(d_fill,
            np.log(2) * d_fill,
            eta_max, facecolor=cell_competition_colour, alpha=A_winner_alpha,
            hatch=A_winner_strong_competition_hatch,
            edgecolor='black',
            linewidth=0
            )

    # Incomplete cell competition regime
    d_fill = np.linspace(d_B, d_intersect, N)
    ax.fill_between(d_fill,
            eta_B,
            eta_max, facecolor=incomplete_competition_colour, alpha=A_winner_alpha,
            hatch=A_winner_incomplete_competition_hatch,
            edgecolor='black',
            linewidth=0
            )

    # Indirect competition regime
    d_fill = np.linspace(0, d_B, N)
    ax.fill_between(d_fill,
            eta_B,
            eta_max, facecolor=indirection_competition_colour, alpha=A_winner_alpha,
            hatch=A_winner_indirect_competition_hatch,
            edgecolor='black',
            linewidth=0
            )

    # Fill in regimes for B winner
    # Complete cell competition regime
    d_fill = np.linspace(0, d_B, N)
    ax.fill_between(d_fill,
            np.log(2) * d_fill,
            np.log(2) * d_B,
            facecolor=cell_competition_colour, alpha=B_winner_alpha,
            hatch=B_winner_strong_competition_hatch,
            edgecolor='black',
            linewidth=0
            )

    # Incomplete competition regime
    d_fill = np.linspace(0, d_B, N)
    ax.fill_between(d_fill,
            np.log(2) * d_B,
            eta_B,
            facecolor=incomplete_competition_colour, alpha=B_winner_alpha,
            hatch=B_winner_incomplete_competition_hatch,
            edgecolor='black',
            linewidth=0
            )

    # Indirect competition regime
    d_fill = np.linspace(d_B, d_intersect, N)
    ax.fill_between(d_fill,
            np.log(2) * d_fill,
            eta_B,
            facecolor=indirection_competition_colour, alpha=B_winner_alpha,
            hatch=B_winner_indirect_competition_hatch,
            edgecolor='black',
            linewidth=0
            )

    # Formatting
    ax.set_xlabel(r'$d_A$')
    ax.set_ylabel(r'$\tilde{\eta}_A$')

    ax.set_xlim([0, d_max])
    ax.set_ylim([0, eta_max])

    ax.set_xticks([0, d_B, d_max])
    ax.set_xticklabels(['$0$', r'$d_B$', '${}$'.format(d_max) ])

    yticks = [0, eta_max]
    yticklabels = ['$0$', '${}$'.format(eta_max)]

    if eta_B < eta_max:
        yticks.append(eta_B)
        yticklabels.append(r'$\eta_B$')

    ax.set_yticks(yticks)
    ax.set_yticklabels(yticklabels)

    if show_legend:
        l = ax.legend(artists, labels, fontsize=10)
        l.set_zorder(100)

def plot_bounds_transformed(ax, eta_B, d_B, N=100, d_max=1, eta_max=1,
        show_legend=False):

    d_A = np.linspace(0, d_max, N)

    artists = []
    labels = []

    # Compute intersection of coexistence line with theta_A = 0.5
    d_intersect = eta_B

    # Plot homotypic viability curve
    # last for zorder
    p, = ax.plot(d_A, d_A, '-',
            color=A_winner_viability_curve, zorder=32)
    artists.append(p)
    labels.append(r"$\lambda_A = \frac{1}{2}$")

    # Plot coexistence line (partially dotted)
    d_A1 = np.linspace(0, d_intersect, N)
    d_A2 = np.linspace(d_intersect, d_max, N)

    p = ax.hlines([eta_B, eta_B], [0, d_intersect], [d_intersect, d_max],
            colors=coexistence_colour, linestyles=['solid', 'dotted'])

    artists.append(p)
    labels.append(r'$\Delta^{\neq}_{A|B} = 0$')

    # Plot neutral competition line (partially dotted)
    eta_neutral_intersect = d_B
    p = ax.vlines([d_B, d_B], [eta_neutral_intersect, 0],
            [eta_max, eta_neutral_intersect], colors=neutral_colour,
            linestyles=['dashed', 'dotted'])
    artists.append(p)
    labels.append(r"$\Delta^{=}_{A|B} = 0$")

    # Plot (d_B, eta_B)
    ax.plot(d_B, eta_B, 'o', color='tab:green', zorder=32)

    # Plot xi_{B|A} = 0.5
    p = ax.vlines([d_intersect, d_intersect], [eta_B, 0],
            [eta_max, eta_B], colors=B_loser_viability_curve,
            linestyles=['solid', 'dotted'])
    #artists.append(p)
    #labels.append(r'$\xi_{B|A}^{\infty} = \frac{1}{2}$')

    # Plot xi_{A|B} = 0.5
    if eta_B > d_B:
        d_A1 = np.linspace(0, d_B, N)
        d_A2 = np.linspace(d_B, 1, N)

        p = ax.hlines([d_B , d_B],
                [0, d_B], [d_B, d_max], colors=A_loser_viability_curve,
                linestyles=['solid', 'dotted'])

        artists.append(p)
        labels.append(r'$\xi_{L|W}^{\infty} = \frac{1}{2}$')

    # Fill in regimes for A winner
    # Complete cell competition regime
    d_fill = np.linspace(d_intersect, d_max, N)
    ax.fill_between(d_fill,
            d_fill,
            eta_max, facecolor=cell_competition_colour, alpha=A_winner_alpha,
            hatch=A_winner_strong_competition_hatch,
            edgecolor='black',
            linewidth=0
            )

    # Incomplete cell competition regime
    d_fill = np.linspace(d_B, d_intersect, N)
    ax.fill_between(d_fill,
            eta_B,
            eta_max, facecolor=incomplete_competition_colour, alpha=A_winner_alpha,
            hatch=A_winner_incomplete_competition_hatch,
            edgecolor='black',
            linewidth=0
            )

    # Indirect competition regime
    d_fill = np.linspace(0, d_B, N)
    ax.fill_between(d_fill,
            eta_B,
            eta_max, facecolor=indirection_competition_colour, alpha=A_winner_alpha,
            hatch=A_winner_indirect_competition_hatch,
            edgecolor='black',
            linewidth=0
            )

    # Fill in regimes for B winner
    # Complete cell competition regime
    d_fill = np.linspace(0, d_B, N)
    ax.fill_between(d_fill,
            d_fill,
            d_B,
            facecolor=cell_competition_colour, alpha=B_winner_alpha,
            hatch=B_winner_strong_competition_hatch,
            edgecolor='black',
            linewidth=0
            )

    # Incomplete competition regime
    d_fill = np.linspace(0, d_B, N)
    ax.fill_between(d_fill,
            d_B,
            eta_B,
            facecolor=incomplete_competition_colour, alpha=B_winner_alpha,
            hatch=B_winner_incomplete_competition_hatch,
            edgecolor='black',
            linewidth=0
            )

    # Indirect competition regime
    d_fill = np.linspace(d_B, d_intersect, N)
    ax.fill_between(d_fill,
            d_fill,
            eta_B,
            facecolor=indirection_competition_colour, alpha=B_winner_alpha,
            hatch=B_winner_indirect_competition_hatch,
            edgecolor='black',
            linewidth=0
            )

    # Formatting
    ax.set_xlabel(r'$d_A$')
    ax.set_ylabel(r'$\tilde{\eta}_A$')

    ax.set_xlim([0, d_max])
    ax.set_ylim([0, eta_max])

    ax.set_xticks([0, d_B, d_max])
    ax.set_xticklabels(['$0$', r'$d_B$', '${}$'.format(d_max) ])

    yticks = [0, eta_max]
    yticklabels = ['$0$', '${}$'.format(eta_max)]

    if eta_B < eta_max:
        yticks.append(eta_B)
        yticklabels.append(r'$\tilde{\eta}_B$')

    ax.set_yticks(yticks)
    ax.set_yticklabels(yticklabels)

    if show_legend:
        l = ax.legend(artists, labels, fontsize=10)
        l.set_zorder(100)
