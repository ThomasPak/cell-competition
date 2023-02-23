import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import pandas as pd

from formatting import set_mpl_customisation
from sys import argv

# Process arguments
if len(argv) != 4:
    print('Usage: {} EXPONENTIAL_DATA_CSV UNIFORM_DATA_CSV IMAGE_PNG'.format(argv[0]))
    exit(1)

exponential_data_csv = argv[1]
uniform_data_csv = argv[2]
image_png = argv[3]

# Customise matplotlib
set_mpl_customisation()

# Import datasets
df_exp = pd.read_csv(exponential_data_csv)
df_unif = pd.read_csv(uniform_data_csv)

# Create subplots
fig = plt.figure(figsize=(8, 10))
gs = fig.add_gridspec(3, 2)
axes = []
axes.append(fig.add_subplot(gs[2, 0:]))
axes.append(fig.add_subplot(gs[0, 0]))
axes.append(fig.add_subplot(gs[0, 1]))
axes.append(fig.add_subplot(gs[1, 0]))
axes.append(fig.add_subplot(gs[1, 1]))

## Exponential cell cycle model

tG1 = 50

ax = axes[0]

grouped = df_exp.groupby("eta")

tdeaths = df_exp["eta"].unique() * tG1
mean_sample_sizes = grouped["effective_g1_sample_size"].mean().to_numpy()

ax.bar(tdeaths, mean_sample_sizes, ec='black', width=15)
ax.axvline(tG1 * np.log(2),
        color='black', linestyle='dashed',
        label=r'$\textrm{Tipping point}$')

ax.set_xticks(np.linspace(0, 200, 5))

ax.set_yscale('log')

ax.set_xlabel('$t_\\dagger$')
ax.set_ylabel(r'$\textrm{Mean sample size}$')

ax.set_title(r'$\textrm{Exponential cell cycle model}$')

ax.legend()

## Uniform cell cycle model

tG1 = 50

for alpha, ax in zip(sorted(df_unif["alpha"].unique()), axes[1:]):

    r = 2 * alpha * tG1
    grouped = df_unif[df_unif["alpha"] == alpha].groupby("eta")

    tdeaths = df_unif["eta"].unique() * tG1
    mean_sample_sizes = grouped["effective_g1_sample_size"].mean().to_numpy()

    ax.bar(tdeaths, mean_sample_sizes, ec='black', width=7.5)
    ax.axvline(tG1, color='black', linestyle='dashed',
            label=r'$\textrm{Tipping point}$')

    if alpha != 0:

        ax.axvspan(tG1 - 0.5 * r, tG1 + 0.5 * r, color='green', alpha=0.2)

        ax.set_xticks([0, tG1 - 0.5 * r, tG1, tG1 + 0.5 * r, 100])
        ax.set_xticklabels(['$0$'] +
                ['${0:g}$'.format(tG1 - 0.5 * r)] +
                ['$50$'] +
                ['${0:g}$'.format(tG1 + 0.5 * r)]
                + ['$100$'])
    else:
        ax.set_xticks([0, tG1, 100])
        ax.set_xticklabels(['$0$'] +
                ['$50$'] +
                ['$100$'])

    ax.set_yscale('log')

    ax.set_xlabel('$t_\\dagger$')
    ax.set_ylabel(r'$\textrm{Mean sample size}$')

    ax.set_title('$\\textrm{{Uniform cell cycle model, }}r = {0:g}$'.format(r))

    ax.legend()

fig.tight_layout()
plt.savefig(image_png)
