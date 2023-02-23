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

# Plot theoretical Eteff
N = 100
maxtdeath = 250

tG1 = 50
Eteff_fun = lambda tdeath: tG1 * (1 - np.exp(- tdeath / tG1))

tdeath = np.linspace(0, maxtdeath, N)
Eteff = Eteff_fun(tdeath)

ax = axes[0]

teff_data = []
positions = []
for eta in df_exp["eta"].unique():
    teff_data.append(df_exp[df_exp["eta"] == eta]["average_time_in_G1"].to_numpy())
    positions.append(int(eta * tG1))

ax.boxplot(teff_data, positions=positions, widths=10, manage_ticks=False,
        zorder=32, whis=(0, 100))

ax.plot(tdeath, Eteff, '-', label=r'$\textrm{Predicted}$')
#ax.axvline(tG1 * np.log(2),
#        color='black', linestyle='dashed',
#        label='Tipping point')

ax.set_xticks(np.arange(0, maxtdeath + tG1, tG1))
ax.set_xlabel('$t_\\dagger$')
ax.set_ylabel(r'$\hat{t}_{\textrm{eff}}: \textrm{Mean G1 duration}$')

ax.set_title(r'$\textrm{Exponential cell cycle model}$')

ax.legend()

## Uniform cell cycle model

# Plot theoretical Eteff
N = 100
maxtdeath = 110

tG1 = 50
def Eteff_fun(tdeath, r):
    
    output = 0 * tdeath + tdeath

    tdeath_s = tdeath[(tdeath >= tG1 - 0.5 * r) & (tdeath <= tG1 + 0.5 * r)]

    if r != 0:
        output[(tdeath >= tG1 - 0.5 * r) & (tdeath <= tG1 + 0.5 * r)] = \
                tG1 - 1/(2 * r) * (tdeath_s - (tG1 + 0.5 * r))**2

    output[tdeath > tG1 + 0.5 * r] = tG1

    return output

tdeath = np.linspace(0, maxtdeath, N)

for alpha, ax in zip(sorted(df_unif["alpha"].unique()), axes[1:]):

    r = 2 * alpha * tG1
    Eteff = Eteff_fun(tdeath, r)

    teff_data = []
    positions = []
    for eta in df_unif["eta"].unique():

        current_data = df_unif[(df_unif["eta"] == eta) & (df_unif["alpha"] ==
            alpha)][ "average_time_in_G1"].to_numpy()

        if np.min(current_data) != np.max(current_data):
            teff_data.append(current_data)
            positions.append(int(eta * tG1))
        else:
            ax.plot(int(eta * tG1), current_data[0], 's', color='orange',
                    zorder=32)

    if len(teff_data) > 0:
        ax.boxplot(teff_data, positions=positions, widths=7, manage_ticks=False,
                zorder=32, whis=(0, 100))

    ax.plot(tdeath, Eteff, '-', label=r'$\textrm{Predicted}$')

    #ax.axvline(tG1, color='black', linestyle='dashed',
    #        label='Tipping point')

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

    ax.set_xlabel('$t_\\dagger$')
    ax.set_ylabel(r'$\hat{t}_{\textrm{eff}}: \textrm{Mean G1 duration}$')
    ax.set_xlim([ax.get_xlim()[0], maxtdeath*1.05])

    ax.set_title('$\\textrm{{Uniform cell cycle model, }}r = {0:g}$'.format(r))

    ax.legend()

fig.tight_layout()
plt.savefig(image_png)
