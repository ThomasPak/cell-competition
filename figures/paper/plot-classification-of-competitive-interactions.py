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
set_mpl_customisation('line')

text_fs = 7

N = 100
beta_B = np.linspace(0, 1, N)

beta_A = 0.6
eta_A = 0.09

#plt.figure(figsize=(6.4, 6.4/15*10))
plt.figure(figsize=(3.543, 2.8))

# Plot region contours

plt.plot(beta_B, eta_A / beta_A * beta_B, 'k', label=r'$\Delta^{\neq}_{A|B} = 0$')
plt.axvline(beta_A, color='k', label=r"$\Delta^{=}_{A|B} = 0$", linestyle='dashed')
plt.plot(beta_A, eta_A, 'go', label=r'$\textrm{Neutral coexistence}$')

# Name parameter regimes
plt.text(0.15, 0.16,
        r'''$\textrm{A direct winner}$
        $\textrm{B direct loser}$''', fontsize=text_fs)
plt.text(0.68, 0.16,
        r'''$\textrm{A indirect winner}$
        $\textrm{B indirect loser}$''', fontsize=text_fs)

rot=90
plt.text(0.54, 0.11,
        r"$\Delta^{=}_{A|B} > 0 $", fontsize=text_fs, color='k', rotation=rot)
plt.text(0.62, 0.11,
        r"$\Delta^{=}_{A|B} < 0 $", fontsize=text_fs, color='k', rotation=rot)

rot = 21
plt.text(0.32, 0.063,
        r"$\Delta^{\neq}_{A|B} > 0 $", fontsize=text_fs, rotation=rot,
        )
plt.text(0.355, 0.04,
        r"$\Delta^{\neq}_{A|B} < 0 $", fontsize=text_fs, rotation=rot,
        )

plt.text(0.28, 0.005,
        r'''$\textrm{A indirect loser}$
        $\textrm{B indirect winner}$''', fontsize=text_fs)
plt.text(0.65, 0.005,
        r'''$\textrm{A direct loser}$
        $\textrm{B direct winner}$''', fontsize=text_fs)

plt.xlim([0, 1])
plt.ylim([0, 0.25])

plt.xticks([0, beta_A, 1], labels=['$0$',
    '$\\beta_B$',
    '$1$'
    ])
plt.yticks([0, eta_A], labels=['$0$',
    '$\\eta_B$',
    ])

plt.xlabel(r'$\beta_A$')
plt.ylabel(r'$\eta_A$')

plt.legend(loc='upper left')

plt.tight_layout()
plt.savefig(image_png)
