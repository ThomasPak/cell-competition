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

N = 100
rho = np.linspace(0, 1, N)
rho_part_1 = np.linspace(0, 1/3, N)
rho_part_2 = np.linspace(1/3, 1, N)

plt.figure(figsize=(4, 4))

# Plot points
A = [0.4, 0.6]
Bp = [0.2, 0.3]
Bpp = [0.4, 0.45]
Cp = [0.4, 0.2]
Cpp = [0.6, 0.3]
D = [0.1, 0.16]
E = [0.25, 0.12]
E2 = [0.25, 0.08]
points = [A, Bp, Bpp, Cp, Cpp, D, E, E2]

for point in points:
    plt.plot(point[0], point[1], color='grey', marker='+')

# Plot region contours

plt.plot(rho, rho, label=r'$\rho$')
plt.plot(rho, (rho + 1)**2 / 4, label=r'$\frac{(1 + \rho)^2}{4}$')
plt.plot(rho, (rho - 1)**2 / 4, label=r'$\frac{(1 - \rho)^2}{4}$')
plt.plot(rho_part_2, (1- rho_part_2)**2, label=r'$(1 - \rho)^2$', color='red')
plt.plot(rho_part_1, (1- rho_part_1)**2, color='red', linestyle='dashed')

# Plot cross sections I and II
rhoI = 0.1
colorI = 'brown'

rhoII = 0.25
colorII = 'black'
plt.axvline(rhoI, linestyle='dotted', color=colorI, label=r'$\textrm{I}$')
plt.axvline(rhoII, linestyle='dotted', color=colorII, label=r'$\textrm{II}$')

# Name parameter regimes

plt.text(0.35, 0.7,     '$A$', fontsize=16)
plt.text(0.15, 0.25,    '$B\'$', fontsize=16)
plt.text(0.365, 0.42,    '$B\'\'$', fontsize=16)
plt.text(0.35, 0.2,     '$C\'$', fontsize=16)
plt.text(0.6, 0.2,      '$C\'\'$', fontsize=16)
plt.text(0.03, 0.12,    '$D$', fontsize=16)
plt.text(0.15, 0.05,    '$E$', fontsize=16)

plt.xlim([0, 1])
plt.ylim([0, 1])

plt.xticks([0, 0.25, 0.5, 0.75, 1])
plt.yticks([0, 0.25, 0.5, 0.75, 1])


plt.xlabel(r'$\rho$')
plt.ylabel(r'$\eta$')

plt.legend()

plt.tight_layout()
plt.savefig(image_png)
