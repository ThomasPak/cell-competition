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

# Etas
etas = [1, 0.5, 0.2, 0.1, 0.05]

theta_fun = lambda beta, eta: 1 - np.exp(- eta / (beta * (1 - beta)))

# Plot
fig = plt.figure(figsize=(4, 4))

N = 100
beta_grid = np.linspace(1 / N, 1 - 1/N, N - 2)

for eta in etas:

    theta = theta_fun(beta_grid, eta)

    plt.plot(beta_grid, theta, label=r'$\eta = {}$'.format(eta))

plt.xlabel(r'$\beta$')
plt.ylabel(r'$\lambda$')

plt.xticks([0, 0.25, 0.5, 0.75, 1])
plt.yticks([0, 0.25, 0.5, 0.75, 1])

plt.legend()

plt.tight_layout()
plt.savefig(image_png)
