import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import pandas as pd

from competition_regimes import plot_bounds_transformed

from formatting import set_mpl_customisation
from sys import argv

# Process arguments
if len(argv) != 2:
    print('Usage: {} IMAGE_PNG'.format(argv[0]))
    exit(1)

image_png = argv[1]

# Customise matplotlib
set_mpl_customisation('combination')

# Set limits
d_max = 1
eta_max = 1

eta_B = 0.6
d_B = 0.25

# Plot
fig, ax = plt.subplots(figsize=(3.544, 2.5))

plot_bounds_transformed(ax, eta_B, d_B, d_max=d_max, eta_max=eta_max,
        show_legend=True)

ax.set_title('$d_B = {}, \\tilde{{\eta}}_B = {}$'.format(d_B, eta_B))

fig.tight_layout()
fig.savefig(image_png)
