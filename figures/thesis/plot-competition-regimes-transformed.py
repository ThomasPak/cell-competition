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
set_mpl_customisation()

# Set limits
d_max = 1
eta_max = 1

eta_B = 0.6
d_B = 0.25

# Plot
fig, ax = plt.subplots(figsize=(4, 3))

plot_bounds_transformed(ax, eta_B, d_B, d_max=d_max, eta_max=eta_max,
        show_legend=True)

ax.set_title('$\\tilde{{\eta}}_B = {}, d_B = {}$'.format(eta_B, d_B))

fig.savefig(image_png, bbox_inches='tight')
