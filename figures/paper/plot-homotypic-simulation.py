import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import pandas as pd

from formatting import set_mpl_customisation
from sys import argv

# Process arguments
if len(argv) != 3:
    print('Usage: {} DATA_CSV IMAGE_PNG'.format(argv[0]))
    exit(1)

data_csv = argv[1]
image_png = argv[2]

# Customise matplotlib
set_mpl_customisation('combination')

# Import dataset
df = pd.read_csv(data_csv)

# minxtract data
t = df['t'].to_numpy()
G1_cell_count = df['G1_cell_count'].to_numpy()
G2_cell_count = df['G2_cell_count'].to_numpy()
total_cell_count = G1_cell_count + G2_cell_count

# Create subplots
fig, ax = plt.subplots(1, figsize=(3.52, 1.6))

# Plot simulations
ax.step(t, G1_cell_count, 'b-', label=r'$\textrm{G1}$', where='post')
ax.step(t, G2_cell_count, 'r-', label=r'$\textrm{G2}$', where='post')
ax.step(t, total_cell_count, 'k-', label=r'$\textrm{Total}$', where='post')
#ax.set_title('$\\eta = {}, \\beta = {}$'.format(eta, beta))
ax.set_xlabel(r'$\textrm{Time}$')
ax.set_ylabel(r'$\textrm{Cell count}$')
ax.legend()

fig.tight_layout()
fig.savefig(image_png)
