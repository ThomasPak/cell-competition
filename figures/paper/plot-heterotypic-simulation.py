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
cell_count_A = df['cell_count_A'].to_numpy()
cell_count_B = df['cell_count_B'].to_numpy()

# Create subplots
fig, ax = plt.subplots(1, figsize=(3.52, 1.6))
#fig, ax = plt.subplots(1, figsize=(7.2, 1.6))

# Plot simulations
ax.step(t, cell_count_A, 'b-', label=r'$\textrm{A}$', where='post')
ax.step(t, cell_count_B, 'r-', label=r'$\textrm{B}$', where='post')
ax.set_xlabel(r'$\textrm{Time}$')
ax.set_ylabel(r'$\textrm{Cell count}$')
ax.legend()

fig.tight_layout()
fig.savefig(image_png)
