import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import pandas as pd

from formatting import set_mpl_customisation
from sys import argv

# Process arguments
if len(argv) != 4:
    print('Usage: {} SIM1_DATA_CSV SIM2_DATA_CSV IMAGE_PNG'.format(argv[0]))
    exit(1)

sim1_data_csv = argv[1]
sim2_data_csv = argv[2]
image_png = argv[3]

# Customise matplotlib
set_mpl_customisation()

# Import dataset
max_df = pd.read_csv(sim1_data_csv)
min_df = pd.read_csv(sim2_data_csv)

# minxtract data
max_t = max_df['t'].to_numpy()
max_G1_cell_count = max_df['G1_cell_count'].to_numpy()
max_G2_cell_count = max_df['G2_cell_count'].to_numpy()
max_total_cell_count = max_G1_cell_count + max_G2_cell_count

min_t = min_df['t'].to_numpy()
min_G1_cell_count = min_df['G1_cell_count'].to_numpy()
min_G2_cell_count = min_df['G2_cell_count'].to_numpy()
min_total_cell_count = min_G1_cell_count + min_G2_cell_count

# Create subplots
fig, axes = plt.subplots(2, figsize=(8, 8))

ax1 = axes[0]
ax2 = axes[1]

# Plot simulations
ax1.step(max_t, max_G1_cell_count, 'b-', label=r'$\textrm{G1}$', where='post')
ax1.step(max_t, max_G2_cell_count, 'r-', label=r'$\textrm{G2}$', where='post')
ax1.step(max_t, max_total_cell_count, 'k-', label=r'$\textrm{Total}$', where='post')
ax1.set_title('$\\eta = 0.2, \\beta = 0.5$')
ax1.set_xlabel('$t$')
ax1.set_ylabel(r'$\textrm{Cell count}$')
ax1.legend()

ax2.step(min_t, min_G1_cell_count, 'b-', label=r'$\textrm{G1}$', where='post')
ax2.step(min_t, min_G2_cell_count, 'r-', label=r'$\textrm{G2}$', where='post')
ax2.step(min_t, min_total_cell_count, 'k-', label=r'$\textrm{Total}$', where='post')
ax2.set_title('$\\eta = 0.05, \\beta = 0.5$')
ax2.set_xlabel('$t$')
ax2.set_ylabel(r'$\textrm{Cell count}$')
ax2.legend()

fig.tight_layout()
fig.savefig(image_png)
