import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import pandas as pd

from formatting import set_mpl_customisation
from sys import argv

# Process arguments
if len(argv) != 5:
    print('Usage: {} SIM1_DATA_CSV SIM2_DATA_CSV SIM3_DATA_CSV IMAGE_PNG'.format(argv[0]))
    exit(1)

sim1_data_csv = argv[1]
sim2_data_csv = argv[2]
sim3_data_csv = argv[3]
image_png = argv[4]

# Customise matplotlib
set_mpl_customisation()

# Import dataset
D_df = pd.read_csv(sim1_data_csv)
E_df = pd.read_csv(sim2_data_csv)
F_df = pd.read_csv(sim3_data_csv) # E low eta called F here for convenience

# Extract data
D_t = D_df['t'].to_numpy()
D_G1_cell_count = D_df['G1_cell_count'].to_numpy()
D_G2_cell_count = D_df['G2_cell_count'].to_numpy()
D_total_cell_count = D_G1_cell_count + D_G2_cell_count

E_t = E_df['t'].to_numpy()
E_G1_cell_count = E_df['G1_cell_count'].to_numpy()
E_G2_cell_count = E_df['G2_cell_count'].to_numpy()
E_total_cell_count = E_G1_cell_count + E_G2_cell_count

F_t = F_df['t'].to_numpy()
F_G1_cell_count = F_df['G1_cell_count'].to_numpy()
F_G2_cell_count = F_df['G2_cell_count'].to_numpy()
F_total_cell_count = F_G1_cell_count + F_G2_cell_count

# Create subplots
#fig, axes = plt.subplots(2, figsize=(8, 8))

#fig, axes = plt.subplots(3, figsize=(8, 10))
fig, axes = plt.subplots(3, figsize=(8, 8))

ax1 = axes[0]
ax2 = axes[1]
ax3 = axes[2]

# Plot simulations
ax1.step(D_t, D_G1_cell_count, 'b-', label=r'$\textrm{G1}$', where='post')
ax1.step(D_t, D_G2_cell_count, 'r-', label=r'$\textrm{G2}$', where='post')
ax1.step(D_t, D_total_cell_count, 'k-', label='Total', where='post')
ax1.set_title('$D: \\rho = 0.1, \\eta = 0.16, \\beta = 0.5$')
ax1.set_xlabel('$t$')
ax1.set_ylabel(r'$\textrm{Cell count}$')
ax1.legend()
ax1.set_xlim([0, D_t[-2]])
ax1.set_ylim([0, 1001])

start = 0
width = D_t[-2]

first_index = np.min(np.where(E_t > start))
last_index = np.min(np.where(E_t > start + width))
max_cell_count_in_window = np.max(E_total_cell_count[first_index:last_index])

ax2.step(E_t, E_G1_cell_count, 'b-', label=r'$\textrm{G1}$', where='post')
ax2.step(E_t, E_G2_cell_count, 'r-', label=r'$\textrm{G2}$', where='post')
ax2.step(E_t, E_total_cell_count, 'k-', label=r'$\textrm{Total}$', where='post')
ax2.set_title('$E: \\rho = 0.25, \\eta = 0.12, \\beta = 0.5$')
ax2.set_xlabel('$t$')
ax2.set_ylabel(r'$\textrm{Cell count}$')
ax2.legend()
ax2.set_xlim([start + 0, start + width])
ax2.set_ylim([ 0, max_cell_count_in_window + 1])

ax3.step(F_t, F_G1_cell_count, 'b-', label=r'$\textrm{G1}$', where='post')
ax3.step(F_t, F_G2_cell_count, 'r-', label=r'$\textrm{G2}$', where='post')
ax3.step(F_t, F_total_cell_count, 'k-', label=r'$\textrm{Total}$', where='post')
ax3.set_title('$E: \\rho = 0.25, \\eta = 0.08, \\beta = 0.5$')
ax3.set_xlabel('$t$')
ax3.set_ylabel(r'$\textrm{Cell count}$')
ax3.legend()
ax3.set_xlim([start + 0, start + width])
ax3.set_ylim([ 0, max_cell_count_in_window + 1])

ax1.set_xlim([0, 1000])
ax2.set_xlim([0, 1000])
ax3.set_xlim([0, 1000])

fig.tight_layout()
fig.savefig(image_png)
