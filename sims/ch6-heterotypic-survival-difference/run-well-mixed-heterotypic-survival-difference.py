import numpy as np
import sys

sys.path.append('../lib')
from g2_death_signal_parameter_sweep import heterotypic_model_2_parameter_sweep

# Exponential cell cycle model
# Fine sweep over beta_B and eta_B
# Coarse sweep over beta_A and eta_A

beta_As = np.array([0.3, 0.5, 0.7])
tG_As = np.array([100])
eta_As = np.array([0.05, 0.1, 0.2])
coef_As = np.array([1])

beta_Bs = np.arange(0.05, 1, 0.05)
tG_Bs = np.array([100])
eta_Bs = np.arange(0.01, 0.26, 0.01)
coef_Bs = np.array([1])

n0s = [(50, 50)]
ccm = 'exponential'

num_iter = 50

tstart = 0
tend = 10000

min_cell_count = 10
max_cell_count = 1000

datafile_prefix = 'heterotypic-survival-difference'

# 23 July 2021
#
# I have gotten interested in terminating the simulation as soon as any cell
# type goes extinct, since the simulation reduces to a homotypic simulation at
# that point, diluting the data from the heterotypic simulation.

min_cell_count_for_clone = {0 : 0, 1 : 0}

# End 23 July 2021

# Process arguments
if len(sys.argv) >= 2:
    seed = int(sys.argv[1])
else:
    seed = None

if len(sys.argv) == 4:
    sim_range = (int(sys.argv[2]), int(sys.argv[3]))
else:
    sim_range = None

# Perform parameter sweep
heterotypic_model_2_parameter_sweep(beta_As, tG_As, eta_As, coef_As,
    beta_Bs, tG_Bs, eta_Bs, coef_Bs, n0s, datafile_prefix, ccm, num_iter,
    tstart, tend, min_cell_count, max_cell_count, seed, sim_range,
    min_cell_count_for_clone=min_cell_count_for_clone)
