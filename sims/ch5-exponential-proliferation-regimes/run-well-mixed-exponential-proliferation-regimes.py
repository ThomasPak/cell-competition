import numpy as np
import sys

sys.path.append('../lib')
from g2_death_signal_parameter_sweep import heterotypic_model_2_parameter_sweep

# Exponential cell cycle model
# Control simulations
# Fine sweep over beta_B and eta_B
# Dummy sweep over beta_A and eta_A

beta_As = np.array([0.5])
tG_As = np.array([100])
eta_As = np.array([0.1])
coef_As = np.array([1])

beta_Bs = np.arange(0.05, 1, 0.05)
tG_Bs = np.array([100])
eta_Bs = np.arange(0.01, 0.26, 0.01)
coef_Bs = np.array([1])

n0s = [(0, 100)]
ccm = 'exponential'

num_iter = 50

tstart = 0
tend = 10000

min_cell_count = 10
max_cell_count = 1000

datafile_prefix = 'heterotypic-model-2-exponential-control'

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
    tstart, tend, min_cell_count, max_cell_count, seed, sim_range)
