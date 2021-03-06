import numpy as np
import sys

sys.path.append('../lib')
from g2_death_signal_parameter_sweep import heterotypic_model_2_parameter_sweep_general

# Exponential cell cycle model
# Fine sweep over beta_B and eta_B
# beta_A, eta_A = (0.2, 0.2), (0.8, 0.2), (0.4, 0.1)
beta_A_eta_As = [(0.2, 0.2), (0.8, 0.2), (0.4, 0.1)]
tG_As = np.array([100])
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

datafile_prefix = 'heterotypic-proliferation-regimes'

# Process arguments
if len(sys.argv) >= 2:
    seed = int(sys.argv[1])
else:
    seed = None

if len(sys.argv) == 4:
    sim_range = (int(sys.argv[2]), int(sys.argv[3]))
else:
    sim_range = None

# Generate parameters
beta_A_data = []
tG_A_data = []
eta_A_data = []
coef_A_data = []

beta_B_data = []
tG_B_data = []
eta_B_data = []
coef_B_data = []

n0_A_data = []
n0_B_data = []

for beta_A, eta_A in beta_A_eta_As:
 for tG_A in tG_As:
  for coef_A in coef_As:
   for beta_B in beta_Bs:
    for tG_B in tG_Bs:
     for eta_B in eta_Bs:
      for coef_B in coef_Bs:
       for n0 in n0s:

        n0_A = n0[0]
        n0_B = n0[1]
        for i in range(num_iter):

         beta_A_data.append(beta_A)
         tG_A_data.append(tG_A)
         eta_A_data.append(eta_A)
         coef_A_data.append(coef_A)

         beta_B_data.append(beta_B)
         tG_B_data.append(tG_B)
         eta_B_data.append(eta_B)
         coef_B_data.append(coef_B)

         n0_A_data.append(n0_A)
         n0_B_data.append(n0_B)

rho_A_data = [0] * len(beta_A_data)
rho_B_data = [0] * len(beta_B_data)

# Perform parameter sweep
heterotypic_model_2_parameter_sweep_general(beta_A_data, tG_A_data, rho_A_data,
        eta_A_data, coef_A_data, beta_B_data, tG_B_data, rho_B_data,
        eta_B_data, coef_B_data, n0_A_data, n0_B_data, datafile_prefix, ccm,
        num_iter, tstart, tend, min_cell_count, max_cell_count, seed,
        sim_range)
