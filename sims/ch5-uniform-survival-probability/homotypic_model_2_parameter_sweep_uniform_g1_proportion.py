import numpy as np
import sys

from heterotypic_model_2_parameter_sweep import heterotypic_model_2_parameter_sweep_general

# Homotypic uniform cell cycle model
# Test g1 proportion with update parameters for D and E

beta_As = np.array([0.3])
tG_As = np.array([100])
rho_As = np.array([0.3])
eta_As = np.array([0.1])
coef_As = np.array([1])

# Uniform ccm parameter sweep
tGs = np.array([100])
coefs = np.array([1])

beta_minimum_dict = {
        # A
        (0.4, 0.6) : [],
        # Bp
        (0.2, 0.3) : [1 - np.sqrt(0.3)],
        # Bpp
        (0.4, 0.45) : [],
        # Cp
        (0.4, 0.2) : [1 - np.sqrt(0.2)],
        # Cpp
        (0.6, 0.3) : [],
        # D
        (0.1, 0.16) : [],
        # E
        (0.25, 0.12) : [],
        # E low eta
        (0.25, 0.08) : [],
        }

num_betas = 10

rho_etas = [
        (0.4, 0.6),
        (0.2, 0.3),
        (0.4, 0.45),
        (0.4, 0.2),
        (0.6, 0.3),
        (0.1, 0.16),
        (0.25, 0.12),
        (0.25, 0.08),
        ]

tG_Bs = tGs
coef_Bs = coefs

n0s = [(0, 100)]
ccm = 'uniform'

num_iter = 100

tstart = 0
tend = np.inf

min_cell_count = 10
max_cell_count = 1000

datafile_prefix = 'homotypic-model-2-uniform-g1-proportion'

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
rho_A_data = []
eta_A_data = []
coef_A_data = []

beta_B_data = []
tG_B_data = []
rho_B_data = []
eta_B_data = []
coef_B_data = []

n0_A_data = []
n0_B_data = []

for rho_A in rho_As:
 for eta_A in eta_As:
  for beta_A in beta_As:
   for tG_A in tG_As:
    for coef_A in coef_As:
     for rho_B, eta_B in rho_etas:
       beta_minimums = beta_minimum_dict[(rho_B, eta_B)]
       betas = np.arange(rho_B, 1, (1 - rho_B) / num_betas)
       betas = np.append(betas, beta_minimums)
       for beta_B in betas:
        for tG_B in tG_Bs:
         for coef_B in coef_Bs:
          for n0 in n0s:

           n0_A = n0[0]
           n0_B = n0[1]
           for i in range(num_iter):

            if rho_A > beta_A  or rho_B > beta_B:
                   continue

            beta_A_data.append(beta_A)
            tG_A_data.append(tG_A)
            rho_A_data.append(rho_A)
            eta_A_data.append(eta_A)
            coef_A_data.append(coef_A)

            beta_B_data.append(beta_B)
            tG_B_data.append(tG_B)
            rho_B_data.append(rho_B)
            eta_B_data.append(eta_B)
            coef_B_data.append(coef_B)

            n0_A_data.append(n0_A)
            n0_B_data.append(n0_B)

# Perform parameter sweep
heterotypic_model_2_parameter_sweep_general(beta_A_data, tG_A_data,
    rho_A_data, eta_A_data, coef_A_data, beta_B_data, tG_B_data, rho_B_data,
    eta_B_data, coef_B_data, n0_A_data, n0_B_data, datafile_prefix, ccm,
    num_iter, tstart, tend, min_cell_count, max_cell_count, seed,
    sim_range)
