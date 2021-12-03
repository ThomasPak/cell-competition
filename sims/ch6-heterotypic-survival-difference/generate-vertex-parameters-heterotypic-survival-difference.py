import numpy as np
from sys import argv

# One argument must be given
if len(argv) != 2:
    print("Usage: {} SIM_NUM".format(argv[0]))
    exit(2)

## 15 June 2021
# Parameter sweep
beta_As = np.array([0.5])
tG_As = np.array([100])
eta_As = np.array([0.1])
coef_As = np.array([1])

beta_Bs = np.array([0.3, 0.5, 0.7])
tG_Bs = np.array([100])
eta_Bs = np.array([0.05, 0.1, 0.2])
coef_Bs = np.array([1])

# Spatial patterns for initial tissue
patterns = ['alternating', 'segregated', 'random', 'balanced', 'island',
        'islandinverted']

# Width of initial tissue
Nxs = np.array([20])

# Number of iterations
num_iter = 15

## End 15 June 2021

## 5 July 2021
## Note that I also changed the max cell count to 1000,
## where previously it was 5000
# Parameter sweep
beta_As = np.array([0.5])
tG_As = np.array([100])
eta_As = np.array([0.1])
coef_As = np.array([1])

beta_Bs = np.arange(0.1, 1, 0.1)
tG_Bs = np.array([100])
eta_Bs = np.arange(0.02, 0.25, 0.02)
coef_Bs = np.array([1])

# Spatial patterns for initial tissue
patterns = ['segregated', 'random']

# Width of initial tissue
Nxs = np.array([10])

# Number of iterations
num_iter = 20
## End 5 July 2021

# This part of the input is the same for all simulations in this suite
base_string = """
output-directory = heterotypic-survival-difference

simulation-time = 10000
dt = 0.005
sampling-timestep-multiple = 200

min-cell-count = 10
max-cell-count = 1000

min-cell-count-for-ancestor = {{0 : 0, 1 : 0}}

seed = {0}

num-cells-across = {1}
num-cells-up = {2}

ancestors = {3}

death-thresholds = {4}
g1-durations = {5}
g2-durations = {6}
normalised-neighbours-g2-constants = {7}

simulation-id = betaA_{8:.2f}_tGA_{9:.2f}_etaA_{10:.2f}_coefA_{11:.2f}_betaB_{12:.2f}_tGB_{13:.2f}_etaB_{14:.2f}_coefB_{15:.2f}_Nx_{1}_Ny_{2}_pattern_{16}_sim_{17}
"""

def generate_ancestors(pattern, Nx, Ny, seed=None):
    N = Nx * Ny
    if pattern == 'alternating':
        return ([0, 1] * N)[:N]
    elif pattern == 'segregated':
        return sorted(([0, 1] * N)[:N])
    elif pattern == 'random':
        rng = np.random.default_rng(seed)
        return rng.permutation(([0, 1] * N)[:N])
    elif pattern == 'balanced':
        return sum(((([0, 1, 1, 0] if j % 4 <= 1 else [1, 0, 0, 1]) * (Nx // 4 + 1))[0:Nx] for j in range(Ny)), [])
    elif pattern == 'island':
        assert Nx == Ny
        a = int(np.ceil(Nx / (4 + 2 * np.sqrt(2))))
        b = Nx - 2 * a
        return [0] * (a * Nx) + ([0] * a + [1] * b + [0] * a) * b + [0] * (a * Nx)
    elif pattern == 'islandinverted':
        assert Nx == Ny
        a = int(np.ceil(Nx / (4 + 2 * np.sqrt(2))))
        b = Nx - 2 * a
        return [1] * (a * Nx) + ([1] * a + [0] * b + [1] * a) * b + [1] * (a * Nx)
    else:
        raise Exception('ancestors pattern not recognised')

# Read given simulation number and initialise simulation counter
param_strings = []
sim_counter = 0

# Generate parameter strings
for beta_A in beta_As:
 for tG_A in tG_As:
  for eta_A in eta_As:
   for coef_A in coef_As:
    for beta_B in beta_Bs:
     for tG_B in tG_Bs:
      for eta_B in eta_Bs:
       for coef_B in coef_Bs:
        for Nx in Nxs:
         for pattern in patterns:
          for i in range(num_iter):

            seed = sim_counter + 1

            ancestors = generate_ancestors(pattern, Nx, Nx, seed)

            death_threshold_A = eta_A * coef_A * tG_A
            death_threshold_B = eta_B * coef_B * tG_B

            g1_duration_A = beta_A * tG_A
            g1_duration_B = beta_B * tG_B

            g2_duration_A = (1 - beta_A) * tG_A
            g2_duration_B = (1 - beta_B) * tG_B

            normalised_neighbours_g2_constant_A = coef_A
            normalised_neighbours_g2_constant_B = coef_B

            death_threshold_mapping = {0 : str(death_threshold_A), 1 : str(death_threshold_B)}
            g1_duration_mapping = {0 : str(g1_duration_A), 1 : str(g1_duration_B)}
            g2_duration_mapping = {0 : str(g2_duration_A), 1 : str(g2_duration_B)}
            normalised_neighbours_g2_constant_mapping = {0 :
                    str(normalised_neighbours_g2_constant_A), 1 :
                    str(normalised_neighbours_g2_constant_B)}

            param_string = base_string.format(
                    seed,
                    Nx,
                    Nx,
                    ','.join(map(str, ancestors)),
                    ','.join(map(death_threshold_mapping.get, ancestors)),
                    ','.join(map(g1_duration_mapping.get, ancestors)),
                    ','.join(map(g2_duration_mapping.get, ancestors)),
                    ','.join(map(normalised_neighbours_g2_constant_mapping.get, ancestors)),
                    beta_A,
                    tG_A,
                    eta_A,
                    coef_A,
                    beta_B,
                    tG_B,
                    eta_B,
                    coef_B,
                    pattern,
                    sim_counter,
                    )

            param_strings.append(param_string)
            sim_counter += 1

if argv[1] == 'all':

    for param_string in param_strings:
        print(param_string)

elif int(argv[1]) < len(param_strings):

    print(param_strings[int(argv[1])])

else:

    # If control arrives here, the given simulation number is invalid
    print("{} is an invalid simulation number".format(argv[1]))
    print("The largest valid simulation number is {}".format(sim_counter - 1))
    exit(2)
