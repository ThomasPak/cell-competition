import numpy as np
from sys import argv

# One argument must be given
if len(argv) != 2:
    print("Usage: {} SIM_NUM".format(argv[0]))
    exit(2)

# This part of the input is the same for all simulations in this suite
base_string = """
output-directory = uniform-survival-probability

simulation-time = 100000
dt = 0.005
sampling-timestep-multiple = 200

min-cell-count = 10
max-cell-count = 1000
normalised-neighbours-g2-constants = 1

num-cells-across = 10
num-cells-up = 10

ancestors = 0

seed = {0}
cell-cycle-models = {1}
death-thresholds = {2}
g1-durations = {3}
g2-durations = {4}
uniform-g1-ranges = {5}
simulation-id = ccm_{1}_rho_{6:.2f}_eta_{7:.2f}_beta_{8:.2f}_sim_{9}
"""

# Define number of iterations
num_iter = 20

# Define values to iterate over for uniform cell cycle model
cell_cycle_model = 'uniform'
tG = 100
param_dict = {
        # A
        (0.4, 0.6) : [0.5],
        # Bp
        (0.2, 0.3) : [0.2, 1 - np.sqrt(0.3), 0.8],
        # Bpp
        (0.4, 0.45) : [0.4, 0.8],
        # Cp
        (0.4, 0.2) : [0.4, 1 - np.sqrt(0.2), 0.9],
        # Cpp
        (0.6, 0.3) : [0.6, 0.9],
        # D
        (0.1, 0.16) : [0.1, 0.5, 0.9],
        # E
        (0.25, 0.12) : [0.25, 0.5, 0.95],
        # E low eta
        (0.25, 0.08) : [0.25, 0.5, 0.95],
        }

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

# Read given simulation number and initialise simulation counter
param_strings = []
sim_counter = 0

# Generate parameter strings
for rho, eta in rho_etas:
    betas = param_dict[(rho, eta)]
    for beta in betas:
        for i in range(num_iter):

            seed = sim_counter + 1
            param_string = base_string.format(seed,
                    cell_cycle_model,
                    eta * tG,
                    beta * tG,
                    tG - beta * tG,
                    2 * rho * tG,
                    rho,
                    eta,
                    beta,
                    sim_counter)

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
