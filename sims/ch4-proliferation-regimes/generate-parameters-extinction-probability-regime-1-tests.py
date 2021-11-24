import numpy as np
from sys import argv

# One argument must be given
if not len(argv) is 2:
    print("Usage: {} SIM_NUM".format(argv[0]))
    exit(2)

# This part of the input is the same for all simulations in this suite
base_string = """
output-directory = extinction-probability-regime-1-tests

simulation-time = 100000
dt = 0.005
sampling-timestep-multiple = 200

num-cells-up = 1

ancestors = 0

g1-durations = 50
g2-durations = 50

base-rates = 1
normalised-neighbours-g2-constants = 0

birth-times = 0

seed = {0}
death-thresholds = {1}
num-cells-across = {2}
simulation-id = theta_{3:.2f}_n0_{2}_sim_{4}
"""

# Define values to iterate over
thetas = np.array([ 1/5, 1/3, 3/7 ])
death_thresholds = 50 * np.log(1 / (1 - thetas))
num_cells_across_array = [ 10, 30, 50 ]
num_iter = 100

# Read given simulation number and initialise simulation counter
param_strings = []
sim_counter = 0

# Generate parameter strings
for theta, death_threshold in zip(thetas, death_thresholds):
    for num_cells_across in num_cells_across_array:
        for i in range(num_iter):

            seed = sim_counter + 1
            param_string = base_string.format(seed, death_threshold,
                    num_cells_across, theta, sim_counter)

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
