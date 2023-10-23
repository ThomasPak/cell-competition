from sys import argv
from scipy.stats import qmc

param_bounds = []

param_bounds.append(('cell-a-target-area', 0.5, 1.5))
param_bounds.append(('cell-b-target-area', 0.5, 1.5))
param_bounds.append(('cell-a-elasticity-parameter', 0.5, 1.5))
param_bounds.append(('cell-b-elasticity-parameter', 0.5, 1.5))
param_bounds.append(('cell-a-contractility-parameter', 0.01, 0.07))
param_bounds.append(('cell-b-contractility-parameter', 0.01, 0.07))
param_bounds.append(('cell-a-line-tension-parameter', 0.06, 0.18))
param_bounds.append(('cell-b-line-tension-parameter', 0.06, 0.18))
param_bounds.append(('cell-a-g1-duration', 0.0, 60.0))
param_bounds.append(('cell-b-g1-duration', 0.0, 60.0))
param_bounds.append(('cell-a-g2-duration', 40.0, 100.0))
param_bounds.append(('cell-b-g2-duration', 40.0, 100.0))
param_bounds.append(('cell-ab-line-tension-parameter', 0.06, 0.18))

seed_offset = 0
num_repeats = 1

if len(argv) == 2 and argv[1] == 'header':
    header_string = 'initial-b-ratio,random-labelling,'
    header_string += ','.join([item[0] for item in param_bounds])
    header_string += ',seed'
    header_string += ',simulation-id'
    print(header_string)
    exit()

sampler = qmc.LatinHypercube(d=13, strength=2, seed=1)
samples = sampler.random(2809)

parameters = []
sim_id = 0

for initial_b_ratio in [0.0, 0.5, 1.0]:

    preprefix = 'initial-b-ratio={},random-labelling=true'.format(initial_b_ratio)

    for sample in samples:

        prefix = preprefix

        for param_bound, value in zip(param_bounds, sample):

            name, low, high = param_bound

            prefix += ',{}={}'.format(name, low + value * (high - low))

        for _ in range(num_repeats):

            parameter_string = prefix
            parameter_string += ',seed={:08d}'.format(sim_id + seed_offset)
            parameter_string += ',simulation-id={:06d}'.format(sim_id)

            print(parameter_string)

            sim_id += 1
