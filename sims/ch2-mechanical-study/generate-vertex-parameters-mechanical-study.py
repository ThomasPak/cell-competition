from itertools import product
from sys import argv

param_dict = {}

param_dict['initial-b-ratio'] = [0.5, 1.0];
param_dict['random-labelling'] = ['true', 'false'];
param_dict['cell-b-target-area'] = [1.0, 0.9, 1.1];
param_dict['cell-b-elasticity-parameter'] = [1.0, 0.9, 1.1];
param_dict['cell-b-contractility-parameter'] = [0.04, 0.03, 0.05];
param_dict['cell-b-line-tension-parameter'] = [0.12, 0.10, 0.14];
param_dict['cell-b-g1-duration'] = [30.0, 20.0, 40.0];
param_dict['cell-b-g2-duration'] = [70.0, 60.0, 80.0];
param_dict['cell-ab-line-tension-parameter'] = [0.12, 0.10, 0.14];

seed_offset = 0
num_repeats = 1

if len(argv) == 2 and argv[1] == 'header':
    header_string = ','.join(param_dict.keys())
    header_string += ',seed'
    header_string += ',simulation-id'
    print(header_string)
    exit()

def make_cartesian_product(param_dict):

    list_of_lists = []

    for name, values in param_dict.items():

        list_of_lists.append(['{}={}'.format(name, value) for value in values])

    return product(*list_of_lists)

sim_id = 0;

for parameter_list in make_cartesian_product(param_dict):

    for _ in range(num_repeats):

        parameter_string = ','.join(parameter_list)
        parameter_string += ',seed={:08d}'.format(sim_id + seed_offset)
        parameter_string += ',simulation-id={:06d}'.format(sim_id)

        print(parameter_string)

        sim_id += 1
