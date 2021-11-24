import pandas as pd
from itertools import product
from sys import argv

df = pd.read_csv('../analysis/mcc-found-parameters-04-11-2021.csv')

param_list = []

for index, row in df.iterrows():

    d = dict(row.drop(['simulation-id', 'seed']))
    l = ['{}={}'.format(key, val).lower() for key, val in d.items()]
    param_list.append(['initial-b-ratio=0.5'] + l)
    param_list.append(['initial-b-ratio=1.0'] + l)

    keys = ['initial-b-ratio'] + [ key for key, val in d.items() ]

seed_offset = 0
num_repeats = 20

if len(argv) == 2 and argv[1] == 'header':
    header_string = ','.join(keys)
    header_string += ',seed'
    header_string += ',simulation-id'
    print(header_string)
    exit()

sim_id = 0;

for parameter_list in param_list:

    for _ in range(num_repeats):

        parameter_string = ','.join(parameter_list)
        parameter_string += ',seed={:08d}'.format(sim_id + seed_offset)
        parameter_string += ',simulation-id={:06d}'.format(sim_id)

        print(parameter_string)

        sim_id += 1
