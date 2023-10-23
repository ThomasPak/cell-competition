import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from sys import argv

# Process arguments
if len(argv) != 5:
    print('Usage: {} DATA_HOMO_A_CSV DATA_HOMO_B_CSV DATA_HETERO_CSV IMAGE_PNG'.format(argv[0]))
    exit(1)

data_homo_A_csv = argv[1]
data_homo_B_csv = argv[2]
data_hetero_csv = argv[3]
image_png = argv[4]

# Read data
dfhomA = pd.read_csv(data_homo_A_csv)
dfhomB = pd.read_csv(data_homo_B_csv)
dfhet = pd.read_csv(data_hetero_csv)

## Preprocess data
df = pd.DataFrame()

# survival frequencies
df['homA_initial_cell_count'] = dfhomA['initial_cell_count']
df['homB_initial_cell_count'] = dfhomB['initial_cell_count']
df['hetA_initial_cell_count'] = dfhet['initial_cell_count_0']
df['hetB_initial_cell_count'] = dfhet['initial_cell_count_1']

df['homA_final_cell_count'] = dfhomA['final_cell_count']
df['homB_final_cell_count'] = dfhomB['final_cell_count']
df['hetA_final_cell_count'] = dfhet['final_cell_count_0']
df['hetB_final_cell_count'] = dfhet['final_cell_count_1']

df['homA_num_apoptosis'] = dfhomA['num_apoptosis']
df['homB_num_apoptosis'] = dfhomB['num_apoptosis']
df['hetA_num_apoptosis'] = dfhet['num_apoptosis_0']
df['hetB_num_apoptosis'] = dfhet['num_apoptosis_1']

df['homA_num_extrusion'] = dfhomA['num_extrusion']
df['homB_num_extrusion'] = dfhomB['num_extrusion']
df['hetA_num_extrusion'] = dfhet['num_extrusion_0']
df['hetB_num_extrusion'] = dfhet['num_extrusion_1']

df['homA_num_deaths'] = df['homA_num_extrusion'] + df['homA_num_apoptosis']
df['homB_num_deaths'] = df['homB_num_extrusion'] + df['homB_num_apoptosis']
df['hetA_num_deaths'] = df['hetA_num_extrusion'] + df['hetA_num_apoptosis']
df['hetB_num_deaths'] = df['hetB_num_extrusion'] + df['hetB_num_apoptosis']

df['homA_num_divisions'] = df['homA_final_cell_count'] - df['homA_initial_cell_count'] + df['homA_num_deaths']
df['homB_num_divisions'] = df['homB_final_cell_count'] - df['homB_initial_cell_count'] + df['homB_num_deaths']
df['hetA_num_divisions'] = df['hetA_final_cell_count'] - df['hetA_initial_cell_count'] + df['hetA_num_deaths']
df['hetB_num_divisions'] = df['hetB_final_cell_count'] - df['hetB_initial_cell_count'] + df['hetB_num_deaths']

df['homA_theta'] = df['homA_num_divisions'] / (df['homA_num_divisions'] + df['homA_num_deaths'])
df['homB_theta'] = df['homB_num_divisions'] / (df['homB_num_divisions'] + df['homB_num_deaths'])
df['hetA_theta'] = df['hetA_num_divisions'] / (df['hetA_num_divisions'] + df['hetA_num_deaths'])
df['hetB_theta'] = df['hetB_num_divisions'] / (df['hetB_num_divisions'] + df['hetB_num_deaths'])

# simulation ids
df['homA_simulation_id'] = dfhomA['simulation-id']
df['homB_simulation_id'] = dfhomB['simulation-id']
df['het_simulation_id'] = dfhet['simulation-id']

# Drop failed simulations
df['homA_State'] = dfhomA['State']
df['homB_State'] = dfhomB['State']
df['het_State'] = dfhet['State']

failed = df[(df['homA_State'] == 'FAILED') |
        (df['homB_State'] == 'FAILED') |
        (df['het_State'] == 'FAILED') |
        (df['homA_State'] == 'TIMEOUT') |
        (df['homB_State'] == 'TIMEOUT') |
        (df['het_State'] == 'TIMEOUT')
        ]

df = df.drop(failed.index)

num_failed = len(failed)

print('**************************************')
print('Dropped {} failed/timeout simulation triplets'.format(num_failed))

# Count cases
counts = []

# Row 1
outcome = df[(df['hetB_theta'] < 0.5) &
                        (df['hetA_theta'] < 0.5) &
                        (df['homB_theta'] < 0.5) &
                        (df['homA_theta'] < 0.5)]
counts.append(len(outcome))

outcome = df[   (df['hetB_theta'] < 0.5) &
                        (df['hetA_theta'] < 0.5) &
                        (df['homB_theta'] < 0.5) &
                        (df['homA_theta'] >= 0.5)]
counts.append(len(outcome))

outcome = df[   (df['hetB_theta'] < 0.5) &
                        (df['hetA_theta'] < 0.5) &
                        (df['homB_theta'] >= 0.5) &
                        (df['homA_theta'] < 0.5)]
counts.append(len(outcome))

outcome = df[   (df['hetB_theta'] < 0.5) &
                        (df['hetA_theta'] < 0.5) &
                        (df['homB_theta'] >= 0.5) &
                        (df['homA_theta'] >= 0.5)]
counts.append(len(outcome))

# Row 2

outcome = df[   (df['hetB_theta'] < 0.5) &
                        (df['hetA_theta'] >= 0.5) &
                        (df['homB_theta'] < 0.5) &
                        (df['homA_theta'] < 0.5)]
counts.append(len(outcome))

outcome = df[   (df['hetB_theta'] < 0.5) &
                        (df['hetA_theta'] >= 0.5) &
                        (df['homB_theta'] < 0.5) &
                        (df['homA_theta'] >= 0.5)]
counts.append(len(outcome))

outcome = df[   (df['hetB_theta'] < 0.5) &
                        (df['hetA_theta'] >= 0.5) &
                        (df['homB_theta'] >= 0.5) &
                        (df['homA_theta'] < 0.5)]
counts.append(len(outcome))

outcome = df[   (df['hetB_theta'] < 0.5) &
                        (df['hetA_theta'] >= 0.5) &
                        (df['homB_theta'] >= 0.5) &
                        (df['homA_theta'] >= 0.5)]
counts.append(len(outcome))

# Row 3

outcome = df[   (df['hetB_theta'] >= 0.5) &
                        (df['hetA_theta'] < 0.5) &
                        (df['homB_theta'] < 0.5) &
                        (df['homA_theta'] < 0.5)]
counts.append(len(outcome))

outcome = df[   (df['hetB_theta'] >= 0.5) &
                        (df['hetA_theta'] < 0.5) &
                        (df['homB_theta'] < 0.5) &
                        (df['homA_theta'] >= 0.5)]
counts.append(len(outcome))

outcome = df[   (df['hetB_theta'] >= 0.5) &
                        (df['hetA_theta'] < 0.5) &
                        (df['homB_theta'] >= 0.5) &
                        (df['homA_theta'] < 0.5)]
counts.append(len(outcome))

outcome = df[   (df['hetB_theta'] >= 0.5) &
                        (df['hetA_theta'] < 0.5) &
                        (df['homB_theta'] >= 0.5) &
                        (df['homA_theta'] >= 0.5)]
counts.append(len(outcome))

# Row 4

outcome = df[   (df['hetB_theta'] >= 0.5) &
                        (df['hetA_theta'] >= 0.5) &
                        (df['homB_theta'] < 0.5) &
                        (df['homA_theta'] < 0.5)]
counts.append(len(outcome))

outcome = df[   (df['hetB_theta'] >= 0.5) &
                        (df['hetA_theta'] >= 0.5) &
                        (df['homB_theta'] < 0.5) &
                        (df['homA_theta'] >= 0.5)]
counts.append(len(outcome))

outcome = df[   (df['hetB_theta'] >= 0.5) &
                        (df['hetA_theta'] >= 0.5) &
                        (df['homB_theta'] >= 0.5) &
                        (df['homA_theta'] < 0.5)]
counts.append(len(outcome))

outcome = df[   (df['hetB_theta'] >= 0.5) &
                        (df['hetA_theta'] >= 0.5) &
                        (df['homB_theta'] >= 0.5) &
                        (df['homA_theta'] >= 0.5)]
counts.append(len(outcome))

counts = np.array(counts)
counts = counts.reshape(4, 4)

# Present results
het_labels = [[r"$\hat{\xi}_{B|A} < 1/2$", r"$\hat{\xi}_{B|A} \geq 1/2$"],
        [r"$\hat{\xi}_{A|B} < 1/2$", r"$\hat{\xi}_{A|B} \geq 1/2$"]]
hom_labels = [[r"$\hat{\lambda}_B < 1/2$", r"$\hat{\lambda}_B \geq 1/2$"],
        [r"$\hat{\lambda}_A < 1/2$", r"$\hat{\lambda}_A \geq 1/2$"]]

row_index = pd.MultiIndex.from_product(het_labels)
column_index = pd.MultiIndex.from_product(hom_labels)

print('**************************************')
print('Printing outcome table:')
print(counts)

print('**************************************')
print('Outcomes satisfying cell competition criteria:')

pd.set_option('display.max_rows', None)

print(df[ ((df['hetB_theta'] < 0.5) &
    (df['hetA_theta'] >= 0.5) &
    (df['homB_theta'] >= 0.5) &
    (df['homA_theta'] >= 0.5)) |
    ((df['hetB_theta'] >= 0.5) &
    (df['hetA_theta'] < 0.5) &
    (df['homB_theta'] >= 0.5) &
    (df['homA_theta'] >= 0.5))][['homA_simulation_id', 'homA_final_cell_count', 'homB_final_cell_count']])

#print(df[['homA_simulation_id', 'homB_simulation_id', 'het_simulation_id', 'homA_theta', 'homB_theta', 'hetA_theta', 'hetB_theta',
#    'homA_final_cell_count', 'homB_final_cell_count', 'hetA_final_cell_count', 'hetB_final_cell_count']])

# Save as figure
fig, ax = plt.subplots(figsize=(8,8))
sns.heatmap(counts, ax=ax, linewidths=1, annot=True, fmt="d")
ax.set_xlabel('Homotypic survival frequency')
ax.set_ylabel('Heterotypic survival frequency')

fig.tight_layout()
plt.savefig(image_png)
