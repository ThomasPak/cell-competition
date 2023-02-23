#!/usr/bin/env python
# coding: utf-8

# We use the Python package pandas to read a comma-separated file containing our data
# Pandas project: https://pandas.pydata.org/
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl
import seaborn as sns

from os.path import basename
from sys import argv

valid_image_types = '''
Valid image types:     psi-eta-histogram
                       b-cells-histogram
                       all-cells-histogram
                       diff-histogram
                       theta-histogram
                       diff-theta-histogram
                       boxplots
                       histograms
                       regime-count
                       mcc-found-parameters
                       pearson-correlation
                       spearman-correlation
'''

# Process arguments
if len(argv) != 5:
    print('Usage: {} DATA_HETERO_CSV DATA_HOMO_CSV IMAGE_TYPE IMAGE_PNG'.format(argv[0]))
    print(valid_image_types)
    exit(1)

data_hetero_csv = argv[1]
data_homo_csv = argv[2]
image_type = argv[3]
image_png = argv[4]

# Process data
dfhet = pd.read_csv(data_hetero_csv)
dfhom = pd.read_csv(data_homo_csv)

dfhet['total-num-cells'] = dfhet['final_cell_count']
dfhet['num-a-cells'] = dfhet['final_cell_count_0']
dfhet['num-b-cells'] = dfhet['final_cell_count_1']
dfhet['b-proportion'] = dfhet['num-b-cells'] / dfhet['total-num-cells']

## Reorder columns so that ab-line tension comes after b-line tension
#dfhet = dfhet.reindex(columns=['simulation-id',
#                                           'seed',
#                                           'initial-b-ratio',
#                                           'random-labelling',
#                                           'cell-b-target-area',
#                                           'cell-b-elasticity-parameter',
#                                           'cell-b-contractility-parameter',
#                                           'cell-b-line-tension-parameter',
#                                           'cell-ab-line-tension-parameter',
#                                           'cell-b-g1-duration',
#                                           'cell-b-g2-duration',
#                                           'cell-b-cell-cycle-duration',
#                                           'num-a-cells',
#                                           'num-b-cells',
#                                           'total-num-cells',
#                                           'b-proportion'])

dfhet['control-num-b-cells'] = dfhom['final_cell_count']
dfhet['control-simulation-id'] = dfhom['simulation-id']

dfhet['control_final_cell_count'] = dfhom['final_cell_count']
dfhet['control_initial_cell_count'] = dfhom['initial_cell_count']
dfhet['control_num_apoptosis'] = dfhom['num_apoptosis']
dfhet['control_num_extrusion'] = dfhom['num_extrusion']

dfhet['cell-b-cell-cycle-duration'] = dfhet['cell-b-g1-duration'] + dfhet['cell-b-g2-duration']

df_all = dfhet.drop('initial-b-ratio', axis=1)
#df_all = dfhet

# Compute competition quotient
df_all['competition-quotient'] = 2 * df_all['num-b-cells'] / (2 * df_all['num-b-cells'] + df_all['control-num-b-cells'])

# Describe summary statistics
print('Describe summary statistics')
print(df_all[['num-a-cells',
        'num-b-cells',
        'total-num-cells',
        'b-proportion',
        'control-num-b-cells',
        'competition-quotient']].describe())

# Print extreme ends of competition quotient
print('extreme ends of competition quotient')
print(df_all[['num-a-cells',
        'num-b-cells',
        'total-num-cells',
        'b-proportion',
        'control-num-b-cells',
        'competition-quotient']].dropna().sort_values(by='competition-quotient'))

# Print Pearson correlation table
print('Pearson correlation table')
print(df_all.drop(['seed', 'simulation-id'], axis=1).corr(method='pearson')[
       ['num-a-cells',
       'num-b-cells',
       'total-num-cells',
       'b-proportion',
       'control-num-b-cells',
       'competition-quotient']])

# Print Spearman correlation table
print('Spearman correlation table')
print(df_all.drop(['seed', 'simulation-id'], axis=1).corr(method='spearman')[
        ['num-a-cells',
        'num-b-cells',
        'total-num-cells',
        'b-proportion',
        'control-num-b-cells',
        'competition-quotient']])

# Set customisation
from formatting import set_mpl_customisation
set_mpl_customisation()

# Create figure
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(8,4))

# Plot histogram of psi and eta
ax1.hist(df_all['b-proportion'], ec='black')
ax2.hist(df_all['competition-quotient'], ec='black')

# Formatting
ax1.set_xlabel(r'$\psi$')
ax1.set_ylabel(r'$\textrm{Count}$')
ax1.set_title(r'$\textrm{Histogram of } \psi$')

ax1.set_xlim([0, 1])
ax1.set_xticks([0, 0.25, 0.5, 0.75, 1])

ax2.set_xlabel(r'$\eta$')
ax2.set_ylabel(r'$\textrm{Count}$')
ax2.set_title(r'$\textrm{Histogram of } \eta$')

ax2.set_xlim([0, 1])
ax2.set_xticks([0, 0.25, 0.5, 0.75, 1])

if image_type == 'psi-eta-histogram':
    fig.tight_layout()
    fig.savefig(image_png)
    exit()

# Plot histograms of B-type cell count
# Create figure
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(8,4))

#_, bins, _ = ax1.hist(df_all['num-b-cells'] * 2, ec='black')
#ax2.hist(df_all['control-num-b-cells'], ec='black', bins=bins)

n1, bins, _ = ax1.hist(df_all['num-b-cells'] * 2, ec='black')
import numpy as np
#bins = np.linspace(0, 420, 11)
n2, _, _ = ax2.hist(df_all['control-num-b-cells'], ec='black', bins=bins)

# Formatting
ax1.set_xlabel(r'$2N_B$')
ax1.set_ylabel(r'$\textrm{Count}$')
ax1.set_title(r'$\textrm{Histogram of } 2N_B$')

ax2.set_xlabel(r'$N_{B,\textrm{control}}$')
ax2.set_ylabel(r'$\textrm{Count}$')
ax2.set_title(r'$\textrm{Histogram of } N_{B,\textrm{control}}$')

if image_type == 'b-cells-histogram':
    fig.tight_layout()
    plt.savefig(image_png)
    exit()

# Plot histograms of total cell count
hist = df_all[['num-a-cells', 'num-b-cells', 'total-num-cells']].hist(figsize=(12, 3), layout=(1,3))

if image_type == 'all-cells-histogram':
    fig.tight_layout()
    plt.savefig(image_png)
    exit()

# Compute heterotypic survival difference
df = df_all
df['num_divisions_0'] = df['final_cell_count_0'] - df['initial_cell_count_0'] + \
        df['num_apoptosis_0'] + df['num_extrusion_0']
df['num_divisions_1'] = df['final_cell_count_1'] - df['initial_cell_count_1'] + \
        df['num_apoptosis_1'] + df['num_extrusion_1']

assert np.all(df['num_divisions_0'] >= 0)
assert np.all(df['num_divisions_1'] >= 0)

df['num_deaths_0'] = df['num_apoptosis_0'] + df['num_extrusion_0']
df['num_deaths_1'] = df['num_apoptosis_1'] + df['num_extrusion_1']

df['theta_0'] = df['num_divisions_0'] / (df['num_divisions_0'] + df['num_deaths_0'])
df['theta_1'] = df['num_divisions_1'] / (df['num_divisions_1'] + df['num_deaths_1'])

df['hetdiff-theta'] = df['theta_1'] - df['theta_0']

# Compute homotypic survival difference B
df['control_num_divisions'] = df['control_final_cell_count'] - \
        df['control_initial_cell_count'] + df['control_num_apoptosis'] + \
        df['control_num_extrusion']

assert np.all(df['control_num_divisions'] >= 0)

df['control_num_deaths'] = df['control_num_apoptosis'] + df['control_num_extrusion']

df['control_theta'] = df['control_num_divisions'] / (df['control_num_divisions'] + df['control_num_deaths'])

df['homdiff-theta-1'] = df['theta_1'] - df['control_theta']

# Compute homotypic survival difference A
dfhomA = df[
   (df['cell-b-target-area'] == 1) &
   (df['cell-b-elasticity-parameter'] == 1) &
   (df['cell-b-contractility-parameter'] == 0.04) &
   (df['cell-b-line-tension-parameter'] == 0.12) &
   (df['cell-ab-line-tension-parameter'] == 0.12) &
   (df['cell-b-g1-duration'] == 30) &
   (df['cell-b-g2-duration'] == 70)]

print('lambda A: {}'.format(dfhomA['control_theta'].mean()))

if len(dfhomA) > 0:
    df['homdiff-theta-0'] = df['theta_0'] - dfhomA['control_theta'].mean()
else:
    df['homdiff-theta-0'] = df['theta_0'] - 0.9942857142857142

# Plot histograms of heterotypic and homotypic survival differences
# Create figure
fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(8,3))

# Create histograms
bins = np.linspace(-1, 1, 40)
ax1.hist(df_all['hetdiff-theta'], ec='black', bins=bins)
ax2.hist(df_all['homdiff-theta-1'], ec='black', bins=bins)
ax3.hist(df_all['homdiff-theta-0'], ec='black', bins=bins)

# Formatting
ax1.set_xlabel(r'$\Delta^{\neq}_{A|B}$')
ax1.set_ylabel(r'$\textrm{Count}$')
ax1.set_title(r'$\textrm{Histogram of } \Delta^{\neq}_{A|B}$')

ax2.set_xlabel(r'$\Delta^=_{A|B}$')
ax2.set_ylabel(r'$\textrm{Count}$')
ax2.set_title(r'$\textrm{Histogram of } \Delta^=_{A|B}$')

ax3.set_xlabel(r'$\Delta^=_{B|A}$')
ax3.set_ylabel(r'$\textrm{Count}$')
ax3.set_title(r'$\textrm{Histogram of } \Delta^=_{B|A}$')

if image_type == 'diff-histogram':
    fig.tight_layout()
    plt.savefig(image_png)
    exit()

# Plot histograms of homotypic and heterotypic survival frequency
# Create figure
fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(8,3))

# Create histograms
bins = np.linspace(0, 1, 20)
ax1.hist(df_all['control_theta'], ec='black', bins=bins)
ax2.hist(df_all['theta_1'], ec='black', bins=bins)
ax3.hist(df_all['theta_0'], ec='black', bins=bins)

# Formatting
ax1.set_xlabel(r'$\lambda_A$')
ax1.set_ylabel(r'$\textrm{Count}$')
ax1.set_title(r'$\textrm{Histogram of } \lambda_A$')

ax2.set_xlabel(r'$\xi_{A|B}$')
ax2.set_ylabel(r'$\textrm{Count}$')
ax2.set_title(r'$\textrm{Histogram of } \xi_{A|B}$')

ax3.set_xlabel(r'$\xi_{B|A}$')
ax3.set_ylabel(r'$\textrm{Count}$')
ax3.set_title(r'$\textrm{Histogram of } \xi_{B|A}$')

if image_type == 'theta-histogram':
    fig.tight_layout()
    plt.savefig(image_png)
    exit()

# Plot histograms of homotypic and heterotypic survival frequency and difference
# Create figure
fig, ((ax1, ax2, ax3), (ax4, ax5, ax6)) = plt.subplots(2, 3, figsize=(8,6))

# Create histograms
Nbins = 30
bins = np.linspace(0, 1, Nbins)
ax1.hist(df_all['control_theta'], ec='black', bins=bins)
ax2.hist(df_all['theta_1'], ec='black', bins=bins)
ax3.hist(df_all['theta_0'], ec='black', bins=bins)

bins = np.linspace(-1, 1, Nbins)
ax4.hist(df_all['hetdiff-theta'], ec='black', bins=bins)
ax5.hist(df_all['homdiff-theta-1'], ec='black', bins=bins)
ax6.hist(df_all['homdiff-theta-0'], ec='black', bins=bins)

# Formatting
ax1.set_xlabel(r'$\hat{\lambda}_A$')
ax1.set_ylabel(r'$\textrm{Count}$')
ax1.set_title(r'$\textrm{Histogram of } \hat{\lambda}_A$')

ax2.set_xlabel(r'$\hat{\xi}_{A|B}$')
ax2.set_ylabel(r'$\textrm{Count}$')
ax2.set_title(r'$\textrm{Histogram of } \hat{\xi}_{A|B}$')

ax3.set_xlabel(r'$\hat{\xi}_{B|A}$')
ax3.set_ylabel(r'$\textrm{Count}$')
ax3.set_title(r'$\textrm{Histogram of } \hat{\xi}_{B|A}$')

# Formatting
ax4.set_xlabel(r'$\hat{\Delta}^{\neq}_{A|B}$')
ax4.set_ylabel(r'$\textrm{Count}$')
ax4.set_title(r'$\textrm{Histogram of } \hat{\Delta}^{\neq}_{A|B}$')

ax5.set_xlabel(r'$\hat{\Delta}^=_{A|B}$')
ax5.set_ylabel(r'$\textrm{Count}$')
ax5.set_title(r'$\textrm{Histogram of } \hat{\Delta}^=_{A|B}$')

ax6.set_xlabel(r'$\hat{\Delta}^=_{B|A}$')
ax6.set_ylabel(r'$\textrm{Count}$')
ax6.set_title(r'$\textrm{Histogram of } \hat{\Delta}^=_{B|A}$')

if image_type == 'diff-theta-histogram':
    fig.tight_layout()
    plt.savefig(image_png)
    exit()

# Plot boxplots
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(8, 4))

ax1.boxplot([df['control_theta'],
    df['theta_1'],
    df['theta_0']],
    whis=(0, 100),
    labels=[r'$\lambda$', r'$\xi_{A|B}$', r'$\xi_{B|A}$'])

#ax1.set_ylim([-0.1, 1.1])

ax2.boxplot(df[['hetdiff-theta', 'homdiff-theta-1', 'homdiff-theta-0']],
        whis=(0, 100),
        labels=[r'$\Delta^{\neq}_{A|B}$', r'$\Delta^{=}_{A|B}$', r'$\Delta^{=}_{B|A}$'])

#ax2.set_ylim([-1, 1])

if image_type == 'boxplots':
    fig.tight_layout()
    plt.savefig(image_png)
    exit()

# Plot more histograms
fig, axes = plt.subplots(3, 2, figsize=(4, 6))

ax1 = axes[0, 0]
ax2 = axes[0, 1]
ax3 = axes[1, 0]
ax4 = axes[1, 1]
ax5 = axes[2, 0]
ax6 = axes[2, 1]

ax1.hist(df_all['control_theta'], ec='black')
ax2.hist(df_all['theta_1'], ec='black')
ax3.hist(df_all['theta_0'], ec='black')
ax4.hist(df_all['hetdiff-theta'], ec='black')
ax5.hist(df_all['homdiff-theta-1'], ec='black')
ax6.hist(df_all['homdiff-theta-0'], ec='black')

#ax2.set_ylim([-1, 1])

if image_type == 'histograms':
    fig.tight_layout()
    plt.savefig(image_png)
    exit()

# Apply competition criteria to cell types and count number of occurrences.
# intrinsically nonviable
ccdf = df[(df['control_theta'] < 0.5) & (df['theta_1'] < 0.5)]
homo_nonviable_hetero_nonviable_num = len(ccdf)

# Loser cell rescue
ccdf = df[(df['control_theta'] < 0.5) & (df['theta_1'] >= 0.5)]
homo_nonviable_hetero_viable_num = len(ccdf)

# non-cell competition
ccdf = df[(df['control_theta'] >= 0.5) & (df['theta_1'] >= 0.5)]
homo_viable_hetero_viable_num = len(ccdf)

# cell competition
ccdf = df[(df['control_theta'] >= 0.5) & (df['theta_1'] < 0.5)]
homo_viable_hetero_nonviable_num = len(ccdf)

dfnew = pd.DataFrame(
        {
            r'$\hat{\lambda} < 1/2$' : [homo_nonviable_hetero_nonviable_num,
            homo_nonviable_hetero_viable_num] ,
            r'$\hat{\lambda} \geq 1/2$' : [homo_viable_hetero_nonviable_num,
            homo_viable_hetero_viable_num]
            })
dfnew.index = [r'$\hat{\xi}_{A|B} < 1/2$', r'$\hat{\xi}_{A|B} \geq 1/2$']

fig, ax = plt.subplots(figsize=(4,5))
sns.heatmap(dfnew, ax=ax, linewidths=1, annot=True, fmt="d")

if image_type == 'regime-count':
    fig.tight_layout()
    plt.savefig(image_png)
    exit()

df['homo_nonviable_hetero_nonviable'] = (df['control_theta'] < 0.5) & (df['theta_1'] < 0.5)
df['homo_nonviable_hetero_viable'] = (df['control_theta'] < 0.5) & (df['theta_1'] >= 0.5)
df['homo_viable_hetero_nonviable'] = (df['control_theta'] >= 0.5) & (df['theta_1'] < 0.5)
df['homo_viable_hetero_viable'] = (df['control_theta'] >= 0.5) & (df['theta_1'] >= 0.5)

grouped = df.groupby(['random-labelling', 'cell-b-target-area',
    'cell-b-elasticity-parameter', 'cell-b-contractility-parameter',
    'cell-b-line-tension-parameter', 'cell-ab-line-tension-parameter',
    'cell-b-g1-duration', 'cell-b-g2-duration', 'cell-b-cell-cycle-duration'])

#from IPython import embed; embed()


print('found cell competition parameters')
ccdfparam = ccdf[['simulation-id',
                                           'seed',
                                           'random-labelling',
                                           'cell-b-target-area',
                                           'cell-b-elasticity-parameter',
                                           'cell-b-contractility-parameter',
                                           'cell-b-line-tension-parameter',
                                           'cell-ab-line-tension-parameter',
                                           'cell-b-g1-duration',
                                           'cell-b-g2-duration',
                                           'cell-b-cell-cycle-duration' ]]
print(ccdfparam)

if image_type == 'mcc-found-parameters':
    ccdfparam.to_csv(image_png)
    exit()

# Plot correlation heatmap
df.rename(columns={
    'random-labelling' : r'$\textrm{Labelling}$',
    'cell-b-target-area' : r'$S^0_A$',
    'cell-b-elasticity-parameter' : r'$K_A$',
    'cell-b-contractility-parameter' : r'$\Gamma_A$',
    'cell-b-line-tension-parameter' : r'$\Lambda_A$',
    'cell-ab-line-tension-parameter': r'$\Lambda_{AB}$',
    'cell-b-g1-duration' : r'$t_{\textrm{G1},A}$',
    'cell-b-g2-duration' : r'$t_{\textrm{G2},A}$',
    'control_theta' : r'$\hat{\lambda}_A$',
    'theta_1' : r'$\hat{\xi}_{A|B}$',
    'theta_0' : r'$\hat{\xi}_{B|A}$',
    'hetdiff-theta' : r'$\hat{\Delta}^{\neq}_{A|B}$',
    'homdiff-theta-1' : r'$\hat{\Delta}^{=}_{A|B}$',
    'homdiff-theta-0' : r'$\hat{\Delta}^{=}_{B|A}$',
    }, inplace=True)

pearson_corr = df[[
    r'$\textrm{Labelling}$',
    r'$S^0_A$',
    r'$K_A$',
    r'$\Gamma_A$',
    r'$\Lambda_A$',
    r'$\Lambda_{AB}$',
    r'$t_{\textrm{G1},A}$',
    r'$t_{\textrm{G2},A}$',
    r'$\hat{\lambda}_A$',
    r'$\hat{\xi}_{A|B}$',
    r'$\hat{\xi}_{B|A}$',
    r'$\hat{\Delta}^{\neq}_{A|B}$',
    r'$\hat{\Delta}^{=}_{A|B}$',
    r'$\hat{\Delta}^{=}_{B|A}$',
    ]].corr(method='pearson')[[
    r'$\hat{\lambda}_A$',
    r'$\hat{\xi}_{A|B}$',
    r'$\hat{\xi}_{B|A}$',
    r'$\hat{\Delta}^{\neq}_{A|B}$',
    r'$\hat{\Delta}^{=}_{A|B}$',
    r'$\hat{\Delta}^{=}_{B|A}$',
        ]].drop([
    r'$\hat{\lambda}_A$',
    r'$\hat{\xi}_{A|B}$',
    r'$\hat{\xi}_{B|A}$',
    r'$\hat{\Delta}^{\neq}_{A|B}$',
    r'$\hat{\Delta}^{=}_{A|B}$',
    r'$\hat{\Delta}^{=}_{B|A}$',
        ])

#pearson_corr = df[[
#    r'$\textrm{Labelling}$',
#    r'$S^0_A$',
#    r'$K_A$',
#    r'$\Gamma_A$',
#    r'$\Lambda_A$',
#    r'$\Lambda_{AB}$',
#    r'$t_{1,A}$',
#    r'$t_{2,A}$',
#    r'$\lambda$',
#    r'$\xi_{A|B}$',
#    r'$\xi_{B|A}$',
#    r'$\Delta^{\neq}_{A|B}$',
#    r'$\Delta^{=}_{A|B}$',
#    r'$\Delta^{=}_{B|A}$',
#    ]].corr(method='pearson')[[
#    r'$\lambda$',
#    r'$\xi_{A|B}$',
#    r'$\xi_{B|A}$',
#    r'$\Delta^{\neq}_{A|B}$',
#    r'$\Delta^{=}_{A|B}$',
#    r'$\Delta^{=}_{B|A}$',
#        ]]
#
fig, ax = plt.subplots(figsize=(4, 3))
#fig, ax = plt.subplots(figsize=(4, 5))

sns.heatmap(pearson_corr, ax=ax, linewidths=1, center=0)

if image_type == 'pearson-correlation':
    fig.tight_layout()
    plt.savefig(image_png)
    exit()

spearman_corr = df[[
    r'$\textrm{Labelling}$',
    r'$S^0_A$',
    r'$K_A$',
    r'$\Gamma_A$',
    r'$\Lambda_A$',
    r'$\Lambda_{AB}$',
    r'$t_{\textrm{G1},A}$',
    r'$t_{\textrm{G2},A}$',
    r'$\hat{\lambda}_A$',
    r'$\hat{\xi}_{A|B}$',
    r'$\hat{\xi}_{B|A}$',
    r'$\hat{\Delta}^{\neq}_{A|B}$',
    r'$\hat{\Delta}^{=}_{A|B}$',
    r'$\hat{\Delta}^{=}_{B|A}$',
    ]].corr(method='spearman')[[
    r'$\hat{\lambda}_A$',
    r'$\hat{\xi}_{A|B}$',
    r'$\hat{\xi}_{B|A}$',
    r'$\hat{\Delta}^{\neq}_{A|B}$',
    r'$\hat{\Delta}^{=}_{A|B}$',
    r'$\hat{\Delta}^{=}_{B|A}$',
        ]].drop([
    r'$\hat{\lambda}_A$',
    r'$\hat{\xi}_{A|B}$',
    r'$\hat{\xi}_{B|A}$',
    r'$\hat{\Delta}^{\neq}_{A|B}$',
    r'$\hat{\Delta}^{=}_{A|B}$',
    r'$\hat{\Delta}^{=}_{B|A}$',
        ])

fig, ax = plt.subplots(figsize=(4, 3))

sns.heatmap(spearman_corr, ax=ax, linewidths=1, center=0)

if image_type == 'spearman-correlation':
    fig.tight_layout()
    plt.savefig(image_png)
    exit()

# Error if reach here
print('Error: {}'.format(valid_image_types))
exit(2)
