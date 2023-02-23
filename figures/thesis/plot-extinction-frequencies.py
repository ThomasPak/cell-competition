import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import pandas as pd

from formatting import set_mpl_customisation
from sys import argv

# Process arguments
if not len(argv) == 4:
    print('Usage: {} DATA_CSV CELL_BASED(vertex,well-mixed,extrusion) IMAGE_PNG'.format(argv[0]))
    exit(1)

data_csv = argv[1]
cell_based = argv[2]
image_png = argv[3]

assert cell_based == 'vertex' or cell_based == 'well-mixed' or cell_based == 'extrusion'

# Customise matplotlib
set_mpl_customisation()

# Bootstrapping parameters
Nboot = 10000
confidence = 0.95
alpha = 1 - confidence

# Seed random number generator for reproducibility (bootstrap)
np.random.seed(1)

# Import dataset
df = pd.read_csv(data_csv)

if cell_based == 'vertex' or cell_based == 'extrusion':

    # Only include completed simulations
    df = df[df['State'] == 'COMPLETED']

    # Assign well-mixed like statuses
    end_time = 100000
    threshold = 0

    df['status'] = -1

    df.loc[(df['final_timestep'] == end_time), 'status'] = 0
    df.loc[(df['final_timestep'] < end_time) & (df['final_cell_count'] == 0), 'status'] = 1
    df.loc[(df['final_timestep'] < end_time) & (df['final_cell_count'] > 0), 'status'] = 2

    df_no_extrusion = pd.DataFrame(df[df['num_extrusion'] == 0])

# Get unique theta and initial cell counts
unique_thetas = df['theta'].unique()
unique_initial_cell_counts = df['initial_cell_count'].unique()

# Status of 1 means extinction
df['went_extinct'] = (df['status'] == 1).to_numpy(dtype=int)

# Group by theta and initial cell count
grouped = df.groupby(['theta', 'initial_cell_count'])
extinction_mean  = grouped['went_extinct'].mean()
extinction_count = grouped['went_extinct'].count()
extinction_sum   = grouped['went_extinct'].sum()

# Functions
def theoretical_extinction_probability(theta, initial_cell_count, Nmax = 50):
    # return ((1 - theta) / theta)**initial_cell_count
    return ( ((1 - theta) / theta)**initial_cell_count \
            - ((1 - theta) / theta)**Nmax ) \
            / (1 - ((1 - theta) / theta)**Nmax)

def fractional_theta(theta):
    if np.isclose(theta, 4/7) or np.isclose(theta, float('{:.2f}'.format(4/7))):
        return 4, 7
    elif np.isclose(theta, 3/5):
        return 3, 5
    elif np.isclose(theta, 2/3) or np.isclose(theta, float('{:.2f}'.format(2/3))):
        return 2, 3
    elif np.isclose(theta, 3/4):
        return 3, 4
    elif np.isclose(theta, 4/5) or np.isclose(theta, float('{:.2f}'.format(4/5))):
        return 4, 5
    else:
        return 0, 0

def bootstrap_errors(num_extinctions, num_samples, extinction_frequencies):

    bootstrap_ci_errors = []

    for k, N, h in zip(num_extinctions, num_samples, extinction_frequencies):

        # reproduce data
        synthetic_data = np.zeros(N)
        synthetic_data[:k] = 1

        # Resample Nboot times and bootstrap
        bootstrap_data = np.random.choice(synthetic_data, (Nboot, N))
        bootstrap_means = np.mean(bootstrap_data, axis=1)
        bootstrap_ci = np.quantile(bootstrap_means, [alpha / 2, 1 - alpha / 2])
        bootstrap_ci_errors.append(np.abs(h - bootstrap_ci))

    return np.array(bootstrap_ci_errors).T

# Create subplots
if cell_based == 'vertex':
    fig = plt.figure(figsize=(8, 8))
    gs = fig.add_gridspec(2, 4)
    axes = []
    axes.append(fig.add_subplot(gs[0, 0:2]))
    axes.append(fig.add_subplot(gs[0, 2:]))
    axes.append(fig.add_subplot(gs[1, 1:-1]))

elif cell_based == 'well-mixed':
    fig = plt.figure(figsize=(8, 10))
    gs = fig.add_gridspec(3, 4)
    axes = []
    axes.append(fig.add_subplot(gs[0, 0:2]))
    axes.append(fig.add_subplot(gs[0, 2:]))
    axes.append(fig.add_subplot(gs[1, 0:2]))
    axes.append(fig.add_subplot(gs[1, 2:]))
    axes.append(fig.add_subplot(gs[2, 1:-1]))

# Plot
if cell_based == 'vertex' or cell_based == 'well-mixed':
    for theta, ax in zip(unique_thetas, axes):

        ii = np.arange(1, np.max(unique_initial_cell_counts) + 1)
        hh = theoretical_extinction_probability(theta, ii)
        ax.plot(ii, hh,
                linestyle='dotted',
                color='black',
                marker='o',
                zorder=-32,
                label=r'$\textrm{Predicted}$')

        num_extinctions = extinction_sum[theta].to_numpy()
        num_samples = extinction_count[theta].to_numpy()
        extinction_frequencies = extinction_mean[theta].to_numpy()

        bootstrap_ci_errors = bootstrap_errors(num_extinctions, num_samples, extinction_frequencies)

        ax.errorbar(unique_initial_cell_counts, extinction_frequencies,
                bootstrap_ci_errors,
                linestyle='none',
                color='orange',
                marker='s',
                label=r'$\textrm{Observed}$')

        # Formatting
        ax.set_xlim([0.9, 5.1])
        ax.set_ylim([-0.1, 1.1])
        ax.set_xticks(np.arange(1, 6, 1))

        # Titel, labels, legend
        num, den = fractional_theta(theta)
        ax.set_title('$\\theta = \\frac{{{}}}{{{}}}$'.format(num, den))
        ax.set_xlabel(r'$\textrm{Initial cell count}$')
        ax.set_ylabel(r'$\hat{h}: \textrm{Extinction frequency}$')
        ax.legend()

    # Save figure with tight layout
    fig.tight_layout()
    plt.savefig(image_png)

if cell_based == 'extrusion':
    ## NO EXTRUSION ##
    # Get unique theta and initial cell counts
    unique_thetas_no_extrusion = df_no_extrusion['theta'].unique()
    unique_initial_cell_counts_no_extrusion = df_no_extrusion['initial_cell_count'].unique()

    # Status of 1 means extinction
    df_no_extrusion['went_extinct'] = (df_no_extrusion['status'] == 1).to_numpy(dtype=int)

    # Group by theta and initial cell count
    grouped_no_extrusion = df_no_extrusion.groupby(['theta', 'initial_cell_count'])
    extinction_mean_no_extrusion  = grouped_no_extrusion['went_extinct'].mean()
    extinction_count_no_extrusion = grouped_no_extrusion['went_extinct'].count()
    extinction_sum_no_extrusion   = grouped_no_extrusion['went_extinct'].sum()

    fig = plt.figure(figsize=(8, 10))
    gs = fig.add_gridspec(3, 2)
    axes1 = []
    axes1.append(fig.add_subplot(gs[0, 0]))
    axes1.append(fig.add_subplot(gs[1, 0]))
    axes1.append(fig.add_subplot(gs[2, 0]))

    axes2 = []
    axes2.append(fig.add_subplot(gs[0, 1]))
    axes2.append(fig.add_subplot(gs[1, 1]))
    axes2.append(fig.add_subplot(gs[2, 1]))
    for theta, ax1, ax2 in zip(unique_thetas, axes1, axes2):

        ii = np.arange(1, np.max(unique_initial_cell_counts) + 1)
        hh = theoretical_extinction_probability(theta, ii)
        ax1.plot(ii, hh,
                linestyle='dotted',
                color='black',
                marker='o',
                zorder=-32,
                label=r'$\textrm{Predicted}$')
        ax2.plot(ii, hh,
                linestyle='dotted',
                color='black',
                marker='o',
                zorder=-32,
                label=r'$\textrm{Predicted}$')

        num_extinctions = extinction_sum[theta].to_numpy()
        num_samples = extinction_count[theta].to_numpy()
        extinction_frequencies = extinction_mean[theta].to_numpy()

        bootstrap_ci_errors = bootstrap_errors(num_extinctions, num_samples, extinction_frequencies)

        num_extinctions_no_extrusion = extinction_sum_no_extrusion[theta].to_numpy()
        num_samples_no_extrusion = extinction_count_no_extrusion[theta].to_numpy()
        extinction_frequencies_no_extrusion = extinction_mean_no_extrusion[theta].to_numpy()

        bootstrap_ci_errors_no_extrusion = bootstrap_errors(num_extinctions_no_extrusion, num_samples_no_extrusion, extinction_frequencies_no_extrusion)

        ax1.errorbar(unique_initial_cell_counts, extinction_frequencies,
                bootstrap_ci_errors,
                linestyle='none',
                color='orange',
                marker='s',
                label=r'$\textrm{Observed}$')

        ax2.errorbar(unique_initial_cell_counts_no_extrusion, extinction_frequencies_no_extrusion,
                bootstrap_ci_errors_no_extrusion,
                linestyle='none',
                color='orange',
                marker='s',
                label=r'$\textrm{Observed}$')

        # Formatting
        ax1.set_xlim([0.9, 5.1])
        ax1.set_ylim([-0.1, 1.1])
        ax1.set_xticks(np.arange(1, 6, 1))

        # Formatting
        ax2.set_xlim([0.9, 5.1])
        ax2.set_ylim([-0.1, 1.1])
        ax2.set_xticks(np.arange(1, 6, 1))

        # Titel, labels, legend
        num, den = fractional_theta(theta)
        ax1.set_title('$\\theta = \\frac{{{}}}{{{}}} \\textrm{{, full dataset}}$'.format(num, den))
        ax1.set_xlabel(r'$\textrm{Initial cell count}$')
        ax1.set_ylabel(r'$\hat{h}: \textrm{Extinction frequency}$')
        ax1.legend()

        num, den = fractional_theta(theta)
        ax2.set_title('$\\theta = \\frac{{{}}}{{{}}} \\textrm{{, no extrusions}}$'.format(num, den))
        ax2.set_xlabel(r'$\textrm{Initial cell count}$')
        ax2.set_ylabel(r'$\hat{h}: \textrm{Extinction frequency}$')
        ax2.legend()

    # Save figure with tight layout
    fig.tight_layout()
    plt.savefig(image_png)
