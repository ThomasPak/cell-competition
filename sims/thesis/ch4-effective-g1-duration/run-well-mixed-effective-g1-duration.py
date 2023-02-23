import numpy as np
import pandas as pd
from scipy.stats import expon, uniform

import sys
sys.path.append('../../../well_mixed')
from well_mixed_death_clock import (WellMixedSimulator,
    WellMixedSimulationData, exponential_ccm, uniform_ccm,
    base_rate_death_signal)

# Exponential cell cycle model
tG1 = 50
tG2 = 50

# Constant base rate death signal
f = base_rate_death_signal
base_rate = 1
Tdeath_fun = lambda eta: eta * base_rate * tG1

# Simulation parameters
tstart = 0
tend = np.inf
max_cell_count = 1000
initial_cell_count = 64
num_eta = 10
num_iter = 100

# Arguments to f and ccm
f_args = (base_rate,)
ccm_args = (tG1,)

# Helper function
def run_g1_truncation_exponential_simulation(eta, seed=None):

    # We create a random_state seeded with seed + 1 to sample the initial
    # conditions in order to avoid correlations with the simulation.
    if not seed is None:
        random_state = np.random.RandomState(seed + 1)
    else:
        random_state = None

    ccm = exponential_ccm

    Tdeath = Tdeath_fun(eta)

    # Initialise simulator
    simulator = WellMixedSimulator(f, ccm, Tdeath, tG2, tstart, tend,
            f_args, ccm_args, max_cell_count)

    # Generate initial conditions
    tau_0 = np.zeros(initial_cell_count)
    tbirth_0 = np.zeros(initial_cell_count)
    tG1_0 = expon.rvs(scale=tG1, size=initial_cell_count, random_state=random_state)
    clone_0 = np.arange(initial_cell_count)

    # Run simulation
    data = simulator.run(tau_0, tbirth_0, tG1_0, clone_0, seed=seed)

    # Return processed data
    return WellMixedSimulationData(data)

if __name__ == '__main__':

    # Exponential ccm parameter sweep
    etas = np.arange(4 / num_eta, 4 + 4 / num_eta, 4 / num_eta)

    # Generate parameters
    eta_data = []
    for eta in etas:
        for i in range(num_iter):
            eta_data.append(eta)

    # If initial seed is given as command-line arguments, create seeds in
    # increments of 2 to avoid correlations between simulations because seed +
    # 1 is used for initial conditions.
    if len(sys.argv) == 2:
        initial_seed = int(sys.argv[1])
        seed_data = np.arange(initial_seed, initial_seed + 2 * len(eta_data), 2)
    else:
        seed_data = [None] * len(eta_data)

    # Run simulations and postprocess data
    status_data = []
    final_timestep_data = []
    final_cell_count_data = []
    num_divisions_data = []
    num_deaths_data = []
    average_time_in_G1_data = []
    effective_g1_sample_size_data = []

    for eta, seed in zip(eta_data, seed_data):

        sim_data = run_g1_truncation_exponential_simulation(eta, seed)

        status = sim_data.get_status()
        t_events = sim_data.get_t_events()
        cell_count = sim_data.get_cell_count()
        num_divisions = sim_data.get_num_divisions()
        num_deaths = sim_data.get_num_deaths()
        effective_time_in_G1 = sim_data.get_effective_time_in_G1()

        if status == 0:
            final_timestep = t_events[-1]
        else:
            final_timestep = t_events[-2]

        final_cell_count = cell_count[-1]

        average_time_in_G1 = np.mean(effective_time_in_G1)
        effective_g1_sample_size = len(effective_time_in_G1)

        status_data.append(status)
        final_timestep_data.append(final_timestep)
        final_cell_count_data.append(final_cell_count)
        num_divisions_data.append(num_divisions)
        num_deaths_data.append(num_deaths)
        average_time_in_G1_data.append(average_time_in_G1)
        effective_g1_sample_size_data.append(effective_g1_sample_size)

    # Create and write dataframe
    df = pd.DataFrame({
        'eta' : eta_data,
        'seed' : seed_data,
        'status' : status_data,
        'final_timestep' : final_timestep_data,
        'final_cell_count' : final_cell_count_data,
        'num_divisions' : num_divisions_data,
        'num_deaths' : num_deaths_data,
        'average_time_in_G1' : average_time_in_G1_data,
        'effective_g1_sample_size' : effective_g1_sample_size_data,
        })

    df.to_csv('exponential-effective-g1-duration-data.csv', index_label='simulation_id')

# Uniform ccm
r_fun = lambda alpha: 2 * alpha * tG1

# Helper function
def run_g1_truncation_uniform_simulation(alpha, eta, seed=None):

    # We create a random_state seeded with seed + 1 to sample the initial
    # conditions in order to avoid correlations with the simulation.
    if not seed is None:
        random_state = np.random.RandomState(seed + 1)
    else:
        random_state = None

    ccm = uniform_ccm

    r = r_fun(alpha)
    Tdeath = Tdeath_fun(eta)

    ccm_args = (tG1,r)

    # Initialise simulator
    simulator = WellMixedSimulator(f, ccm, Tdeath, tG2, tstart, tend,
            f_args, ccm_args, max_cell_count)

    # Generate initial conditions
    tau_0 = np.zeros(initial_cell_count)
    tbirth_0 = np.zeros(initial_cell_count)
    tG1_0 = uniform.rvs(loc=tG1 - 0.5 * r, scale=r, size=initial_cell_count,
            random_state=random_state)
    clone_0 = np.arange(initial_cell_count)

    # Run simulation
    data = simulator.run(tau_0, tbirth_0, tG1_0, clone_0, seed=seed)

    # Return processed data
    return WellMixedSimulationData(data)

if __name__ == '__main__':

    # Uniform ccm parameter sweep
    alphas = [0.3, 0.5, 0.7, 1.0]
    etas = np.arange(2 / num_eta, 2 + 2 / num_eta, 2 / num_eta)

    # Generate parameters
    alpha_data = []
    eta_data = []
    for alpha in alphas:
        for eta in etas:
            for i in range(num_iter):
                alpha_data.append(alpha)
                eta_data.append(eta)

    # If initial seed is given as command-line arguments, create seeds in
    # increments of 2 to avoid correlations between simulations because seed +
    # 1 is used for initial conditions.
    if len(sys.argv) == 2:
        initial_seed = int(sys.argv[1])
        seed_data = np.arange(initial_seed, initial_seed + 2 * len(eta_data), 2)
    else:
        seed_data = [None] * len(eta_data)

    # Run simulations and postprocess data
    status_data = []
    final_timestep_data = []
    final_cell_count_data = []
    num_divisions_data = []
    num_deaths_data = []
    average_time_in_G1_data = []
    effective_g1_sample_size_data = []

    for alpha, eta, seed in zip(alpha_data, eta_data, seed_data):

        sim_data = run_g1_truncation_uniform_simulation(alpha, eta, seed)

        status = sim_data.get_status()
        t_events = sim_data.get_t_events()
        cell_count = sim_data.get_cell_count()
        num_divisions = sim_data.get_num_divisions()
        num_deaths = sim_data.get_num_deaths()
        effective_time_in_G1 = sim_data.get_effective_time_in_G1()

        if status == 0:
            final_timestep = t_events[-1]
        else:
            final_timestep = t_events[-2]

        final_cell_count = cell_count[-1]

        average_time_in_G1 = np.mean(effective_time_in_G1)
        effective_g1_sample_size = len(effective_time_in_G1)

        status_data.append(status)
        final_timestep_data.append(final_timestep)
        final_cell_count_data.append(final_cell_count)
        num_divisions_data.append(num_divisions)
        num_deaths_data.append(num_deaths)
        average_time_in_G1_data.append(average_time_in_G1)
        effective_g1_sample_size_data.append(effective_g1_sample_size)

    # Create and write dataframe
    df = pd.DataFrame({
        'alpha' : alpha_data,
        'eta' : eta_data,
        'seed' : seed_data,
        'status' : status_data,
        'final_timestep' : final_timestep_data,
        'final_cell_count' : final_cell_count_data,
        'num_divisions' : num_divisions_data,
        'num_deaths' : num_deaths_data,
        'average_time_in_G1' : average_time_in_G1_data,
        'effective_g1_sample_size' : effective_g1_sample_size_data,
        })

    df.to_csv('uniform-effective-g1-duration-data.csv', index_label='simulation_id')
