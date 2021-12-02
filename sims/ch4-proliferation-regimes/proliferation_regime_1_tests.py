import numpy as np
import pandas as pd
from scipy.stats import expon, uniform

import sys
sys.path.append('../../well_mixed')
from well_mixed_death_clock import (WellMixedSimulator,
    WellMixedSimulationData, exponential_ccm, uniform_ccm,
    base_rate_death_signal)

# Exponential cell cycle model
tG1 = 50
tG2 = 50

# Constant base rate death signal
f = base_rate_death_signal
base_rate = 1

# Simulation parameters
tstart = 0
tend = np.inf
num_iter = 1000

# Arguments to f and ccm
f_args = (base_rate,)
ccm_args = (tG1,)

# Helper function
def run_proliferation_regime_exponential_simulation(Tdeath, initial_cell_count, seed=None, max_cell_count=np.inf):

    # We create a random_state seeded with seed + 1 to sample the initial
    # conditions in order to avoid correlations with the simulation.
    if not seed is None:
        random_state = np.random.RandomState(seed + 1)
    else:
        random_state = None

    ccm = exponential_ccm

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
    thetas = np.array([ 1/5, 1/4, 1/3, 2/5, 3/7 ])
    Tdeaths = base_rate * tG1 * np.log(1 / (1 - thetas))
    initial_cell_counts = [ 10, 20, 30, 40, 50 ]

    # Generate parameters
    theta_data = []
    Tdeath_data = []
    initial_cell_count_data = []
    for theta, Tdeath in zip(thetas, Tdeaths):
        for initial_cell_count in initial_cell_counts:
            for i in range(num_iter):
                theta_data.append(theta)
                Tdeath_data.append(Tdeath)
                initial_cell_count_data.append(initial_cell_count)

    # If initial seed is given as command-line arguments, create seeds in
    # increments of 2 to avoid correlations between simulations because seed +
    # 1 is used for initial conditions.
    if len(sys.argv) == 2:
        initial_seed = int(sys.argv[1])
        seed_data = np.arange(initial_seed, initial_seed + 2 * len(Tdeath_data), 2)
    else:
        seed_data = [None] * len(Tdeath_data)

    # Run simulations and postprocess data
    status_data = []
    final_timestep_data = []
    final_cell_count_data = []
    num_divisions_data = []
    num_deaths_data = []

    for Tdeath, initial_cell_count, seed in \
            zip(Tdeath_data, initial_cell_count_data, seed_data):

        sim_data = run_proliferation_regime_exponential_simulation(Tdeath, initial_cell_count, seed)

        status = sim_data.get_status()
        t_events = sim_data.get_t_events()
        cell_count = sim_data.get_cell_count()
        num_divisions = sim_data.get_num_divisions()
        num_deaths = sim_data.get_num_deaths()

        if status == 0:
            final_timestep = t_events[-1]
        else:
            final_timestep = t_events[-2]

        final_cell_count = cell_count[-1]

        status_data.append(status)
        final_timestep_data.append(final_timestep)
        final_cell_count_data.append(final_cell_count)
        num_divisions_data.append(num_divisions)
        num_deaths_data.append(num_deaths)

    # Create and write dataframe
    df = pd.DataFrame({
        'theta' : theta_data,
        'Tdeath' : Tdeath_data,
        'initial_cell_count' : initial_cell_count_data,
        'seed' : seed_data,
        'status' : status_data,
        'final_timestep' : final_timestep_data,
        'final_cell_count' : final_cell_count_data,
        'num_divisions' : num_divisions_data,
        'num_deaths' : num_deaths_data,
        })

    df.to_csv('mc-proliferation-regime-1-data.csv', index_label='simulation_id')
