import numpy as np
import pandas as pd
from scipy.stats import expon, uniform

import sys
sys.path.append('../../../well_mixed')
from well_mixed_death_clock import (WellMixedSimulator,
    WellMixedSimulationData, exponential_ccm, uniform_ccm,
    normalised_g2_death_signal)

# Cell cycle parameters
tG1_fun = lambda beta, tG: beta * tG
tG2_fun = lambda beta, tG: (1 - beta) * tG

# normalised G2 death signal
f = normalised_g2_death_signal
coef = 1
Tdeath_fun = lambda eta, tG: eta * coef * tG

# Simulation parameters
tstart = 0
tend = np.inf
min_cell_count = 10
max_cell_count = 1000
num_iter = 100
initial_cell_count = 100
num_beta = 10

# Arguments to f
f_args = (coef,)

# Helper function
def run_g1_proportion_range_exponential_simulation(tG, eta, beta, seed=None):

    # We create a random_state seeded with seed + 1 to sample the initial
    # conditions in order to avoid correlations with the simulation.
    if not seed is None:
        random_state = np.random.RandomState(seed + 1)
    else:
        random_state = None

    tG1 = tG1_fun(beta, tG)
    tG2 = tG2_fun(beta, tG)
    Tdeath = Tdeath_fun(eta, tG)

    ccm = exponential_ccm
    ccm_args = (tG1,)

    # Initialise simulator
    simulator = WellMixedSimulator(f, ccm, Tdeath, tG2, tstart, tend,
            f_args, ccm_args, max_cell_count, min_cell_count)

    # Generate initial conditions
    tau_0 = np.zeros(initial_cell_count)
    tbirth_0 = uniform.rvs(loc= - (tG1 + tG2), scale = tG1 + tG2,
            size=initial_cell_count, random_state=random_state)
    clone_0 = np.arange(initial_cell_count)

    # Sample G1 durations until birth invariant is satisfied.
    tG1_0 = []
    for tbirth in tbirth_0:

        candidate_tG1 = - np.inf

        while not - tbirth - tG2 < candidate_tG1:
            candidate_tG1 = expon.rvs(scale=tG1, random_state=random_state)

        tG1_0.append(candidate_tG1)

    tG1_0 = np.array(tG1_0)

    # Run simulation
    data = simulator.run(tau_0, tbirth_0, tG1_0, clone_0, seed=seed)

    # Return processed data
    return WellMixedSimulationData(data)

if __name__ == '__main__':

    # Exponential ccm parameter sweep
    tGs = np.array([100])
    etas = np.array([1, 1/2, 1/5, 1/10, 1/20])
    betas = np.arange(1 / num_beta, 1, 1 / num_beta)

    # Generate parameters
    tG_data = []
    eta_data = []
    beta_data = []
    for tG in tGs:
        for eta in etas:
            for beta in betas:
                for i in range(num_iter):
                    tG_data.append(tG)
                    eta_data.append(eta)
                    beta_data.append(beta)

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

    for tG, eta, beta, seed in zip(tG_data, eta_data, beta_data, seed_data):

        sim_data = run_g1_proportion_range_exponential_simulation(tG, eta, beta, seed)

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
        'tG' : tG_data,
        'eta' : eta_data,
        'beta' : beta_data,
        'seed' : seed_data,
        'status' : status_data,
        'final_timestep' : final_timestep_data,
        'final_cell_count' : final_cell_count_data,
        'num_divisions' : num_divisions_data,
        'num_deaths' : num_deaths_data,
        })

    df.to_csv('exponential-survival-probability-data.csv', index_label='simulation_id')
