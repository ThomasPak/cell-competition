import numpy as np
import pandas as pd
from scipy.stats import expon, uniform

import sys
sys.path.append('../scripts')
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

# Set death clock to tick at same rate as ergodic assumption
# Keep drawing G1 durations until birth and death invariants are respected
def run_g1_proportion_range_exponential_simulation_ergodic_initial_conditions(tG, eta, beta, seed=None):

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
    tbirth_0 = uniform.rvs(loc= - (tG1 + tG2), scale = tG1 + tG2,
            size=initial_cell_count, random_state=random_state)
    clone_0 = np.arange(initial_cell_count)

    # Sample G1 durations until birth invariant is satisfied.
    tG1_0 = []
    tau_0 = []
    for tbirth in tbirth_0:

        while True:

            candidate_tG1 = expon.rvs(scale=tG1, random_state=random_state)

            candidate_tau = coef * (1 - beta) * (- tbirth)

            if (- tbirth - tG2 < candidate_tG1) and \
                (not (- tbirth < candidate_tG1) or (candidate_tau < Tdeath)):

                tG1_0.append(candidate_tG1)
                tau_0.append(candidate_tau)

                break

    tG1_0 = np.array(tG1_0)
    tau_0 = np.array(tau_0)

    # Run simulation
    data = simulator.run(tau_0, tbirth_0, tG1_0, clone_0, seed=seed)

    # Return processed data
    return WellMixedSimulationData(data)

# Every cell is born at t=0 and has a fixed G1 duration equal to tG1
def run_g1_proportion_range_exponential_simulation_fixed_initial_conditions(tG, eta, beta, seed=None):

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
    tbirth_0 = np.zeros(initial_cell_count)
    clone_0 = np.arange(initial_cell_count)

    # Let tG1 be fixed
    tG1_0 = np.full(initial_cell_count, tG1)
    tau_0 = 0 * tG1_0

    # Run simulation
    data = simulator.run(tau_0, tbirth_0, tG1_0, clone_0, seed=seed)

    # Return processed data
    return WellMixedSimulationData(data)

# Every cell is born at the start of the simulation, with G1 durations drawn
# from the exponential distribution
def run_g1_proportion_range_exponential_simulation_zero_initial_conditions(tG, eta, beta, seed=None):

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
    tau_0 = np.zeros(initial_cell_count);
    tG1_0 = expon.rvs(scale=tG1, random_state=random_state, size=initial_cell_count)
    tbirth_0 = np.zeros(initial_cell_count)
    clone_0 = np.arange(initial_cell_count)

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

    df.to_csv('mc-g1-proportion-range-exponential-data.csv', index_label='simulation_id')

# Uniform ccm parameter sweep
r_fun = lambda alpha, tG: 2 * alpha * tG

# Helper function
def run_g1_proportion_range_uniform_simulation(tG, alpha, eta, beta, seed=None):

    # We create a random_state seeded with seed + 1 to sample the initial
    # conditions in order to avoid correlations with the simulation.
    if not seed is None:
        random_state = np.random.RandomState(seed + 1)
    else:
        random_state = None

    r = r_fun(alpha, tG)
    tG1 = tG1_fun(beta, tG)
    tG2 = tG2_fun(beta, tG)
    Tdeath = Tdeath_fun(eta, tG)

    ccm = uniform_ccm
    ccm_args = (tG1,r)

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
            candidate_tG1 = uniform.rvs(loc=tG1 - 0.5 * r, scale=r,
                    random_state=random_state)

        tG1_0.append(candidate_tG1)

    tG1_0 = np.array(tG1_0)

    # Run simulation
    data = simulator.run(tau_0, tbirth_0, tG1_0, clone_0, seed=seed)

    # Return processed data
    return WellMixedSimulationData(data)

if __name__ == '__main__':

    # Uniform ccm parameter sweep
    beta_minimum_dict = {
            # A
            (2/5, 3/5) : [],
            # Bp
            (1/5, 3/10) : [1 - np.sqrt(3/10)],
            # Bpp
            (2/5, 9/20) : [],
            # C
            (1/10, 3/20) : [],
            # Dp
            (2/5, 1/5) : [1 - np.sqrt(1/5)],
            # Dpp
            (3/5, 3/10) : [],
            # E
            (1/5, 1/10) : [],
            }

    tGs = np.array([100])
    alpha_etas = [(2/5, 3/5), (1/5, 3/10), (2/5, 9/20),
            (1/10, 3/20), (2/5, 1/5), (3/5, 3/10), (1/5, 1/10)]

    # Generate parameters
    tG_data = []
    alpha_data = []
    eta_data = []
    beta_data = []
    for tG in tGs:
        for alpha, eta in alpha_etas:
            beta_minimums = beta_minimum_dict[(alpha, eta)]
            betas = np.arange(alpha, 1, (1 - alpha) / num_beta)
            betas = np.append(betas, beta_minimums)

            for beta in betas:
                for i in range(num_iter):
                    tG_data.append(tG)
                    alpha_data.append(alpha)
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

    for tG, alpha, eta, beta, seed in \
            zip(tG_data, alpha_data, eta_data, beta_data, seed_data):

        sim_data = run_g1_proportion_range_uniform_simulation(tG, alpha, eta, beta, seed)

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
        'alpha' : alpha_data,
        'eta' : eta_data,
        'beta' : beta_data,
        'seed' : seed_data,
        'status' : status_data,
        'final_timestep' : final_timestep_data,
        'final_cell_count' : final_cell_count_data,
        'num_divisions' : num_divisions_data,
        'num_deaths' : num_deaths_data,
        })

    df.to_csv('mc-g1-proportion-range-uniform-data.csv', index_label='simulation_id')
