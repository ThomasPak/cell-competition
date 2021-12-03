import numpy as np
import pandas as pd
from scipy.stats import expon, uniform

import sys
sys.path.append('../../well_mixed')
from well_mixed_death_clock import (WellMixedSimulator,
    WellMixedSimulationData, exponential_ccm_heterotypic,
    uniform_ccm_heterotypic, normalised_g2_death_signal_heterotypic)

def run_heterotypic_model_2_exponential_test(tG_A, beta_A, coef_A, eta_A, n0_A,
    tG_B, beta_B, coef_B, eta_B, n0_B, min_cell_count, max_cell_count,
    tstart, tend, seed=None, min_cell_count_for_clone={}, max_cell_count_for_clone={},
    init_cond=None):

    # Cell cycle
    tG1_A = beta_A * tG_A
    tG2_A = (1 - beta_A) * tG_A

    tG1_B = beta_B * tG_B
    tG2_B = (1 - beta_B) * tG_B

    ccm_args = (tG1_A, tG1_B)

    def tG2(clone):

        assert np.all(np.logical_or(clone == 0, clone == 1))

        return tG2_A * (clone == 0) + tG2_B * (clone == 1)

    # Death signal
    f = normalised_g2_death_signal_heterotypic
    f_args = (coef_A, coef_B)

    # Death threshold
    Tdeath_A = eta_A * coef_A * tG_A
    Tdeath_B = eta_B * coef_B * tG_B

    def Tdeath(clone):

        assert np.all(np.logical_or(clone == 0, clone == 1))

        return Tdeath_A * (clone == 0) + Tdeath_B * (clone == 1)

    # Random state
    if not seed is None:
        random_state = np.random.RandomState(seed + 1)
    else:
        random_state = None

    # Ergodic death signal
    ergodic_g = (n0_A * (1 - beta_A) + n0_B * (1 - beta_B)) / (n0_A + n0_B)

    # Initial conditions
    tbirth_0_A = uniform.rvs(loc= - (tG1_A + tG2_A), scale=tG1_A + tG2_A,
            size=n0_A, random_state=random_state)
    clone_0_A = np.array([0] * n0_A)

    tG1_0_A = []
    tau_0_A = []
    for tbirth in tbirth_0_A:

            candidate_tG1 = - np.inf
            candidate_tau = 0

            if init_cond == 'ergodic':
                while (not - tbirth - tG2_A < candidate_tG1) or (not ((not (- tbirth < candidate_tG1)) or (candidate_tau < Tdeath_A))):
                    candidate_tG1 = expon.rvs(scale=tG1_A, random_state=random_state)
                    candidate_tau = coef_A * ergodic_g * candidate_tG1
            elif init_cond == None:
                while not - tbirth - tG2_A < candidate_tG1:
                    candidate_tG1 = expon.rvs(scale=tG1_A, random_state=random_state)

            tG1_0_A.append(candidate_tG1)
            tau_0_A.append(candidate_tau)

    if init_cond == 'birth':

        tbirth_0_A = np.array(tbirth_0_A) * 0
        tau_0_A = np.array(tau_0_A) * 0
        tG1_0_A = expon.rvs(scale=tG1_A, random_state=random_state, size=len(tG1_0_A))

    tG1_0_A = np.array(tG1_0_A)
    tau_0_A = np.array(tau_0_A)

    tbirth_0_B = uniform.rvs(loc= - (tG1_B + tG2_B), scale=tG1_B + tG2_B,
            size=n0_B, random_state=random_state)
    clone_0_B = np.array([1] * n0_B)

    tG1_0_B = []
    tau_0_B = []
    for tbirth in tbirth_0_B:

            candidate_tG1 = - np.inf
            candidate_tau = 0

            if init_cond == 'ergodic':
                while (not - tbirth - tG2_B < candidate_tG1) or (not ((not (- tbirth < candidate_tG1)) or (candidate_tau < Tdeath_B))):
                    candidate_tG1 = expon.rvs(scale=tG1_B, random_state=random_state)
                    candidate_tau = coef_B * ergodic_g * candidate_tG1
            elif init_cond == None:
                while not - tbirth - tG2_B < candidate_tG1:
                    candidate_tG1 = expon.rvs(scale=tG1_B, random_state=random_state)

            tG1_0_B.append(candidate_tG1)
            tau_0_B.append(candidate_tau)

    if init_cond == 'birth':

        tbirth_0_B = np.array(tbirth_0_B) * 0
        tau_0_B = np.array(tau_0_B) * 0
        tG1_0_B = expon.rvs(scale=tG1_B, random_state=random_state, size=len(tG1_0_B))

    tG1_0_B = np.array(tG1_0_B)
    tau_0_B = np.array(tau_0_B)

    tau_0 = np.concatenate([tau_0_A, tau_0_B])
    tbirth_0 = np.concatenate([tbirth_0_A, tbirth_0_B])
    clone_0 = np.concatenate([clone_0_A, clone_0_B])
    tG1_0 = np.concatenate([tG1_0_A, tG1_0_B])

    # Run simulation
    simulator = WellMixedSimulator(f, exponential_ccm_heterotypic, Tdeath, tG2,
            tstart, tend, f_args, ccm_args, max_cell_count, min_cell_count,
            min_cell_count_for_clone=min_cell_count_for_clone,
            max_cell_count_for_clone=max_cell_count_for_clone)

    raw_data = simulator.run(tau_0, tbirth_0, tG1_0, clone_0, seed=seed)

    # Return processed data
    return WellMixedSimulationData(raw_data)

def run_heterotypic_model_2_uniform_test(tG_A, beta_A, rho_A, coef_A, eta_A,
        n0_A, tG_B, beta_B, rho_B, coef_B, eta_B, n0_B, min_cell_count,
        max_cell_count, tstart, tend, seed=None, min_cell_count_for_clone={},
        max_cell_count_for_clone={}):

    # Cell cycle
    tG1_A = beta_A * tG_A
    tG2_A = (1 - beta_A) * tG_A
    r_A = 2 * tG_A * rho_A

    tG1_B = beta_B * tG_B
    tG2_B = (1 - beta_B) * tG_B
    r_B = 2 * tG_B * rho_B

    ccm_args = (tG1_A, r_A, tG1_B, r_B)

    def tG2(clone):

        assert np.all(np.logical_or(clone == 0, clone == 1))

        return tG2_A * (clone == 0) + tG2_B * (clone == 1)

    # Death signal
    f = normalised_g2_death_signal_heterotypic
    f_args = (coef_A, coef_B)

    # Death threshold
    Tdeath_A = eta_A * coef_A * tG_A
    Tdeath_B = eta_B * coef_B * tG_B

    def Tdeath(clone):

        assert np.all(np.logical_or(clone == 0, clone == 1))

        return Tdeath_A * (clone == 0) + Tdeath_B * (clone == 1)

    # Random state
    if not seed is None:
        random_state = np.random.RandomState(seed + 1)
    else:
        random_state = None

    # Initial conditions
    tau_0_A = np.zeros(n0_A)
    tbirth_0_A = uniform.rvs(loc= - (tG1_A + tG2_A), scale=tG1_A + tG2_A,
            size=n0_A, random_state=random_state)
    clone_0_A = np.array([0] * n0_A)

    tG1_0_A = []
    for tbirth in tbirth_0_A:

            candidate_tG1 = - np.inf

            while not - tbirth - tG2_A < candidate_tG1:
                candidate_tG1 = uniform.rvs(loc=tG1_A - 0.5 * r_A, scale=r_A,
                        random_state=random_state)

            tG1_0_A.append(candidate_tG1)

    tG1_0_A = np.array(tG1_0_A)

    tau_0_B = np.zeros(n0_B)
    tbirth_0_B = uniform.rvs(loc= - (tG1_B + tG2_B), scale=tG1_B + tG2_B,
            size=n0_B, random_state=random_state)
    clone_0_B = np.array([1] * n0_B)

    tG1_0_B = []
    for tbirth in tbirth_0_B:

            candidate_tG1 = - np.inf

            while not - tbirth - tG2_B < candidate_tG1:
                candidate_tG1 = uniform.rvs(loc=tG1_B - 0.5 * r_B, scale=r_B,
                        random_state=random_state)

            tG1_0_B.append(candidate_tG1)

    tG1_0_B = np.array(tG1_0_B)

    tau_0 = np.concatenate([tau_0_A, tau_0_B])
    tbirth_0 = np.concatenate([tbirth_0_A, tbirth_0_B])
    clone_0 = np.concatenate([clone_0_A, clone_0_B])
    tG1_0 = np.concatenate([tG1_0_A, tG1_0_B])

    # Run simulation
    simulator = WellMixedSimulator(f, uniform_ccm_heterotypic, Tdeath, tG2,
            tstart, tend, f_args, ccm_args, max_cell_count, min_cell_count,
            min_cell_count_for_clone=min_cell_count_for_clone,
            max_cell_count_for_clone=max_cell_count_for_clone)

    raw_data = simulator.run(tau_0, tbirth_0, tG1_0, clone_0, seed=seed)

    # Return processed data
    return WellMixedSimulationData(raw_data)

if __name__ == '__main__':

    ## Exponential
    # Cell cycle parameters
    tG_A = 100
    tG_B = 100

    beta_A = 0.3
    beta_B = 0.7

    # Death signal
    coef_A = 1
    coef_B = 1

    # Death threshold
    eta_A = .1
    eta_B = .3

    # Initial conditions
    #n0 = 500
    #n0_A = 450

    n0 = 50
    n0_A = 25

    n0_B = n0 - n0_A

    # Simulation parameters
    tstart = 0
    tend = 100
    #tend = 1000
    #tend = 1500
    #tend = np.inf
    #tend = 3000
    #seed = None
    seed = 1

    #max_cell_count = np.inf
    min_cell_count = 10
    #max_cell_count = 1000
    max_cell_count = 100

    min_cell_count_for_clone = {0 : 22}
    max_cell_count_for_clone = {1 : 35}
    #max_cell_count_for_clone = {}

    # Run simulation
    data = run_heterotypic_model_2_exponential_test(tG_A, beta_A, coef_A, eta_A, n0_A,
            tG_B, beta_B, coef_B, eta_B, n0_B, min_cell_count, max_cell_count,
            tstart, tend, seed, min_cell_count_for_clone, max_cell_count_for_clone)

    # Process data
    clone_cell_count = data.get_clone_cell_count()

    cell_count_A = clone_cell_count[:,0]
    cell_count_B = clone_cell_count[:,1]

    G1_cell_count = data.get_G1_cell_count()
    G2_cell_count = data.get_G2_cell_count()

    g = G2_cell_count / (G1_cell_count + G2_cell_count)

    ergodic_g = (cell_count_A * (1 - beta_A) + cell_count_B * (1 - beta_B)) / \
        (cell_count_A + cell_count_B)

    t = data.get_t_events()

    # Compute survival differential
    num_divisions_A = data.get_num_divisions_for_clone(0)
    num_divisions_B = data.get_num_divisions_for_clone(1)
    num_deaths_A = data.get_num_deaths_for_clone(0)
    num_deaths_B = data.get_num_deaths_for_clone(1)

    theta_A = num_divisions_A / (num_divisions_A + num_deaths_A)
    theta_B = num_divisions_B / (num_divisions_B + num_deaths_B)

    diff_theta = theta_A - theta_B

    print('theta_A, theta_B: {}, {}'.format(theta_A, theta_B))
    print('diff_theta: {}'.format(diff_theta))

    # Plot data
    import matplotlib.pyplot as plt
    plt.rc('text', usetex=True)

    fig, ax1 = plt.subplots(figsize=(10, 5))

    # Cell counts
    ax1.plot(t, cell_count_A, label='A')
    ax1.plot(t, cell_count_B, label='B')

    ax1.set_xlabel('t')
    ax1.set_ylabel('cell count')

    ax1.set_title(r'''$A: t_G = {}, \eta = {}, \beta = {}, n_0 = {}$
    $B: t_G = {}, \eta = {}, \beta = {}, n_0 = {}$'''.format(
        tG_A,
        eta_A,
        beta_A,
        n0_A,
        tG_B,
        eta_B,
        beta_B,
        n0_B))

    ax1.legend(loc='best')

    # Death signal
    ax2 = ax1.twinx()

    ax2.plot(t, g, 'k', label='g(t)')
    ax2.plot(t, ergodic_g, 'r', label='ergodic')

    ax2.set_ylabel('death signal')

    ax2.set_ylim([0, 1])

    ax2.legend(loc='best')

    # Save figure
    fig.tight_layout()
    fig.savefig('heterotypic_model_2_exponential.png')

    ## Uniform
    # Cell cycle parameters
    tG_A = 100
    tG_B = 100
    #tG_A = 1

    beta_A = 0.3
    beta_B = 0.7

    rho_A = 0.3
    rho_B = 0.3

    # Death signal
    coef_A = 1
    coef_B = 1

    # Death threshold
    eta_A = 0.15
    eta_B = 0.3

    # Initial conditions
    n0 = 100
    n0_A = 30

    n0 = 50
    n0_A = 25

    n0_B = n0 - n0_A

    # Simulation parameters
    tstart = 0
    tend = 100
    #tend = 1000
    #tend = 1500
    #tend = np.inf
    #tend = 3000
    #tend = 500
    #seed = None
    seed = 1

    #max_cell_count = np.inf
    min_cell_count = 10
    #max_cell_count = 1000
    max_cell_count = 100

    min_cell_count_for_clone = {}
    max_cell_count_for_clone = {1 : 42, 0 : 37}

    # Run simulation
    data = run_heterotypic_model_2_uniform_test(tG_A, beta_A, rho_A, coef_A, eta_A,
            n0_A, tG_B, beta_B, rho_B, coef_B, eta_B, n0_B, min_cell_count,
            max_cell_count, tstart, tend, seed, min_cell_count_for_clone,
            max_cell_count_for_clone)

    # Process data
    clone_cell_count = data.get_clone_cell_count()

    cell_count_A = clone_cell_count[:,0]
    cell_count_B = clone_cell_count[:,1]

    G1_cell_count = data.get_G1_cell_count()
    G2_cell_count = data.get_G2_cell_count()

    g = G2_cell_count / (G1_cell_count + G2_cell_count)

    ergodic_g = (cell_count_A * (1 - beta_A) + cell_count_B * (1 - beta_B)) / \
        (cell_count_A + cell_count_B)

    t = data.get_t_events()

    # Compute survival differential
    num_divisions_A = data.get_num_divisions_for_clone(0)
    num_divisions_B = data.get_num_divisions_for_clone(1)
    num_deaths_A = data.get_num_deaths_for_clone(0)
    num_deaths_B = data.get_num_deaths_for_clone(1)

    theta_A = num_divisions_A / (num_divisions_A + num_deaths_A)
    theta_B = num_divisions_B / (num_divisions_B + num_deaths_B)

    diff_theta = theta_A - theta_B

    print('theta_A, theta_B: {}, {}'.format(theta_A, theta_B))
    print('diff_theta: {}'.format(diff_theta))

    # Plot data
    import matplotlib.pyplot as plt
    plt.rc('text', usetex=True)

    fig, ax1 = plt.subplots(figsize=(10, 5))

    # Cell counts
    ax1.plot(t, cell_count_A, label='A')
    ax1.plot(t, cell_count_B, label='B')

    ax1.set_xlabel('t')
    ax1.set_ylabel('cell count')

    ax1.set_title(r'''$A: t_G = {}, \eta = {}, \beta = {}, \rho = {}, n_0 = {}$
    $B: t_G = {}, \eta = {}, \beta = {}, \rho = {}, n_0 = {}$'''.format(
        tG_A,
        eta_A,
        beta_A,
        rho_A,
        n0_A,
        tG_B,
        eta_B,
        beta_B,
        rho_B,
        n0_B))

    ax1.legend(loc='best')

    # Death signal
    ax2 = ax1.twinx()

    ax2.plot(t, g, 'k', label='g(t)')
    ax2.plot(t, ergodic_g, 'r', label='ergodic')

    ax2.set_ylabel('death signal')

    ax2.set_ylim([0, 1])

    ax2.legend(loc='best')

    # Save figure
    fig.tight_layout()
    fig.savefig('heterotypic_model_2_uniform.png')
