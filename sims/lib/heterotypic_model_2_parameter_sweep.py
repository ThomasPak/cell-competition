import numpy as np
import pandas as pd
import sys

from heterotypic_model_2_test import (run_heterotypic_model_2_exponential_test,
        run_heterotypic_model_2_uniform_test)

from well_mixed_death_clock import normalised_g2_death_signal_heterotypic

def heterotypic_model_2_parameter_sweep(beta_As, tG_As, eta_As, coef_As,
        beta_Bs, tG_Bs, eta_Bs, coef_Bs, n0s, datafile_prefix, ccm, num_iter,
        tstart, tend, min_cell_count, max_cell_count, seed=None, sim_range=None,
        min_cell_count_for_clone={}, max_cell_count_for_clone={}):
    """
    Implicitly regularises uniform ccm by setting rho_A = beta_A and rho_B =
    beta_B.
    """

    # Generate parameters
    beta_A_data = []
    tG_A_data = []
    eta_A_data = []
    coef_A_data = []

    beta_B_data = []
    tG_B_data = []
    eta_B_data = []
    coef_B_data = []

    n0_A_data = []
    n0_B_data = []

    for beta_A in beta_As:
     for tG_A in tG_As:
      for eta_A in eta_As:
       for coef_A in coef_As:
        for beta_B in beta_Bs:
         for tG_B in tG_Bs:
          for eta_B in eta_Bs:
           for coef_B in coef_Bs:
            for n0 in n0s:
             n0_A = n0[0]
             n0_B = n0[1]
             for i in range(num_iter):
              beta_A_data.append(beta_A)
              tG_A_data.append(tG_A)
              eta_A_data.append(eta_A)
              coef_A_data.append(coef_A)

              beta_B_data.append(beta_B)
              tG_B_data.append(tG_B)
              eta_B_data.append(eta_B)
              coef_B_data.append(coef_B)

              n0_A_data.append(n0_A)
              n0_B_data.append(n0_B)

    rho_A_data = beta_A_data
    rho_B_data = beta_B_data

    heterotypic_model_2_parameter_sweep_general(beta_A_data, tG_A_data,
            rho_A_data, eta_A_data, coef_A_data, beta_B_data, tG_B_data,
            rho_B_data, eta_B_data, coef_B_data, n0_A_data, n0_B_data,
            datafile_prefix, ccm, num_iter, tstart, tend, min_cell_count,
            max_cell_count, seed, sim_range,
            min_cell_count_for_clone=min_cell_count_for_clone,
            max_cell_count_for_clone=max_cell_count_for_clone)

def heterotypic_model_2_parameter_sweep_general(beta_A_data, tG_A_data,
    rho_A_data, eta_A_data, coef_A_data, beta_B_data, tG_B_data, rho_B_data,
    eta_B_data, coef_B_data, n0_A_data, n0_B_data, datafile_prefix, ccm,
    num_iter, tstart, tend, min_cell_count, max_cell_count, seed=None,
    sim_range=None, min_cell_count_for_clone={}, max_cell_count_for_clone={}):

    assert (len(beta_A_data) == len(tG_A_data) == len(rho_A_data) ==
            len(eta_A_data) == len(coef_A_data) == len(beta_B_data) ==
            len(tG_B_data) == len(rho_B_data) == len(eta_B_data) ==
            len(coef_B_data) == len(n0_A_data) == len(n0_B_data))

    # If initial seed is given, create seeds in increments of 2 to avoid
    # correlations between simulations because seed + 1 is used for initial
    # conditions.
    if not seed is None:
        initial_seed = seed
        seed_data = np.arange(initial_seed, initial_seed + 2 * len(beta_A_data), 2)
    else:
        seed_data = [None] * len(beta_A_data)

    # Default simulation range is 0 and len(data) - 1
    if not sim_range is None:
        start_sim_num = sim_range[0]
        end_sim_num = sim_range[1]
    else:
        start_sim_num = 0
        end_sim_num = len(beta_A_data)

    print('Total number of simulations: {}'.format(len(beta_A_data)))
    print('Current simulation range: [{}, {})'.format(start_sim_num, end_sim_num))

    # Run simulations and postprocess data
    status_data = []
    status_info_data = []
    final_timestep_data = []
    final_cell_count_A_data = []
    final_cell_count_B_data = []
    num_divisions_A_data = []
    num_divisions_B_data = []
    num_deaths_A_data = []
    num_deaths_B_data = []
    ergodic_rms_A_data = []
    ergodic_rms_B_data = []

    for i, (beta_A, tG_A, rho_A, eta_A, coef_A, beta_B, tG_B, rho_B, eta_B,
            coef_B, n0_A, n0_B, seed) in enumerate(zip(
            beta_A_data, tG_A_data, rho_A_data, eta_A_data, coef_A_data,
            beta_B_data, tG_B_data, rho_B_data, eta_B_data, coef_B_data,
            n0_A_data, n0_B_data, seed_data)):

        if not (start_sim_num <= i and i < end_sim_num):
            continue

        if ccm == 'exponential':
            sim_data = run_heterotypic_model_2_exponential_test(tG_A, beta_A, coef_A, eta_A, n0_A,
                    tG_B, beta_B, coef_B, eta_B, n0_B, min_cell_count, max_cell_count,
                    tstart, tend, seed, min_cell_count_for_clone=min_cell_count_for_clone,
                    max_cell_count_for_clone=max_cell_count_for_clone)

        elif ccm == 'uniform':
            sim_data = run_heterotypic_model_2_uniform_test(tG_A, beta_A, rho_A,
                    coef_A, eta_A, n0_A, tG_B, beta_B, rho_B, coef_B, eta_B, n0_B,
                    min_cell_count, max_cell_count, tstart, tend, seed,
                    min_cell_count_for_clone=min_cell_count_for_clone,
                    max_cell_count_for_clone=max_cell_count_for_clone)

        else:
            raise Exception('Cannot recognise ccm: "{}"'.format(ccm))

        def f(*args):
            signal = normalised_g2_death_signal_heterotypic(*args)

            if len(signal) > 0:
                return signal[0]
            else:
                return 0

        status = sim_data.get_status()
        status_info = sim_data.get_status_info()
        t_events = sim_data.get_t_events()
        cell_count_A = sim_data.get_cell_count_for_clone(0)
        cell_count_B = sim_data.get_cell_count_for_clone(1)
        num_divisions_A = sim_data.get_num_divisions_for_clone(0)
        num_divisions_B = sim_data.get_num_divisions_for_clone(1)
        num_deaths_A = sim_data.get_num_deaths_for_clone(0)
        num_deaths_B = sim_data.get_num_deaths_for_clone(1)
        ergodic_rms_A = sim_data.get_ergodic_rms(f, coef_A, beta_A)
        ergodic_rms_B = sim_data.get_ergodic_rms(f, coef_B, beta_B)

        if status == 0:
            final_timestep = t_events[-1]
        else:
            final_timestep = t_events[-2]

        final_cell_count_A = cell_count_A[-1]
        final_cell_count_B = cell_count_B[-1]

        status_data.append(status)
        status_info_data.append(status_info)
        final_timestep_data.append(final_timestep)
        final_cell_count_A_data.append(final_cell_count_A)
        final_cell_count_B_data.append(final_cell_count_B)
        num_divisions_A_data.append(num_divisions_A)
        num_divisions_B_data.append(num_divisions_B)
        num_deaths_A_data.append(num_deaths_A)
        num_deaths_B_data.append(num_deaths_B)
        ergodic_rms_A_data.append(ergodic_rms_A)
        ergodic_rms_B_data.append(ergodic_rms_B)

    # Create and write dataframe
    df = pd.DataFrame({
        # Input parameters
        'beta_A' : beta_A_data[start_sim_num:end_sim_num],
        'tG_A' : tG_A_data[start_sim_num:end_sim_num],
        'rho_A' : rho_A_data[start_sim_num:end_sim_num],
        'eta_A' : eta_A_data[start_sim_num:end_sim_num],
        'coef_A' : coef_A_data[start_sim_num:end_sim_num],
        'beta_B' : beta_B_data[start_sim_num:end_sim_num],
        'tG_B' : tG_B_data[start_sim_num:end_sim_num],
        'rho_B' : rho_B_data[start_sim_num:end_sim_num],
        'eta_B' : eta_B_data[start_sim_num:end_sim_num],
        'coef_B' : coef_B_data[start_sim_num:end_sim_num],
        'n0_A' : n0_A_data[start_sim_num:end_sim_num],
        'n0_B' : n0_B_data[start_sim_num:end_sim_num],
        'seed' : seed_data[start_sim_num:end_sim_num],
        # Output data
        'status' : status_data,
        'status_info' : status_info_data,
        'final_timestep' : final_timestep_data,
        'final_cell_count_A' : final_cell_count_A_data,
        'final_cell_count_B' : final_cell_count_B_data,
        'num_divisions_A' : num_divisions_A_data,
        'num_divisions_B' : num_divisions_B_data,
        'num_deaths_A' : num_deaths_A_data,
        'num_deaths_B' : num_deaths_B_data,
        'ergodic_rms_A' : ergodic_rms_A_data,
        'ergodic_rms_B' : ergodic_rms_B_data,
        }, index=range(start_sim_num, end_sim_num))

    df.to_csv('{1}-data-{2:0{0}d}-to-{3:0{0}d}.csv'.format(len(str(len(beta_A_data))),
        datafile_prefix, start_sim_num, end_sim_num),
        index_label='simulation_id')

if __name__ == '__main__':

    ## 28 May 2021
    # Parameter sweep
    beta_As = np.array([0.3, 0.5, 0.7])
    tG_As = np.array([100])
    eta_As = np.array([0.05, 0.01, 0.2])
    coef_As = np.array([1])

    beta_Bs = np.array([0.3, 0.5, 0.7])
    tG_Bs = np.array([100])
    eta_Bs = np.array([0.05, 0.01, 0.2])
    coef_Bs = np.array([1])

    n0s = [(50, 50)]

    num_iter = 10

    # Simulation parameters
    tstart = 0
    tend = np.inf
    min_cell_count = 10
    max_cell_count = 1000

    ## End 28 May 2021

    ## 2 June 2021
    num_iter = 100
    ## End 2 June 2021

    ## 9 June 2021 - Fine parameter sweep

    beta_As = np.array([0.5])
    tG_As = np.array([100])
    eta_As = np.array([0.1])
    coef_As = np.array([1])

    beta_Bs = np.arange(0.05, 1, 0.05)
    tG_Bs = np.array([100])
    eta_Bs = np.arange(0.01, 0.26, 0.01)
    coef_Bs = np.array([1])

    n0s = [(50, 50)]

    num_iter = 50

    tend = 10000
    ## End 9 June 2021

    ## 16 June 2021 - Fine parameter sweep expanded
    beta_As = np.array([0.3, 0.5, 0.7])
    eta_As = np.array([0.05, 0.1, 0.2])
    ## End 16 June 2021

    ## 18 June 2021 - Fine parameter sweep with num_iter = 1
    beta_As = np.array([0.5])
    tG_As = np.array([100])
    eta_As = np.array([0.1])
    coef_As = np.array([1])

    beta_Bs = np.arange(0.05, 1, 0.05)
    tG_Bs = np.array([100])
    eta_Bs = np.arange(0.01, 0.26, 0.01)
    coef_Bs = np.array([1])

    n0s = [(50, 50)]

    num_iter = 1

    tend = 10000
    ## End 18 June 2021

    ## 18 June 2021 - Fine parameter sweep varying tG,A and tG,B
    tG_As = np.array([50, 100, 200])
    tG_Bs = np.array([50, 100, 200])

    num_iter = 50
    ## End 9 June 2021

    ## 22 June 2021 note ##
    ## New variable 'ccm' to choose the cell cycle model.
    ## End 22 June 2021 note ##

    ## 22 June 2021 - Fine parameter sweep expanded for uniform cell cycle model
    beta_As = np.array([0.3, 0.5, 0.7])
    tG_As = np.array([100])
    eta_As = np.array([0.05, 0.1, 0.2])
    coef_As = np.array([1])

    beta_Bs = np.arange(0.05, 1, 0.05)
    tG_Bs = np.array([100])
    eta_Bs = np.arange(0.01, 0.26, 0.01)
    coef_Bs = np.array([1])

    n0s = [(50, 50)]
    ccm = 'uniform'

    num_iter = 50
    tend = 10000

    datafile_prefix = 'heterotypic-model-2-uniform-ccm'

    ## End 22 June 2021

    if len(sys.argv) >= 2:
        seed = int(sys.argv[1])
    else:
        seed = None

    if len(sys.argv) == 4:
        sim_range = (int(sys.argv[2]), int(sys.argv[3]))
    else:
        sim_range = None

    heterotypic_model_2_parameter_sweep(beta_As, tG_As, eta_As, coef_As,
        beta_Bs, tG_Bs, eta_Bs, coef_Bs, n0s, datafile_prefix, ccm, num_iter,
        tstart, tend, min_cell_count, max_cell_count, seed, sim_range)
