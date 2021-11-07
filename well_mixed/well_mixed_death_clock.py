import numpy as np
import pandas as pd
import warnings
from scipy.stats import expon, uniform
from scipy.optimize import root_scalar

try:
    from scipy.integrate import solve_ivp
except ImportError:
    print("Warning, solve_ivp could not be imported.  Use f_is_stepwise_constant = True")
    def solve_ivp(*args, **kwargs):
        raise NotImplementedError

def interpolate_F_inverse(tau, t_grid, F_grid):

    if F_grid[-1] < tau:
        return np.inf

    ## Left bound
    # Largest i such that F_i <= tau
    i = np.max(np.nonzero(F_grid <= tau))

    # Smallest j such that F_j == F_i
    if tau == F_grid[i]:
        j = np.min(np.nonzero(F_grid == tau))
        return t_grid[j]

    ## Right bound
    k = i + 1

    return np.interp(tau, [F_grid[i], F_grid[k]], [t_grid[i], t_grid[k]])

def F_sister_bias(t, N, beta, tG):

    output = 0
    output += (t // tG) * (1 - beta) * tG

    remainder = np.mod(t, tG)

    if remainder < beta * tG:
        output += remainder * (N - 1) * (1 - beta) / N

    else:
        output += (N - 1) / N * (1 - beta) * beta * tG
        output += (remainder - beta * tG) * \
                (((N - 1) * (1 - beta) + 1) / N)

    return output

def inverse_F_sister_bias(tau, N, beta, tG):

    output = 0
    output += (tau // ((1 - beta) * tG)) * tG

    remainder = np.mod(tau, (1 - beta) * tG)

    if remainder < (N - 1) / N * (1 - beta) * beta * tG:
        output += remainder / ((N - 1) / N * (1 - beta))

    else:
        output += beta * tG
        output += (remainder - (N - 1) / N * (1 - beta) * beta * tG) / \
                (((N - 1) * (1 - beta) + 1) / N)

    return output

inverse_F_sister_bias = np.vectorize(inverse_F_sister_bias)

def average_f_g1_proportion(avg_tG1_fun, gamma_fun, tG, beta, Tdeath, c):

    ## Compute gamma (cell cycle model parameters)
    gamma = gamma_fun(beta, tG)

    # Define functions
    def eta_fun(f, tG, beta, Tdeath, c):
        if c * beta * tG * f == 0:
            output =  np.inf
        else:
            output =  Tdeath / (c * beta * tG * f)
        return output

    def g_fun(f, gamma, tG, beta, Tdeath, c):
        output = f - (1 - beta) * tG / \
        (avg_tG1_fun(eta_fun(f, tG, beta, Tdeath, c), gamma)
                    + (1 - beta) * tG )
        return output

    # Solution always lies between 0 and 1
    bracket = [0, 1]

    # Initialise output
    output = []

    sol = root_scalar(g_fun,
            args=(gamma, tG, beta, Tdeath, c),
            bracket=bracket,
            x0 = 1 - beta)

    if sol.converged:
        return sol.root
    else:
        return np.nan

average_f_g1_proportion = np.vectorize(average_f_g1_proportion)

def exponential_ccm(random_state=None, clone=None, tG1_param=50):

    return expon.rvs(scale=tG1_param, random_state=random_state)

def exponential_ccm_heterotypic(random_state=None, clone=None, tG1_param_clone_0=50, tG1_param_clone_1=50):

    assert clone == 0 or clone == 1

    if clone == 0:
        return expon.rvs(scale=tG1_param_clone_0, random_state=random_state)
    elif clone == 1:
        return expon.rvs(scale=tG1_param_clone_1, random_state=random_state)
    else:
        raise Exception('This should not be reached. Something has gone horribly wrong.')

def exponential_cdf(t, tG1_param=50):

    return 1 - np.exp(- t / tG1_param)

def uniform_ccm(random_state=None, clone=None, tG1_param=50, r=20):

    assert 0.5 * r / tG1_param <= 1

    return uniform.rvs(loc=tG1_param - 0.5 * r, scale=r, random_state=random_state)

def uniform_ccm_heterotypic(random_state=None, clone=None, tG1_param_clone_0=50,
        r_clone_0=20, tG1_param_clone_1=50, r_clone_1=20):

    assert 0.5 * r_clone_0 / tG1_param_clone_0 <= 1
    assert 0.5 * r_clone_1 / tG1_param_clone_1 <= 1
    assert clone == 0 or clone == 1

    if clone == 0:
        return uniform.rvs(loc=tG1_param_clone_0 - 0.5 * r_clone_0,
                scale=r_clone_0, random_state=random_state)
    elif clone == 1:
        return uniform.rvs(loc=tG1_param_clone_1 - 0.5 * r_clone_1,
                scale=r_clone_1, random_state=random_state)
    else:
        raise Exception('This should not be reached. Something has gone horribly wrong.')

def uniform_cdf(t, tG1_param=50, r=20):

    assert 0.5 * r / tG1_param <= 1

    if t < tG1_param - 0.5 * r:
        return 0

    if tG1_param + 0.5 * r < t:
        return 1

    return (t - (tG1_param - 0.5 * r)) / r

def base_rate_death_signal(t, tau, tbirth, tG1, clone, isinG1, base_rate=1):

    return base_rate * np.ones(tau.shape)

def base_rate_death_signal_heterotypic(t, tau, tbirth, tG1, clone, isinG1,
        base_rate_clone_0=1, base_rate_clone_1=1):

    assert np.all(np.logical_or(clone == 0, clone == 1))

    return base_rate_clone_0 * (clone == 0) + base_rate_clone_1 * (clone == 1)

def normalised_g2_death_signal(t, tau, tbirth, tG1, clone, isinG1, coef=1):

    # All cells in G1
    if np.all(isinG1):
        return np.zeros(tau.shape)

    # All cells in G2
    if np.all(np.logical_not(isinG1)):
        return coef * np.ones(tau.shape)

    # Neither of these scenarios
    return coef * np.sum(np.logical_not(isinG1)) / (tau.size - 1) * np.ones(tau.shape)

def normalised_g2_death_signal_heterotypic(t, tau, tbirth, tG1, clone, isinG1,
        coef_clone_0=1, coef_clone_1=1):

    assert np.all(np.logical_or(clone == 0, clone == 1))

    coef = coef_clone_0 * (clone == 0) + coef_clone_1 * (clone == 1)

    # All cells in G1
    if np.all(isinG1):
        return np.zeros(tau.shape)

    # All cells in G2
    if np.all(np.logical_not(isinG1)):
        return coef * np.ones(tau.shape)

    # Neither of these scenarios
    return coef * np.sum(np.logical_not(isinG1)) / (tau.size - 1) * np.ones(tau.shape)

def g2_death_signal(t, tau, tbirth, tG1, clone, isinG1, coef=1):

    return coef * np.sum(np.logical_not(isinG1)) * np.ones(tau.shape)

def g2_death_signal_heterotypic(t, tau, tbirth, tG1, clone, isinG1,
        coef_clone_0=1, coef_clone_1=1):

    assert np.all(np.logical_or(clone == 0, clone == 1))

    coef = coef_clone_0 * (clone == 0) + coef_clone_1 * (clone == 1)

    return coef * np.sum(np.logical_not(isinG1)) * np.ones(tau.shape)

class WellMixedSimulator(object):

    def __init__(self, f=base_rate_death_signal, ccm=exponential_ccm,
            Tdeath=100, tG2=50, tstart=0, tend=500,
            f_args=(),
            ccm_args=(),
            max_cell_count=np.inf,
            min_cell_count=0,
            f_is_stepwise_constant=True,
            min_cell_count_for_clone={},
            max_cell_count_for_clone={},
            apoptosis_at_checkpoint=False,
            switch_apoptosis_time=None,
            ):

        # Some assertions
        assert callable(f)
        assert callable(ccm)
        if not callable(Tdeath):
            assert Tdeath >= 0
            Tdeath = lambda clone, Tdeath=Tdeath: Tdeath * np.ones(clone.shape)
        if not callable(tG2):
            assert tG2 >= 0
            tG2 = lambda clone, tG2=tG2: tG2 * np.ones(clone.shape)
        assert tend >= tstart

        self.f = f
        self.ccm = ccm
        self.Tdeath = Tdeath
        self.tG2 = tG2
        self.tstart = tstart
        self.tend = tend
        self.f_args = f_args
        self.ccm_args = ccm_args
        self.min_cell_count = min_cell_count
        self.max_cell_count = max_cell_count
        self.f_is_stepwise_constant = f_is_stepwise_constant
        self.max_cell_count_for_clone = max_cell_count_for_clone
        self.min_cell_count_for_clone = min_cell_count_for_clone
        self.apoptosis_at_checkpoint = apoptosis_at_checkpoint
        self.switch_apoptosis_time = switch_apoptosis_time
        self.switched = False

        # Save parameters in dict for output
        self.param = {
            'f' : f,
            'ccm' : ccm,
            'Tdeath' : Tdeath,
            'tG2' : tG2,
            'tstart' : tstart,
            'tend' : tend,
            'f_args' : f_args,
            'ccm_args' : ccm_args,
            'min_cell_count' : min_cell_count,
            'max_cell_count' : max_cell_count,
            'f_is_stepwise_constant' : f_is_stepwise_constant,
            'min_cell_count_for_clone' : min_cell_count_for_clone,
            'max_cell_count_for_clone' : max_cell_count_for_clone,
            'apoptosis_at_checkpoint' : apoptosis_at_checkpoint,
            'switch_apoptosis_time' : switch_apoptosis_time,
            }

        # Define division, transition, and death events
        if not self.f_is_stepwise_constant:

            if self.apoptosis_at_checkpoint:
                raise NotImplementedError

            if self.switch_apoptosis_time:
                raise NotImplementedError

            division_event = lambda t, tau, tbirth, tG1, clone, isinG1, *args: \
                    np.amax(t - tbirth - tG1 - self.tG2(clone))
            division_event.terminal = True
            division_event.direction = 1

            def transition_event(t, tau, tbirth, tG1, clone, isinG1, *args):
                if np.any(isinG1):
                    return np.amax(t - tbirth[isinG1] - tG1[isinG1])
                else:
                    return -1
            transition_event.terminal = True
            transition_event.direction = 1

            def death_event(t, tau, tbirth, tG1, clone, isinG1, *args):
                if np.any(isinG1):
                    return np.amax(tau[isinG1] - self.Tdeath(clone)[isinG1])
                else:
                    return -1
            death_event.terminal = True
            death_event.direction = 1

            self.events = [ division_event, transition_event, death_event ]

    def run(self,
            tau_0=np.zeros(4),
            tbirth_0=np.random.rand(4) * -100,
            tG1_0=expon.rvs(scale=50, size=4) + 50,
            clone_0=np.arange(4),
            seed=None,
            ):
        """
        Runs Monte Carlo Death Clock simulation with given initial conditions
        and returns dictionary with following items:

        (t1, t2, ..., tN are the time points of the N division/transition/death events.)

        t_events: [ tstart, t1, t2, ..., tN, tend ]

        t_grid: [ np.array(time points for ODE solution) per time interval ]
        tau: [ num-cells-by-num-timepoints np.array(tau) per time interval ]
        cell_indices: [ np.array(indices of cells) per time interval ]
        isinG1: [ np.array(True if in G1, else False) per time interval ]

        division: [ np.array(indices of cells undergoing division) per event ]
        transition: [ np.array(indices of cells undergoing transition) per event ]
        death: [ np.array(indices of cells undergoing death) per event ]

        tbirth: np.array(birth times per cell)
        tG1: np.array(G1 durations per cell)
        clone: np.array(clones per cell)

        status: 0 means end of simulation time reached
                1 means zero cell count reached
                2 means max cell count reached
                3 means min cell count reached
                4 means clone-specific max cell count reached (see status_info for clone)
                5 means clone-specific min cell count reached (see status_info for clone)

        status_info: status = 0: None
                     status = 1: None
                     status = 2: None
                     status = 3: None
                     status = 4: clone whose max cell count was reached
                     status = 5: clone whose min cell count was reached

        param: dictionary containing parameters
               f: death clock signal function
               ccm: cell cycle model function
               Tdeath: death threshold function
               tG2: G2 duration function
               tstart: start time
               tend: end time
               f_args: additional args to f
               ccm_args: additional args to ccm
               min_cell_count: self-explanatory
               max_cell_count: self-explanatory
               f_is_stepwise_constant: whether f is stepwise constant between events
               min_cell_count_for_clone: dict containing clone-specific minimum cell count
               max_cell_count_for_clone: dict containing clone-specific maximum cell count

        init_cond: dictionary containing initial conditions
                   tau_0: tau
                   tbirth_0: birth times
                   tG1_0: G1 durations
                   clone_0: clones
                   seed: seed for random number generator
        """

        self.tau_0 = np.array(tau_0, dtype=float)
        self.tbirth_0 = np.array(tbirth_0, dtype=float)
        self.tG1_0 = np.array(tG1_0, dtype=float)
        self.clone_0 = np.array(clone_0, dtype=int)

        self.check_initial_conditions()

        # Create random state if seed is not None
        if not seed is None:
            self.random_state = np.random.RandomState(seed)
        else:
            self.random_state = None

        # Save initial conditions in dict for output
        self.init_cond = {
            'tau_0' : np.array(tau_0, dtype=float),
            'tbirth_0' : np.array(tbirth_0, dtype=float),
            'tG1_0' : np.array(tG1_0, dtype=float),
            'clone_0' : np.array(clone_0, dtype=int),
            'seed' : seed,
            }

        # Index cells
        self.cell_indices_now = np.arange(len(self.tbirth_0))
        self.last_cell_index = self.cell_indices_now[-1]

        # Initialise simulation time
        self.t_now = self.tstart

        # Initialise state variables
        self.tau_now = np.array(self.tau_0)
        self.tbirth_now = np.array(self.tbirth_0)
        self.tG1_now = np.array(self.tG1_0)
        self.clone_now = np.array(self.clone_0)

        # Helper state variable
        self.isinG1_now = self.t_now - self.tbirth_now < self.tG1_now

        # Initialise output data
        self.t_events_data = [self.t_now]

        self.t_grid_data = []
        self.tau_data = []

        self.division_data = []
        self.transition_data = []
        self.death_data = []

        self.cell_indices_data = []
        self.isinG1_data = []

        self.tbirth_data = np.array(self.tbirth_now)
        self.tG1_data = np.array(self.tG1_now)
        self.clone_data = np.array(self.clone_now)

        # Initialise status and status_info
        self.status = None
        self.status_info = None

        # Simulation loop
        while True:


            t, tau, event_occurred = self.solve_until_next_event()


            # Switch apoptosis mode if required
            if self.switch_apoptosis_time:
                if t[-1] >= self.switch_apoptosis_time and not self.switched:
                    self.apoptosis_at_checkpoint = not self.apoptosis_at_checkpoint
                    self.switched = True

            # Save data
            self.t_events_data.append(t[-1])
            self.t_grid_data.append(t)
            self.tau_data.append(tau)
            self.cell_indices_data.append(np.array(self.cell_indices_now))
            self.isinG1_data.append(np.array(self.isinG1_now))

            # If division, transition or death event occurred
            if event_occurred:

                # Update tau and t
                self.tau_now = np.array(self.tau_data[-1][:,-1])
                self.t_now = self.t_events_data[-1]

                # If event happened or past end simulation time, break loop and
                # do NOT record event.  In other words, discrete events
                # happening at exactly tend will not be recorded.
                if self.t_now >= self.tend:
                    break


                self.do_cell_transitions_divisions_deaths()


                # If there are no cells remaining, tack on last data item and
                # break
                if len(self.tau_now) == 0:

                    self.t_events_data.append(self.tend)
                    self.t_grid_data.append(np.array([self.t_now, self.tend]))
                    self.tau_data.append(np.array([]))
                    self.cell_indices_data.append(np.array([]))
                    self.isinG1_data.append(np.array([]))

                    self.status = 1
                    break

                # If the maximum number of cells is hit, tack on last data item
                # and break
                if len(self.tau_now) >= self.max_cell_count:

                    self.t_events_data.append(self.tend)
                    self.t_grid_data.append(np.array([self.t_now, self.tend]))
                    self.tau_data.append(np.array([self.tau_now, self.tau_now]).T)
                    self.cell_indices_data.append(np.array(self.cell_indices_now))
                    self.isinG1_data.append(np.array(self.isinG1_now))

                    self.status = 2
                    break

                # If the minimum number of cells is hit, tack on last data item
                # and break
                if len(self.tau_now) <= self.min_cell_count:

                    self.t_events_data.append(self.tend)
                    self.t_grid_data.append(np.array([self.t_now, self.tend]))
                    self.tau_data.append(np.array([self.tau_now, self.tau_now]).T)
                    self.cell_indices_data.append(np.array(self.cell_indices_now))
                    self.isinG1_data.append(np.array(self.isinG1_now))

                    self.status = 3
                    break

                # If the maximum number of cells for a clone is hit, tack on
                # last data item and break
                for clone, max_cell_count in self.max_cell_count_for_clone.items():
                    if np.sum(self.clone_now == clone) >= max_cell_count:

                        self.t_events_data.append(self.tend)
                        self.t_grid_data.append(np.array([self.t_now, self.tend]))
                        self.tau_data.append(np.array([self.tau_now, self.tau_now]).T)
                        self.cell_indices_data.append(np.array(self.cell_indices_now))
                        self.isinG1_data.append(np.array(self.isinG1_now))

                        self.status = 4
                        self.status_info = clone
                        break

                # Break out of while loop
                if self.status == 4:
                    break

                # If the minimum number of cells for a clone is hit, tack on
                # last data item and break
                for clone, min_cell_count in self.min_cell_count_for_clone.items():
                    if np.sum(self.clone_now == clone) <= min_cell_count:

                        self.t_events_data.append(self.tend)
                        self.t_grid_data.append(np.array([self.t_now, self.tend]))
                        self.tau_data.append(np.array([self.tau_now, self.tau_now]).T)
                        self.cell_indices_data.append(np.array(self.cell_indices_now))
                        self.isinG1_data.append(np.array(self.isinG1_now))

                        self.status = 5
                        self.status_info = clone
                        break

                # Break out of while loop
                if self.status == 5:
                    break

            # Else simulation has terminated
            else:
                self.status = 0
                break

        return {
                't_events' : self.t_events_data,
                't_grid' : self.t_grid_data,
                'tau' : self.tau_data,
                'cell_indices' : self.cell_indices_data,
                'isinG1' : self.isinG1_data,
                'division' : self.division_data,
                'transition' : self.transition_data,
                'death' : self.death_data,
                'tbirth' : self.tbirth_data,
                'tG1' : self.tG1_data,
                'clone' : self.clone_data,
                'status' : self.status,
                'status_info' : self.status_info,
                'param' : self.param,
                'init_cond' : self.init_cond,
                }

    def sample_g1_duration(self, clone):
        return self.ccm(self.random_state, clone, *self.ccm_args)

    def check_initial_conditions(self):

        if not len(self.tau_0) == len(self.tbirth_0) == len(self.tG1_0) == len(self.clone_0):
            raise ValueError("tau_0, tbirth_0, tG1_0, and clone_0 must have the same length")

        if not len(self.tau_0) < self.max_cell_count:
            raise ValueError("The initial cell count ({}) must be smaller than "
            "the maximum cell count ({})".format(len(self.tau_0), self.max_cell_count))

        if not len(self.tau_0) > 0:
            raise ValueError("The initial cell count ({}) must be larger than "
            "0".format(len(self.tau_0)))

        if not len(self.tau_0) > self.min_cell_count:
            raise ValueError("The initial cell count ({}) must be larger than "
            "the minimum cell count ({})".format(len(self.tau_0), self.min_cell_count))

        for clone, min_cell_count in self.min_cell_count_for_clone.items():
            if not np.sum(self.clone_0 == clone) > min_cell_count:
                raise ValueError("The initial cell count of clone {0} ({1}) must "
                "be larger than the minimum cell count for clone {0} ({2})".format(
                    clone, np.sum(self.clone_0 == clone), min_cell_count))

        for clone, max_cell_count in self.max_cell_count_for_clone.items():
            if not np.sum(self.clone_0 == clone) < max_cell_count:
                raise ValueError("The initial cell count of clone {0} ({1}) must "
                "be smaller than the maximum cell count for clone {0} ({2})".format(
                    clone, np.sum(self.clone_0 == clone), max_cell_count))

        if not np.all(self.tbirth_0 <= self.tstart):
            raise ValueError("Initial birth times must be at or before start of simulation")

        if not self.apoptosis_at_checkpoint:
            if not np.all(np.logical_or(self.tau_0 < self.Tdeath(self.clone_0),
                np.logical_and(self.tau_0 >= self.Tdeath(self.clone_0), self.tstart
                - self.tbirth_0 >= self.tG1_0))):
                raise ValueError("Death invariant is violated in initial conditions")

        if not np.all(self.tstart - self.tbirth_0 < self.tG1_0 + self.tG2(self.clone_0)):
            raise ValueError("Birth invariant is violated in initial conditions")

    def do_cell_divisions(self):

        division_indices = np.nonzero(np.isclose(self.t_now - self.tbirth_now, self.tG1_now + self.tG2(self.clone_now)))

        self.division_data.append(self.cell_indices_now[division_indices])

        for division_index in division_indices[0]:

            # Reset death clock of daughter cells
            self.tau_now[division_index] = 0
            self.tau_now = np.append(self.tau_now, 0)

            # Set time of birth on daughter cells
            self.tbirth_now[division_index] = self.t_now
            self.tbirth_now = np.append(self.tbirth_now, self.t_now)

            # Draw random G1 duration for daughter cells
            self.tG1_now[division_index] = self.sample_g1_duration(self.clone_now[division_index])
            self.tG1_now = np.append(self.tG1_now, self.sample_g1_duration(self.clone_now[division_index]))

            # Set clone of new daughter cell
            self.clone_now = np.append(self.clone_now, self.clone_now[division_index])

            # Both cell starts in G1
            self.isinG1_now[division_index] = True
            self.isinG1_now = np.append(self.isinG1_now, True)

            # Generate new indices for cells
            self.cell_indices_now[division_index] = self.last_cell_index + 1
            self.cell_indices_now = np.append(self.cell_indices_now, self.last_cell_index + 2)

            # Save static data
            self.tbirth_data = np.append(self.tbirth_data, [self.t_now, self.t_now])
            self.tG1_data = np.append(self.tG1_data,
                    [ self.tG1_now[division_index], self.tG1_now[-1] ])
            self.clone_data = np.append(self.clone_data,
                    [ self.clone_now[-1], self.clone_now[-1] ])

            # Update last cell index
            self.last_cell_index += 2

    def do_cell_transitions(self):

        if self.apoptosis_at_checkpoint:
            transition_indices = np.nonzero(
                    np.logical_and(np.isclose(self.t_now - self.tbirth_now, self.tG1_now),
                        self.tau_now < self.Tdeath(self.clone_now)))
        else:
            transition_indices = np.nonzero(np.isclose(self.t_now - self.tbirth_now, self.tG1_now))

        self.transition_data.append(self.cell_indices_now[transition_indices])

        self.isinG1_now[transition_indices] = False

    def do_cell_deaths(self):

        if self.apoptosis_at_checkpoint:
            death_indices = np.nonzero(
                    np.logical_and(np.isclose(self.t_now - self.tbirth_now, self.tG1_now),
                        self.tau_now >= self.Tdeath(self.clone_now)))
        else:
            death_indices = np.nonzero(
                    np.logical_and(
                        np.isclose(self.tau_now, self.Tdeath(self.clone_now)),
                        self.isinG1_now))

        self.death_data.append(self.cell_indices_now[death_indices])

        # Traverse from last index to first, else indexing is incorrect
        for death_index in death_indices[0][::-1]:

            # Remove dead cell
            self.tau_now = np.delete(self.tau_now, death_index)
            self.tbirth_now = np.delete(self.tbirth_now, death_index)
            self.tG1_now = np.delete(self.tG1_now, death_index)
            self.clone_now = np.delete(self.clone_now, death_index)
            self.isinG1_now = np.delete(self.isinG1_now, death_index)
            self.cell_indices_now = np.delete(self.cell_indices_now, death_index)

    def do_cell_transitions_divisions_deaths(self):

        self.do_cell_transitions()
        self.do_cell_divisions()
        self.do_cell_deaths()

    def solve_until_next_event_for_constant_f(self):

        # Compute constant f
        f = self.f(self.t_now, self.tau_now, self.tbirth_now, self.tG1_now,
                self.clone_now, self.isinG1_now, *self.f_args)

        # Find time until next event
        time_until_next_division = np.min(self.tbirth_now + self.tG1_now + self.tG2(self.clone_now) - self.t_now)

        if np.any(self.isinG1_now):

            tbirth_filtered = self.tbirth_now[self.isinG1_now]
            tG1_filtered = self.tG1_now[self.isinG1_now]
            tau_filtered = self.tau_now[self.isinG1_now]
            f_filtered = f[self.isinG1_now]
            Tdeath_filtered = self.Tdeath(self.clone_now)[self.isinG1_now]

            time_until_next_transition = np.min(tbirth_filtered + tG1_filtered - self.t_now)
            with np.errstate(divide='ignore'):
                time_until_next_death = np.min((Tdeath_filtered - tau_filtered) / f_filtered)

        else:
            time_until_next_transition = np.inf
            time_until_next_death = np.inf

        # If apoptosis at checkpoint, then death events are at the same time as
        # transition events
        if self.apoptosis_at_checkpoint:
            time_until_next_death = np.inf

        time_until_next_event = \
            np.min([time_until_next_division, time_until_next_transition, time_until_next_death])

        assert time_until_next_event >= 0

        # Event occurred if the time of the next event is before termination
        event_occurred = self.t_now + time_until_next_event < self.tend

        # Get timestep
        if event_occurred:
            timestep = time_until_next_event
        else:
            timestep = self.tend - self.t_now

        # Compute and return results
        output_t = np.array([self.t_now, self.t_now + timestep])
        output_tau = np.array([self.tau_now, self.tau_now + f * timestep]).T

        return output_t, output_tau, event_occurred

    def solve_until_next_event_for_nonconstant_f(self):

        # Solve ODE
        output = solve_ivp(self.f, [self.t_now, self.tend], self.tau_now,
                args=(self.tbirth_now, self.tG1_now, self.clone_now, self.isinG1_now) + self.f_args,
                events=self.events)

        if not output.success:
            raise Exception('An error occured in scipy.integrate.solve_ivp: "{}"'.format(
                output.message))

        if output.status == 1:
            # Sanity check: only one event should have happened
            assert np.sum([len(events) for events in output.y_events]) == 1

        return output.t, output.y, output.status == 1

    def solve_until_next_event(self):

        if self.f_is_stepwise_constant:
            return self.solve_until_next_event_for_constant_f()

        else:
            return self.solve_until_next_event_for_nonconstant_f()

class WellMixedSimulationData(object):

    def __init__(self, data):

        # Store raw data
        self.data = data

        # Initialise all data members to None for lazy initialisation
        # global data
        self.status = None
        self.status_info = None
        self.unique_clones = None
        self.num_divisions = None
        self.num_transitions = None
        self.num_deaths = None
        self.total_cell_count = None # sum of initial cells and cells born during simulation
        self.num_divisions_for_clone = {}
        self.num_transitions_for_clone = {}
        self.num_deaths_for_clone = {}
        self.f = None

        # timeseries data (per time interval)
        self.t_events = None
        self.cell_count = None
        self.clone_cell_count = None
        self.G1_cell_count = None
        self.G2_cell_count = None
        self.timeseries_df = None
        self.G1_cell_count_for_clone = {}
        self.G2_cell_count_for_clone = {}
        self.cell_count_for_clone = {}

        # fine timeseries data (per ODE time point)
        self.t_grid = None
        self.tau = {}

        # cellwise data
        self.tbirth = None                  # Time of birth
        self.tG1 = None                     # G1 duration
        self.clone = None                   # clone
        self.died = None                    # whether cell died during simulation
        self.tdeath = None                  # time of death, inf if not died
        self.divided = None                 # whether cell divided during simulation
        self.tdivision = None               # time of division, inf if not divided
        self.transitioned = None            # whether cell transitioned to G2
        self.ttransition = None             # time of transition to G2, inf if not transitioned
        self.t_last_alive = None            # time that cell was last alive
                                            # (until it divided, died or simulation terminated)
        self.max_age = None                 # t_last_alive - tbirth
        self.time_in_G1 = None              # time spent in G1 (from birth until
                                            # transition/death/termination)
        self.time_in_G2 = None              # time spent in G2 (from transition until
                                            # termination/divions, 0 if not transitioned)
        self.last_tau = None                # Last death clock value
        self.average_f = None               # Average death clock signal: last death
                                            # clock value divided by time spent in G1
        self.effective_time_in_G1 = None    # Effective time spent in G1 (from birth until
                                            # transition/death/termination) until first
                                            # cell that did neither die nor transition
        self.cellwise_df = None

    def __str__(self):

        if self.get_status() == 0:
            status_str = 'End of simulation reached'

        elif self.get_status() == 1:
            status_str = 'Extinction reached'

        elif self.get_status() == 2:
            status_str = 'Maximum cell count reached'

        elif self.get_status() == 3:
            status_str = 'Minimum cell count reached'

        else:
            raise Exception('Never reached')

        unique_clones = self.get_unique_clones()
        num_divisions = self.get_num_divisions()
        num_transitions = self.get_num_transitions()
        num_deaths = self.get_num_deaths()
        total_cell_count = self.get_total_cell_count()
        timeseries_df = self.get_timeseries_df()
        timeseries_df_all = timeseries_df.describe(include='all')
        cellwise_df = self.get_cellwise_df()
        with warnings.catch_warnings():
            # cellwise_df may contain 'inf' values (e.g. for tdeath) This
            # creates warnings when computing summary statistics, which we
            # suppressed by this construction
            warnings.simplefilter('ignore')
            cellwise_df_all = cellwise_df.describe(include='all')

        return """
Global data
-----------

status:\t{}
unique_clones:\t\t{}
num_divisions:\t\t{}
num_transitions:\t{}
num_deaths:\t\t{}
total_cell_count:\t{}

Timeseries data
---------------

{}

{}

Cellwise data
-------------

{}

{}
""".format(status_str,
unique_clones,
num_divisions,
num_transitions,
num_deaths,
total_cell_count,
timeseries_df,
timeseries_df_all,
cellwise_df,
cellwise_df_all
)

    def get_unique_clones(self):

        if self.unique_clones is None:
            self.unique_clones = np.unique(self.data['clone'])

        return self.unique_clones

    def get_num_divisions(self):

        if self.num_divisions is None:
            self.num_divisions = np.sum(self.get_divided())

        return self.num_divisions

    def get_num_divisions_for_clone(self, clone):

        if not clone in self.num_divisions_for_clone:
            self.num_divisions_for_clone[clone] = np.sum(
                    np.logical_and(self.get_clone() == clone, self.get_divided()))

        return self.num_divisions_for_clone[clone]

    def get_num_transitions(self):
        # Note that this is not the sum(self.get_transitioned()), because that
        # function also counts cells that are in G2 at the start of the
        # simulation as transitioned

        if self.num_transitions is None:
            self.num_transitions = sum(array.size for array in self.data['transition'])

        return self.num_transitions

    def get_num_transitions_for_clone(self, clone):

        if not clone in self.num_transitions_for_clone:
            self.num_transitions_for_clone[clone] = sum(
                    np.sum(self.get_clone()[array] == clone) for array in self.data['transition'])

        return self.num_transitions_for_clone[clone]

    def get_num_deaths(self):

        if self.num_deaths is None:
            self.num_deaths = np.sum(self.get_died())

        return self.num_deaths

    def get_num_deaths_for_clone(self, clone):

        if not clone in self.num_deaths_for_clone:
            self.num_deaths_for_clone[clone] = np.sum(
                    np.logical_and(self.get_clone() == clone, self.get_died()))

        return self.num_deaths_for_clone[clone]

    def get_status(self):

        if self.status is None:
            self.status = self.data['status']

        return self.status

    def get_status_info(self):

        if self.status_info is None:
            self.status_info = self.data['status_info']

        return self.status_info

    def get_t_events(self):

        if self.t_events is None:
            self.t_events = np.array(self.data['t_events'])

        return self.t_events

    def get_cell_count(self):

        if self.cell_count is None:
            self.cell_count = np.array(
                [len(cell_indices) for cell_indices in self.data['cell_indices']])

            # Repeat last cell count at termination time
            self.cell_count = np.append(self.cell_count, self.cell_count[-1])

        return self.cell_count

    def get_cell_count_for_clone(self, clone):

        if not clone in self.cell_count_for_clone:
            cell_count = []
            for cell_indices in self.data['cell_indices']:

                cell_count.append(np.sum(self.get_clone()[cell_indices] == clone))

            # Repeat last cell count at termination time
            cell_count.append(cell_count[-1])

            # Save as numpy array
            self.cell_count_for_clone[clone] = np.array(cell_count)

        return self.cell_count_for_clone[clone]

    def get_total_cell_count(self):

        if self.total_cell_count is None:
            self.total_cell_count = len(self.get_tbirth())

        return self.total_cell_count

    def get_clone_cell_count(self):

        if self.clone_cell_count is None:

            unique_clones = self.get_unique_clones()
            self.clone_cell_count = np.zeros((len(self.get_t_events()), len(unique_clones)))

            for i, cell_indices in enumerate(self.data['cell_indices']):

                if len(cell_indices) == 0:
                    continue

                clones = self.get_clone()[cell_indices]

                unique_clones_now, counts = np.unique(clones, return_counts=True)
                for clone, count in zip(unique_clones_now, counts):

                    idx = np.searchsorted(unique_clones, clone)
                    self.clone_cell_count[i, idx] = count

            # Repeat last cell count at termination time
            self.clone_cell_count[-1,:] = self.clone_cell_count[-2, :]

        return self.clone_cell_count

    def get_G1_cell_count(self):

        if self.G1_cell_count is None:

            self.G1_cell_count = np.array([np.sum(isinG1)
                for isinG1 in self.data['isinG1']])

            # Repeat last cell count at termination time
            self.G1_cell_count = np.append(self.G1_cell_count, self.G1_cell_count[-1])

        return self.G1_cell_count

    def get_G1_cell_count_for_clone(self, clone):

        if not clone in self.G1_cell_count_for_clone:
            G1_cell_count = []
            for cell_indices, isinG1 in zip(self.data['cell_indices'], self.data['isinG1']):

                G1_cell_count.append(np.sum(np.logical_and(self.get_clone()[cell_indices] == clone,
                    isinG1)))

            # Repeat last cell count at termination time
            G1_cell_count.append(G1_cell_count[-1])

            # Save as numpy array
            self.G1_cell_count_for_clone[clone] = np.array(G1_cell_count)

        return self.G1_cell_count_for_clone[clone]

    def get_G2_cell_count(self):

        if self.G2_cell_count is None:

            self.G2_cell_count = np.array([np.sum(np.logical_not(isinG1))
                for isinG1 in self.data['isinG1']])

            # Repeat last cell count at termination time
            self.G2_cell_count = np.append(self.G2_cell_count, self.G2_cell_count[-1])

        return self.G2_cell_count

    def get_G2_cell_count_for_clone(self, clone):

        if not clone in self.G2_cell_count_for_clone:
            G2_cell_count = []
            for cell_indices, isinG1 in zip(self.data['cell_indices'], self.data['isinG1']):

                G2_cell_count.append(np.sum(np.logical_and(self.get_clone()[cell_indices] == clone,
                    np.logical_not(isinG1))))

            # Repeat last cell count at termination time
            G2_cell_count.append(G2_cell_count[-1])

            # Save as numpy array
            self.G2_cell_count_for_clone[clone] = np.array(G2_cell_count)

        return self.G2_cell_count_for_clone[clone]

    def get_timeseries_df(self):

        if self.timeseries_df is None:

            timeseries_dict = {
                't_events' : self.get_t_events(),
                'cell_count' : self.get_cell_count(),
                }

            for clone, cell_count in \
                    zip(self.get_unique_clones(), self.get_clone_cell_count().T):
                timeseries_dict['clone_{}'.format(clone)] = cell_count

            self.timeseries_df = pd.DataFrame(timeseries_dict)

        return self.timeseries_df

    def get_t_grid(self):

        if self.t_grid is None:

            self.t_grid = np.concatenate(self.data['t_grid'])

        return self.t_grid

    def get_tau_for_cell_index(self, cell_index):

        if not cell_index in self.tau:

            cur_tau = np.array([], dtype=float)

            last_G1_tau = 0

            # Loop over time intervals
            for t_grid, tau, cell_indices, isinG1 in \
                    zip(self.data['t_grid'], self.data['tau'],
                            self.data['cell_indices'], self.data['isinG1']):

                # Cell is alive in this time interval
                if cell_index in cell_indices:

                    matches = np.nonzero(cell_indices == cell_index)[0]
                    assert len(matches) == 1
                    local_idx = matches[0]

                    if isinG1[local_idx]:

                        cur_tau = np.append(cur_tau, tau[local_idx])
                        last_G1_tau = cur_tau[-1]

                    else:

                        cur_tau = np.append(cur_tau, t_grid * 0 + last_G1_tau)

                else:

                    cur_tau = np.append(cur_tau, t_grid * 0)

            self.tau[cell_index] = cur_tau

        return self.tau[cell_index]

    def get_last_tau(self):

        if self.last_tau is None:

            self.last_tau = np.zeros(self.get_total_cell_count(), dtype=float)

            # Loop over time intervals
            for t_grid, tau, cell_indices, isinG1 in \
                    zip(self.data['t_grid'], self.data['tau'],
                            self.data['cell_indices'], self.data['isinG1']):

                cur_last_taus = tau[:, -1]

                self.last_tau[cell_indices[isinG1]] = cur_last_taus[isinG1]

        return self.last_tau

    def get_average_f(self):

        if self.average_f is None:

            # If time_in_G1 is zero, then replace by inf so that average f
            # computes to zero.
            t = np.array(self.get_time_in_G1())
            t[t == 0] = np.inf

            self.average_f = self.get_last_tau() / t

        return self.average_f

    def get_tbirth(self):

        if self.tbirth is None:
            self.tbirth = self.data['tbirth']

        return self.tbirth

    def get_tG1(self):

        if self.tG1 is None:
            self.tG1 = self.data['tG1']

        return self.tG1

    def get_clone(self):

        if self.clone is None:
            self.clone = self.data['clone']

        return self.clone

    def get_died(self):

        if self.died is None:
            self.died = np.zeros(self.get_total_cell_count(), dtype=bool)
            for death_indices in self.data['death']:
                self.died[death_indices] = True

        return self.died

    def get_tdeath(self):

        if self.tdeath is None:
            self.tdeath = np.ones(self.get_total_cell_count(), dtype=float) * np.inf
            for time, death_indices in zip(self.get_t_events()[1:-1], self.data['death']):
                self.tdeath[death_indices] = time

        return self.tdeath

    def get_divided(self):

        if self.divided is None:
            self.divided = np.zeros(self.get_total_cell_count(), dtype=bool)
            for division_indices in self.data['division']:
                self.divided[division_indices] = True

        return self.divided

    def get_tdivision(self):

        if self.tdivision is None:
            self.tdivision = np.ones(self.get_total_cell_count(), dtype=float) * np.inf
            for time, division_indices in zip(self.get_t_events()[1:-1], self.data['division']):
                self.tdivision[division_indices] = time

        return self.tdivision

    def get_transitioned(self):

        if self.transitioned is None:
            self.transitioned = np.zeros(self.get_total_cell_count(), dtype=bool)
            for transition_indices in self.data['transition']:
                self.transitioned[transition_indices] = True

            # Cells that were already in G2 at start of simulation should be
            # considered to have transitioned
            tstart = self.get_t_events()[0]
            cell_indices_start = self.data['cell_indices'][0]
            tG1_start = self.get_tG1()[cell_indices_start]
            tbirth_start = self.get_tbirth()[cell_indices_start]

            in_G2_at_start = tstart - tbirth_start >= tG1_start

            self.transitioned[cell_indices_start[in_G2_at_start]] = True

        return self.transitioned

    def get_ttransition(self):

        if self.ttransition is None:
            self.ttransition = np.ones(self.get_total_cell_count(), dtype=float) * np.inf
            for time, transition_indices in zip(self.get_t_events()[1:-1], self.data['transition']):
                self.ttransition[transition_indices] = time

            # Cells that were already in G2 at start of simulation should be
            # considered to have spent tG1 in G1
            tstart = self.get_t_events()[0]
            cell_indices_start = self.data['cell_indices'][0]
            tG1_start = self.get_tG1()[cell_indices_start]
            tbirth_start = self.get_tbirth()[cell_indices_start]

            in_G2_at_start = tstart - tbirth_start >= tG1_start

            self.ttransition[cell_indices_start[in_G2_at_start]] = \
                tbirth_start[in_G2_at_start] + tG1_start[in_G2_at_start]

        return self.ttransition

    def get_t_last_alive(self):

        if self.t_last_alive is None:

            # In case of normal termination, tend is last time point in
            # t_events, else the second to last one
            if self.get_status() == 0:
                tend = self.get_t_events()[-1]
            else:
                tend = self.get_t_events()[-2]

            self.t_last_alive = np.minimum(tend,
                    np.minimum(self.get_tdeath(), self.get_tdivision()))

        return self.t_last_alive

    def get_max_age(self):

        if self.max_age is None:
            self.max_age = self.get_t_last_alive() - self.get_tbirth()

        return self.max_age

    def get_time_in_G1(self):

        if self.time_in_G1 is None:
            self.time_in_G1 = np.minimum(self.get_ttransition(), self.get_t_last_alive()) \
                    - self.get_tbirth()

        return self.time_in_G1

    def get_time_in_G2(self):

        if self.time_in_G2 is None:

            self.time_in_G2 = np.zeros(self.get_total_cell_count(), dtype=float)

            transitioned = self.get_transitioned()
            t_last_alive = self.get_t_last_alive()
            ttransition = self.get_ttransition()

            self.time_in_G2[transitioned] = t_last_alive[transitioned] \
                    - ttransition[transitioned]

        return self.time_in_G2

    def get_effective_time_in_G1(self):

        if self.effective_time_in_G1 is None:
            died_or_transitioned = np.logical_or(self.get_died(), self.get_transitioned())
            indices = np.nonzero(np.logical_not(died_or_transitioned))[0]
            if len(indices) > 0:
                self.effective_time_in_G1 = self.get_time_in_G1()[:indices[0]]
            else:
                self.effective_time_in_G1 = self.get_time_in_G1()

        return self.effective_time_in_G1

    def get_cellwise_df(self):

        if self.cellwise_df is None:

            cellwise_dict = {
                'tbirth' : self.get_tbirth(),
                'tG1' : self.get_tG1(),
                'clone' : self.get_clone(),
                'died' : self.get_died(),
                'tdeath' : self.get_tdeath(),
                'divided' : self.get_divided(),
                'tdivision' : self.get_tdivision(),
                'transitioned' : self.get_transitioned(),
                'ttransition' : self.get_ttransition(),
                't_last_alive' : self.get_t_last_alive(),
                'max_age' : self.get_max_age(),
                'time_in_G1' : self.get_time_in_G1(),
                'time_in_G2' : self.get_time_in_G2(),
                }

            self.cellwise_df = pd.DataFrame(cellwise_dict)

        return self.cellwise_df

    def get_death_clock_signal(self, f, f_args=None):

        if f_args is None:
            f_args = self.data['param']['f_args']

        cur_f = []

        for t_grid, tau, cell_indices, isinG1 in \
                zip(self.data['t_grid'], self.data['tau'],
                        self.data['cell_indices'], self.data['isinG1']):

            tbirth = self.get_tbirth()[cell_indices]
            tG1 = self.get_tG1()[cell_indices]
            clone = self.get_clone()[cell_indices]

            for t, cur_tau in zip(t_grid, tau.T):

                cur_f.append(f(t, cur_tau, tbirth, tG1, clone, isinG1, *f_args))

        return np.array(cur_f)

    def get_integrated_death_clock_signal(self, f):

        t_grid = self.get_t_grid()
        cur_f = self.get_death_clock_signal(f)

        assert len(t_grid) == len(cur_f)

        # Apply trapezoid rule to integrate signal
        F = np.concatenate([[0.0], 0.5 * (cur_f[1:] + cur_f[:-1]) * (t_grid[1:] - t_grid[:-1])])
        F = np.cumsum(F)

        return F

    def get_ergodic_rms(self, f, c, beta):

        t_grid = self.get_t_grid()
        cur_f = self.get_death_clock_signal(f)

        square_error = (cur_f / c - (1 - beta))**2

        assert len(t_grid) == len(cur_f)

        # Apply trapezoid rule to integrate square error
        F = np.sum(0.5 * (square_error[1:] + square_error[:-1]) * (t_grid[1:] - t_grid[:-1]))

        rms = np.sqrt(F / (t_grid[-1] - t_grid[0]))

        return rms

    def get_f(self):

        if self.f is None:
            self.f = self.data['param']['f']

        return self.f

    def get_hindsight_survival_probability(self, t, f, cdf, ccm_args=None, Tdeath=None):

        # Make flat array of t
        t = np.array(t).flatten()

        # Check validity of t
        t_grid = self.get_t_grid()
        assert np.all(t_grid[0] <= t)
        assert np.all(t <= t_grid[-1])

        # Assert non-strict monotonicity of F
        F_grid = self.get_integrated_death_clock_signal(f)
        assert np.all(np.diff(F_grid) >= 0)

        # If t_grid[-1] is infinity, drop last element from t_grid and F_grid
        if t_grid[-1] == np.inf:
            t_grid = t_grid[:-1]
            F_grid = F_grid[:-1]

        # Retrieve parameters if needed
        if Tdeath is None:
            # None argument because we assume homotypic Tdeath
            Tdeath = self.data['param']['Tdeath'](np.array([None]))

        if ccm_args is None:
            ccm_args = self.data['param']['ccm_args']

        # Initialise output
        output = []

        for cur_t in t:

            ## Step 1: Find death time
            F0 = np.interp(cur_t, t_grid, F_grid)
            F1 = Tdeath + F0

            t1 = interpolate_F_inverse(F1, t_grid, F_grid)

            # Step 2: Compute time until death
            time_until_death = t1 - cur_t

            # Sanity check
            assert time_until_death >= 0

            # Step 3: if G1 duration is smaller than time until death, then the cell
            # survives.  Hence the survival probability is the cumulative distribution
            # function of the time until death
            output.append(cdf(time_until_death, *ccm_args))

        return np.array(output)

    def plot_cell_history(self, ax, modulo=np.inf):

        tbirth = self.get_tbirth()
        tG1 = self.get_tG1()
        transitioned = self.get_transitioned()
        ttransition = self.get_ttransition()
        t_last_alive = self.get_t_last_alive()

        first_G1 = True
        first_G2 = True

        for cell_idx, (tbirth, transitioned, ttransition, t_last_alive) in enumerate(zip(
            self.get_tbirth(),
            self.get_transitioned(),
            self.get_ttransition(),
            self.get_t_last_alive())):

            cell_idx = np.mod(cell_idx, modulo)

            if first_G1:
                first_G1 = False
                G1_label='G1'
            else:
                G1_label=None

            if transitioned:

                if first_G2:
                    first_G2 = False
                    G2_label='G2'
                else:
                    G2_label=None

                ax.plot([tbirth, ttransition], [cell_idx, cell_idx], 'b', label=G1_label)
                ax.plot([ttransition, t_last_alive], [cell_idx, cell_idx], 'r', label=G2_label)

            else:

                ax.plot([tbirth, t_last_alive], [cell_idx, cell_idx], 'b', label=G1_label)

        ax.set_xlabel('t')
        ax.set_ylabel('cell_idx mod {}'.format(modulo))
        ax.legend()

    def plot_cell_counts(self, ax):

        t_events = self.get_t_events()
        cell_count = self.get_cell_count()
        G1_cell_count = self.get_G1_cell_count()
        G2_cell_count = self.get_G2_cell_count()

        ax.step(t_events, cell_count, where='post', label='total')
        ax.step(t_events, G1_cell_count, where='post', label='G1')
        ax.step(t_events, G2_cell_count, where='post', label='G2')

        ax.set_xlabel('t')
        ax.set_ylabel('cell count')

        ax.legend()

if __name__ == "__main__":

    # Test 1: single cell with constant base rate and fixed G1, no divisions or deaths
    print("Test 1: no divisions or deaths")
    f = base_rate_death_signal
    ccm = lambda random_state, clone: 1

    Tdeath = 2
    tG2 = 1
    tstart = 0
    tend = 1

    tau_0 = [0]
    tbirth_0 = [0]
    tG1_0 = [1]
    clone_0 = [0]

    simulator = WellMixedSimulator(f, ccm, Tdeath, tG2, tstart, tend)
    data = simulator.run(tau_0, tbirth_0, tG1_0, clone_0)
    mc_data = WellMixedSimulationData(data)
    for key, value in data.items():
        print('{}: {}'.format(key, value))
    print()

    assert mc_data.get_num_divisions() == 0
    assert mc_data.get_num_transitions() == 0
    assert mc_data.get_num_deaths() == 0
    assert mc_data.get_status() == 0

    # Test 2: same, but with 1 division
    print("Test 2: 1 division")
    tbirth_0 = [-1.5]

    simulator = WellMixedSimulator(f, ccm, Tdeath, tG2, tstart, tend)
    data = simulator.run(tau_0, tbirth_0, tG1_0, clone_0)
    mc_data = WellMixedSimulationData(data)
    for key, value in data.items():
        print('{}: {}'.format(key, value))
    print()

    assert mc_data.get_num_divisions() == 1
    assert mc_data.get_num_transitions() == 0
    assert mc_data.get_num_deaths() == 0
    assert mc_data.get_status() == 0

    # Test 3: same, but with 3 divisions and 2 transitions, also test different
    # clone number
    print("Test 3: 3 divisions and 2 transitions")
    tstart, tend = 0, 3
    simulator = WellMixedSimulator(f, ccm, Tdeath, tG2, tstart, tend)
    data = simulator.run(tau_0, tbirth_0, tG1_0, [1])
    mc_data = WellMixedSimulationData(data)
    for key, value in data.items():
        print('{}: {}'.format(key, value))
    print()

    assert mc_data.get_num_divisions() == 3
    assert mc_data.get_num_transitions() == 2
    assert mc_data.get_num_deaths() == 0
    assert mc_data.get_status() == 0
    assert np.all(mc_data.get_tbirth() == [-1.5, 0.5, 0.5, 2.5, 2.5, 2.5, 2.5])
    assert np.all(mc_data.get_tG1() == [1.0] * 7)
    assert np.all(mc_data.get_clone() == [1] * 7)
    assert np.all(mc_data.get_died() == [False] * 7)
    assert np.all(mc_data.get_tdeath() == [np.inf] * 7)
    assert np.all(mc_data.get_divided() == [True, True, True, False, False, False, False])
    assert np.all(mc_data.get_tdivision() == [0.5, 2.5, 2.5, np.inf, np.inf, np.inf, np.inf])
    assert np.all(mc_data.get_transitioned() == [True, True, True, False, False, False, False])
    assert np.all(mc_data.get_ttransition() == [-0.5, 1.5, 1.5, np.inf, np.inf, np.inf, np.inf])
    assert np.all(mc_data.get_t_last_alive() == [0.5, 2.5, 2.5, 3, 3, 3, 3])
    assert np.all(mc_data.get_max_age() == [2, 2, 2, 0.5, 0.5, 0.5, 0.5])
    assert np.all(mc_data.get_time_in_G1() == [1, 1, 1, 0.5, 0.5, 0.5, 0.5,  ])
    assert np.all(mc_data.get_time_in_G2() == [1, 1, 1, 0.0, 0.0, 0.0, 0.0,  ])
    assert np.all(mc_data.get_effective_time_in_G1() == [1, 1, 1])

    # Test 4: same, but with 1 death
    print("Test 4: 1 death")
    tstart, tend = 0, 1
    tbirth_0 = [0]
    tau_0 = [1.5]
    simulator = WellMixedSimulator(f, ccm, Tdeath, tG2, tstart, tend)
    data = simulator.run(tau_0, tbirth_0, tG1_0, clone_0)
    mc_data = WellMixedSimulationData(data)
    for key, value in data.items():
        print('{}: {}'.format(key, value))
    print()

    assert mc_data.get_num_divisions() == 0
    assert mc_data.get_num_transitions() == 0
    assert mc_data.get_num_deaths() == 1
    assert mc_data.get_status() == 1
    assert np.all(mc_data.get_tbirth() == [0])
    assert np.all(mc_data.get_tG1() == [1.0])
    assert np.all(mc_data.get_clone() == [0])
    assert np.all(mc_data.get_died() == [True])
    assert np.all(mc_data.get_tdeath() == [0.5] )
    assert np.all(mc_data.get_divided() == [False] )
    assert np.all(mc_data.get_tdivision() == [np.inf])
    assert np.all(mc_data.get_transitioned() == [False])
    assert np.all(mc_data.get_ttransition() == [np.inf])
    assert np.all(mc_data.get_t_last_alive() == [0.5])
    assert np.all(mc_data.get_max_age() == [0.5])
    assert np.all(mc_data.get_time_in_G1() == [0.5])

    # Test 5: same, but with transition
    print("Test 5: 1 transition")
    tstart, tend = 0, 1
    tbirth_0 = [-0.5]
    tau_0 = [1.5]
    simulator = WellMixedSimulator(f, ccm, Tdeath, tG2, tstart, tend)
    data = simulator.run(tau_0, tbirth_0, tG1_0, clone_0)
    mc_data = WellMixedSimulationData(data)
    for key, value in data.items():
        print('{}: {}'.format(key, value))
    print()

    assert mc_data.get_num_divisions() == 0
    assert mc_data.get_num_transitions() == 1
    assert mc_data.get_num_deaths() == 0
    assert mc_data.get_status() == 0

    # Test 6: two cells, with constant base rate and fixed G1, no division or deaths
    print("Test 6: 2 cells, no divisions or deaths")
    tbirth_0 = [0, 0]
    tau_0 = [0, 0]
    tG1_0 = [1, 1]
    clone_0 = [0, 1]
    simulator = WellMixedSimulator(f, ccm, Tdeath, tG2, tstart, tend)
    data = simulator.run(tau_0, tbirth_0, tG1_0, clone_0)
    mc_data = WellMixedSimulationData(data)
    for key, value in data.items():
        print('{}: {}'.format(key, value))
    print()

    assert mc_data.get_num_divisions() == 0
    assert mc_data.get_num_transitions() == 0
    assert mc_data.get_num_deaths() == 0
    assert mc_data.get_status() == 0

    # Test 7: same, but with 1 division
    print("Test 7: 2 cells, 1 division")
    tbirth_0 = [-1.5, 0]
    tau_0 = [0, 0]
    tG1_0 = [1, 1]
    clone_0 = [0, 1]
    simulator = WellMixedSimulator(f, ccm, Tdeath, tG2, tstart, tend)
    data = simulator.run(tau_0, tbirth_0, tG1_0, clone_0)
    mc_data = WellMixedSimulationData(data)
    for key, value in data.items():
        print('{}: {}'.format(key, value))
    print()

    assert mc_data.get_num_divisions() == 1
    assert mc_data.get_num_transitions() == 0
    assert mc_data.get_num_deaths() == 0
    assert mc_data.get_status() == 0

    # Test 8: same, but with 2 divisions
    print("Test 8: 2 cells, 2 divisions, not simultaneous")
    tbirth_0 = [-1.5, -1.4]
    tau_0 = [0, 0]
    tG1_0 = [1, 1]
    clone_0 = [0, 1]
    simulator = WellMixedSimulator(f, ccm, Tdeath, tG2, tstart, tend)
    data = simulator.run(tau_0, tbirth_0, tG1_0, clone_0)
    mc_data = WellMixedSimulationData(data)
    for key, value in data.items():
        print('{}: {}'.format(key, value))
    print()

    assert mc_data.get_num_divisions() == 2
    assert mc_data.get_num_transitions() == 0
    assert mc_data.get_num_deaths() == 0
    assert mc_data.get_status() == 0

    # Test 9: same, but with 1 death
    print("Test 9: 2 cells, 1 death")
    tbirth_0 = [0, 0]
    tau_0 = [1.5, 0]
    tG1_0 = [1, 1]
    clone_0 = [0, 1]
    simulator = WellMixedSimulator(f, ccm, Tdeath, tG2, tstart, tend)
    data = simulator.run(tau_0, tbirth_0, tG1_0, clone_0)
    mc_data = WellMixedSimulationData(data)
    for key, value in data.items():
        print('{}: {}'.format(key, value))
    print()

    assert mc_data.get_num_divisions() == 0
    assert mc_data.get_num_transitions() == 0
    assert mc_data.get_num_deaths() == 1
    assert mc_data.get_status() == 0

    # Test 10: same, but with 2 deaths
    print("Test 10: 2 cells, 2 deaths, not simultaneous")
    tbirth_0 = [0, 0]
    tau_0 = [1.5, 1.6]
    tG1_0 = [1, 1]
    clone_0 = [0, 1]
    simulator = WellMixedSimulator(f, ccm, Tdeath, tG2, tstart, tend)
    data = simulator.run(tau_0, tbirth_0, tG1_0, clone_0)
    mc_data = WellMixedSimulationData(data)
    for key, value in data.items():
        print('{}: {}'.format(key, value))
    print()

    assert mc_data.get_num_divisions() == 0
    assert mc_data.get_num_transitions() == 0
    assert mc_data.get_num_deaths() == 2
    assert mc_data.get_status() == 1

    # Test 11: same, but with 2 deaths, simultaneously
    print("Test 11: 2 cells, 2 deaths, simultaneous")
    tbirth_0 = [0, 0]
    tau_0 = [1.5, 1.5]
    tG1_0 = [1, 1]
    clone_0 = [0, 1]
    simulator = WellMixedSimulator(f, ccm, Tdeath, tG2, tstart, tend)
    data = simulator.run(tau_0, tbirth_0, tG1_0, clone_0)
    mc_data = WellMixedSimulationData(data)
    for key, value in data.items():
        print('{}: {}'.format(key, value))
    print()

    assert mc_data.get_num_divisions() == 0
    assert mc_data.get_num_transitions() == 0
    assert mc_data.get_num_deaths() == 2
    assert mc_data.get_status() == 1

    # Test 12: same, but with 1 division, 1 death
    print("Test 12: 2 cells, 1 division, 1 death, simultaneous")
    tbirth_0 =  [-1.5, 0]
    tau_0 = [0, 1.5]
    tG1_0 = [1, 1]
    clone_0 = [0, 1]
    simulator = WellMixedSimulator(f, ccm, Tdeath, tG2, tstart, tend)
    data = simulator.run(tau_0, tbirth_0, tG1_0, clone_0)
    mc_data = WellMixedSimulationData(data)
    for key, value in data.items():
        print('{}: {}'.format(key, value))
    print()

    assert mc_data.get_num_divisions() == 1
    assert mc_data.get_num_transitions() == 0
    assert mc_data.get_num_deaths() == 1
    assert mc_data.get_status() == 0

    # Test 13: same, but with 1 division, 1 death
    print("Test 13: 1 cells, 1 division, 1 death, non-simultaneous")
    tbirth_0 =  [-1.4, 0]
    tau_0 = [0, 1.5]
    tG1_0 = [1, 1]
    clone_0 = [0, 1]
    simulator = WellMixedSimulator(f, ccm, Tdeath, tG2, tstart, tend)
    data = simulator.run(tau_0, tbirth_0, tG1_0, clone_0)
    mc_data = WellMixedSimulationData(data)
    for key, value in data.items():
        print('{}: {}'.format(key, value))
    print()

    assert mc_data.get_num_divisions() == 1
    assert mc_data.get_num_transitions() == 0
    assert mc_data.get_num_deaths() == 1
    assert mc_data.get_status() == 0

    # Test 14: same, but with 2 divisions
    print("Test 14: 2 cells, 2 divisions, negative control for death")
    tbirth_0 = [-1.5, -1.4]
    tau_0 = [1.9, 1.8]
    tG1_0 = [1, 1]
    clone_0 = [0, 1]
    simulator = WellMixedSimulator(f, ccm, Tdeath, tG2, tstart, tend)
    data = simulator.run(tau_0, tbirth_0, tG1_0, clone_0)
    mc_data = WellMixedSimulationData(data)
    for key, value in data.items():
        print('{}: {}'.format(key, value))
    print()

    assert mc_data.get_num_divisions() == 2
    assert mc_data.get_num_transitions() == 0
    assert mc_data.get_num_deaths() == 0
    assert mc_data.get_status() == 0

    # Test 15: same as 14, but with f_is_stepwise_constant = False
    print("Test 15: 2 cells, 2 divisions, negative control for death")
    tbirth_0 = [-1.5, -1.4]
    tau_0 = [1.9, 1.8]
    tG1_0 = [1, 1]
    clone_0 = [0, 1]
    simulator = WellMixedSimulator(f, ccm, Tdeath, tG2, tstart, tend, f_is_stepwise_constant=False)
    data = simulator.run(tau_0, tbirth_0, tG1_0, clone_0)
    mc_data = WellMixedSimulationData(data)
    for key, value in data.items():
        print('{}: {}'.format(key, value))
    print()

    assert mc_data.get_num_divisions() == 2
    assert mc_data.get_num_transitions() == 0
    assert mc_data.get_num_deaths() == 0
    assert mc_data.get_status() == 0

    # Test 16: same as 14, but with ccm that takes argument
    print("Test 16: 2 cells, 2 divisions, cell cycle model with arg")
    tbirth_0 = [-1.5, -1.4]
    tau_0 = [1.9, 1.8]
    tG1_0 = [1, 1]
    clone_0 = [0, 1]
    ccm = lambda random_state, clone, x: x
    simulator = WellMixedSimulator(f, ccm, Tdeath, tG2, tstart, tend, ccm_args=(.9,))
    data = simulator.run(tau_0, tbirth_0, tG1_0, clone_0)
    mc_data = WellMixedSimulationData(data)
    for key, value in data.items():
        print('{}: {}'.format(key, value))
    print()

    assert mc_data.get_num_divisions() == 2
    assert mc_data.get_num_transitions() == 0
    assert mc_data.get_num_deaths() == 0
    assert mc_data.get_status() == 0

    # Test 17: same as 14, but with f that takes argument
    print("Test 17: 2 cells, 2 divisions, death clock signal with arg")
    tbirth_0 = [-1.5, -1.4]
    tau_0 = [1.9, 1.8]
    tG1_0 = [1, 1]
    clone_0 = [0, 1]
    ccm = lambda random_state, clone: 1
    simulator = WellMixedSimulator(f, ccm, Tdeath, tG2, tstart, tend, f_args=(.9,))
    data = simulator.run(tau_0, tbirth_0, tG1_0, clone_0)
    mc_data = WellMixedSimulationData(data)
    for key, value in data.items():
        print('{}: {}'.format(key, value))
    print()

    assert mc_data.get_num_divisions() == 2
    assert mc_data.get_num_transitions() == 0
    assert mc_data.get_num_deaths() == 0
    assert mc_data.get_status() == 0

    # Test 18: same as 17, but with f_is_stepwise_constant = False
    print("Test 18: 2 cells, 2 divisions, death clock signal with arg")
    tbirth_0 = [-1.5, -1.4]
    tau_0 = [1.9, 1.8]
    tG1_0 = [1, 1]
    clone_0 = [0, 1]
    ccm = lambda random_state, clone: 1
    simulator = WellMixedSimulator(f, ccm, Tdeath, tG2, tstart, tend, f_args=(.9,),
            f_is_stepwise_constant=False)
    data = simulator.run(tau_0, tbirth_0, tG1_0, clone_0)
    mc_data = WellMixedSimulationData(data)
    for key, value in data.items():
        print('{}: {}'.format(key, value))
    print()

    assert mc_data.get_num_divisions() == 2
    assert mc_data.get_num_transitions() == 0
    assert mc_data.get_num_deaths() == 0
    assert mc_data.get_status() == 0

    # Test 19: exponential cell cycle model
    print("Test 19: exponential cell cycle model")
    from time import time
    t1 = time()
    simulator = WellMixedSimulator(tstart=0, tend=np.inf, max_cell_count=50)
    data = simulator.run()
    t2 = time()

    for key in data:
        print('{}: {}'.format(key, data[key]))
    print()

    mc_data = WellMixedSimulationData(data)

    assert mc_data.get_status() == 2

    cell_count = mc_data.get_cell_count()
    total_cell_count = mc_data.get_total_cell_count()
    final_cell_count = cell_count[-1]

    cell_count_alt = np.array([len(cell_indices) for cell_indices in data['cell_indices']])
    total_cell_count_alt = len(data['tbirth'])
    final_cell_count_alt = cell_count_alt[-1]

    t = mc_data.get_t_events()
    t_alt = np.array([t_grid[0] for t_grid in data['t_grid']])

    print("cell_count.shape: {}".format(cell_count.shape))
    print("cell_count_alt.shape: {}".format(cell_count_alt.shape))

    print("total_cell_count: {}".format(total_cell_count))
    print("total_cell_count_alt: {}".format(total_cell_count_alt))

    print("final_cell_count: {}".format(final_cell_count))
    print("final_cell_count_alt: {}".format(final_cell_count_alt))

    print("t.shape: {}".format(t.shape))
    print("t_alt.shape: {}".format(t_alt.shape))

    print("t t_alt: cell_count cell_count_alt")
    for cur_t, cur_t_alt, cur_cell_count, cur_cell_count_alt in \
            zip(t, t_alt, cell_count, cell_count_alt):
        print("{:.2f} {:.2f}: {} {}".format(
            cur_t,
            cur_t_alt,
            cur_cell_count,
            cur_cell_count_alt
            ))

    print("t: cell_count")
    for cur_t, cur_cell_count, in zip(t, cell_count):
        print("{:.2f}: {}".format(
            cur_t,
            cur_cell_count
            ))

    clone_header_str = "t:"
    clones = mc_data.get_unique_clones()
    for clone in clones:
        clone_header_str += " " + str(clone)
    print(clone_header_str)

    clone_cell_count = mc_data.get_clone_cell_count()
    clone_cell_count_0 = mc_data.get_cell_count_for_clone(0)
    clone_cell_count_1 = mc_data.get_cell_count_for_clone(1)
    assert np.all(clone_cell_count[:, 0] == clone_cell_count_0)
    assert np.all(clone_cell_count[:, 1] == clone_cell_count_1)
    for cur_t, cur_clone_cell_count in zip(t, clone_cell_count):
        print("{:.2f}: {}".format(cur_t, cur_clone_cell_count))

    print("Computation time: {}s".format(t2 - t1))

    print("Computation time per final cell: {}".format((t2 - t1) / cell_count[-1]))
    print("Computation time per total cell: {}".format((t2 - t1) / total_cell_count))

    # Test 20: uniform cell cycle model
    print("Test 20: uniform cell cycle model")
    tG1_0 = uniform.rvs(loc=50 - 0.5 * 20, scale=20, size=4) + 50
    t1 = time()
    simulator = WellMixedSimulator(tstart=0, tend=np.inf, max_cell_count=50,
            ccm=uniform_ccm)
    data = simulator.run(tG1_0=tG1_0)
    t2 = time()

    for key in data:
        print('{}: {}'.format(key, data[key]))
    print()

    mc_data = WellMixedSimulationData(data)

    assert mc_data.get_status() == 2

    cell_count = mc_data.get_cell_count()
    total_cell_count = mc_data.get_total_cell_count()
    final_cell_count = cell_count[-1]

    cell_count_alt = np.array([len(cell_indices) for cell_indices in data['cell_indices']])
    total_cell_count_alt = len(data['tbirth'])
    final_cell_count_alt = cell_count_alt[-1]

    t = mc_data.get_t_events()
    t_alt = np.array([t_grid[0] for t_grid in data['t_grid']])

    print("cell_count.shape: {}".format(cell_count.shape))
    print("cell_count_alt.shape: {}".format(cell_count_alt.shape))

    print("total_cell_count: {}".format(total_cell_count))
    print("total_cell_count_alt: {}".format(total_cell_count_alt))

    print("final_cell_count: {}".format(final_cell_count))
    print("final_cell_count_alt: {}".format(final_cell_count_alt))

    print("t.shape: {}".format(t.shape))
    print("t_alt.shape: {}".format(t_alt.shape))

    print("t t_alt: cell_count cell_count_alt")
    for cur_t, cur_t_alt, cur_cell_count, cur_cell_count_alt in \
            zip(t, t_alt, cell_count, cell_count_alt):
        print("{:.2f} {:.2f}: {} {}".format(
            cur_t,
            cur_t_alt,
            cur_cell_count,
            cur_cell_count_alt
            ))

    print("t: cell_count")
    for cur_t, cur_cell_count, in zip(t, cell_count):
        print("{:.2f}: {}".format(
            cur_t,
            cur_cell_count
            ))

    clone_header_str = "t:"
    clones = mc_data.get_unique_clones()
    for clone in clones:
        clone_header_str += " " + str(clone)
    print(clone_header_str)

    clone_cell_count = mc_data.get_clone_cell_count()
    clone_cell_count_0 = mc_data.get_cell_count_for_clone(0)
    clone_cell_count_1 = mc_data.get_cell_count_for_clone(1)
    assert np.all(clone_cell_count[:, 0] == clone_cell_count_0)
    assert np.all(clone_cell_count[:, 1] == clone_cell_count_1)
    for cur_t, cur_clone_cell_count in zip(t, clone_cell_count):
        print("{:.2f}: {}".format(cur_t, cur_clone_cell_count))

    print("Computation time: {}s".format(t2 - t1))

    print("Computation time per final cell: {}".format((t2 - t1) / cell_count[-1]))
    print("Computation time per total cell: {}".format((t2 - t1) / total_cell_count))

    # Test 21: exponential cell cycle model with normalised g2 signal
    print("Test 21: exponential cell cycle model with normalised g2 signal")
    t1 = time()
    simulator = WellMixedSimulator(tstart=0, tend=np.inf, max_cell_count=50,
            f=normalised_g2_death_signal, Tdeath=25)
    data = simulator.run()
    t2 = time()

    for key in data:
        print('{}: {}'.format(key, data[key]))
    print()

    mc_data = WellMixedSimulationData(data)

    cell_count = mc_data.get_cell_count()
    total_cell_count = mc_data.get_total_cell_count()
    final_cell_count = cell_count[-1]

    cell_count_alt = np.array([len(cell_indices) for cell_indices in data['cell_indices']])
    total_cell_count_alt = len(data['tbirth'])
    final_cell_count_alt = cell_count_alt[-1]

    t = mc_data.get_t_events()
    t_alt = np.array([t_grid[0] for t_grid in data['t_grid']])

    print("cell_count.shape: {}".format(cell_count.shape))
    print("cell_count_alt.shape: {}".format(cell_count_alt.shape))

    print("total_cell_count: {}".format(total_cell_count))
    print("total_cell_count_alt: {}".format(total_cell_count_alt))

    print("final_cell_count: {}".format(final_cell_count))
    print("final_cell_count_alt: {}".format(final_cell_count_alt))

    print("t.shape: {}".format(t.shape))
    print("t_alt.shape: {}".format(t_alt.shape))

    print("t: cell_count")
    for cur_t, cur_cell_count, in zip(t, cell_count):
        print("{:.2f}: {}".format(
            cur_t,
            cur_cell_count
            ))

    clone_header_str = "t:"
    clones = mc_data.get_unique_clones()
    for clone in clones:
        clone_header_str += " " + str(clone)
    print(clone_header_str)

    clone_cell_count = mc_data.get_clone_cell_count()
    clone_cell_count_0 = mc_data.get_cell_count_for_clone(0)
    clone_cell_count_1 = mc_data.get_cell_count_for_clone(1)
    assert np.all(clone_cell_count[:, 0] == clone_cell_count_0)
    assert np.all(clone_cell_count[:, 1] == clone_cell_count_1)
    for cur_t, cur_clone_cell_count in zip(t, clone_cell_count):
        print("{:.2f}: {}".format(cur_t, cur_clone_cell_count))

    print("Computation time: {}s".format(t2 - t1))

    print("Computation time per final cell: {}".format((t2 - t1) / cell_count[-1]))
    print("Computation time per total cell: {}".format((t2 - t1) / total_cell_count))

    # Test 22: exponential cell cycle model with absolute g2 signal
    print("Test 22: exponential cell cycle model with absolute g2 signal")
    t1 = time()
    simulator = WellMixedSimulator(tstart=0, tend=2000, max_cell_count=50,
            f=g2_death_signal, Tdeath=300)
    data = simulator.run()
    t2 = time()

    for key in data:
        print('{}: {}'.format(key, data[key]))
    print()

    mc_data = WellMixedSimulationData(data)

    cell_count = mc_data.get_cell_count()
    total_cell_count = mc_data.get_total_cell_count()
    final_cell_count = cell_count[-1]

    cell_count_alt = np.array([len(cell_indices) for cell_indices in data['cell_indices']])
    total_cell_count_alt = len(data['tbirth'])
    final_cell_count_alt = cell_count_alt[-1]

    t = mc_data.get_t_events()
    t_alt = np.array([t_grid[0] for t_grid in data['t_grid']])

    print("cell_count.shape: {}".format(cell_count.shape))
    print("cell_count_alt.shape: {}".format(cell_count_alt.shape))

    print("total_cell_count: {}".format(total_cell_count))
    print("total_cell_count_alt: {}".format(total_cell_count_alt))

    print("final_cell_count: {}".format(final_cell_count))
    print("final_cell_count_alt: {}".format(final_cell_count_alt))

    print("t.shape: {}".format(t.shape))
    print("t_alt.shape: {}".format(t_alt.shape))

    print("t: cell_count")
    for cur_t, cur_cell_count, in zip(t, cell_count):
        print("{:.2f}: {}".format(
            cur_t,
            cur_cell_count
            ))

    clone_header_str = "t:"
    clones = mc_data.get_unique_clones()
    for clone in clones:
        clone_header_str += " " + str(clone)
    print(clone_header_str)

    clone_cell_count = mc_data.get_clone_cell_count()
    clone_cell_count_0 = mc_data.get_cell_count_for_clone(0)
    clone_cell_count_1 = mc_data.get_cell_count_for_clone(1)
    assert np.all(clone_cell_count[:, 0] == clone_cell_count_0)
    assert np.all(clone_cell_count[:, 1] == clone_cell_count_1)
    for cur_t, cur_clone_cell_count in zip(t, clone_cell_count):
        print("{:.2f}: {}".format(cur_t, cur_clone_cell_count))

    print("Computation time: {}s".format(t2 - t1))

    print("Computation time per final cell: {}".format((t2 - t1) / cell_count[-1]))
    print("Computation time per total cell: {}".format((t2 - t1) / total_cell_count))

    # Print data
    print(mc_data)

    # Test 23: G1 proportion analysis function
    avg_tG1_fun = lambda eta, tG1: tG1 * (1 +
            (np.exp(-eta) - 1) / (np.exp(eta) - 1))

    gamma_fun = lambda beta, tG: beta * tG

    tG = 100
    Tdeath = 50
    c = 1

    N_beta = 10
    beta = np.arange( 1 / N_beta, 1, 1 / N_beta)

    f = average_f_g1_proportion(avg_tG1_fun, gamma_fun, tG, beta, Tdeath, c)

    print("f: {}".format(f))
    print("1 - beta: {}".format(1 - beta))

    theta_fun = lambda tG, beta, Tdeath, c, f: 1 - np.exp(- Tdeath / (c * f * beta * tG))
    theta_approx_fun = lambda tG, beta, Tdeath, c: theta_fun(tG, beta, Tdeath, c, 1 - beta)
    theta_approx_fun2 = lambda tG, beta, Tdeath, c: 1 - np.exp( - Tdeath / (c * beta * (1 - beta) * tG))

    print("theta(f): {}".format(theta_fun(tG, beta, Tdeath, c, f)))
    print("theta(beta): {}".format(theta_approx_fun(tG, beta, Tdeath, c)))
    print("theta2(beta): {}".format(theta_approx_fun2(tG, beta, Tdeath, c)))

    # Test 24: same as test 2, but with max_cell_count = 2
    print("Test 24: 1 division")
    f = base_rate_death_signal
    ccm = lambda random_state, clone: 1

    Tdeath = 2
    tG2 = 1
    tstart = 0
    tend = 1

    tau_0 = [0]
    tbirth_0 = [-1.5]
    tG1_0 = [1]
    clone_0 = [0]

    simulator = WellMixedSimulator(f, ccm, Tdeath, tG2, tstart, tend, max_cell_count=2)
    data = simulator.run(tau_0, tbirth_0, tG1_0, clone_0)
    mc_data = WellMixedSimulationData(data)
    for key, value in data.items():
        print('{}: {}'.format(key, value))
    print()

    assert mc_data.get_num_divisions() == 1
    assert mc_data.get_num_transitions() == 0
    assert mc_data.get_num_deaths() == 0
    assert mc_data.get_status() == 2
    assert np.all(mc_data.get_t_last_alive() == [0.5, 0.5, 0.5])

    # Test 25: test cumulative density functions
    assert np.isclose(exponential_cdf(0), 0)
    assert np.isclose(exponential_cdf(50), 1 - 1 / np.e)
    assert np.isclose(exponential_cdf(100), 1 - 1 / np.e / np.e)
    assert np.isclose(exponential_cdf(100, 100), 1 - 1 / np.e)
    assert np.isclose(exponential_cdf(200, 100), 1 - 1 / np.e / np.e)

    assert np.isclose(uniform_cdf(0), 0)
    assert np.isclose(uniform_cdf(50), 0.5)
    assert np.isclose(uniform_cdf(100), 1)
    assert np.isclose(uniform_cdf(100, 100), 0.5)
    assert np.isclose(uniform_cdf(200, 100), 1)

    assert np.isclose(uniform_cdf(55), 0.75)
    assert np.isclose(uniform_cdf(550, 500, 200), 0.75)

    # Test 26: interpolate_F_inverse
    t_grid = np.array([0, 1, 1, 2, 2, 3, 3, 4, 4, 5])
    F_grid = np.array([0, 1, 1, 1, 1, 1, 1, 2, 2, 2])
    assert np.isclose(interpolate_F_inverse(0, t_grid, F_grid), 0)
    assert np.isclose(interpolate_F_inverse(0.5, t_grid, F_grid), 0.5)
    assert np.isclose(interpolate_F_inverse(1, t_grid, F_grid), 1)
    assert np.isclose(interpolate_F_inverse(1.5, t_grid, F_grid), 3.5)
    assert np.isclose(interpolate_F_inverse(2, t_grid, F_grid), 4)
    assert np.isclose(interpolate_F_inverse(3, t_grid, F_grid), np.inf)

    # Test 27: zero G2 duration
    print("Test 27: zero G2 duration")
    f = base_rate_death_signal
    ccm = lambda random_state, clone: 1

    Tdeath = 2
    tG2 = 0
    tstart = 0
    tend = 2

    tau_0 = [0]
    tbirth_0 = [0]
    tG1_0 = [1]
    clone_0 = [0]

    simulator = WellMixedSimulator(f, ccm, Tdeath, tG2, tstart, tend)
    data = simulator.run(tau_0, tbirth_0, tG1_0, clone_0)
    mc_data = WellMixedSimulationData(data)
    for key, value in data.items():
        print('{}: {}'.format(key, value))
    print()

    assert mc_data.get_num_divisions() == 1
    assert mc_data.get_num_transitions() == 1
    assert mc_data.get_num_deaths() == 0
    assert mc_data.get_status() == 0

    # Test 28: instant death
    print("Test 28: instant death")
    base_rate_death_signal
    f = lambda a, b, c, d, e, f: base_rate_death_signal(a, b, c, d, e, f, base_rate=1000)
    ccm = lambda random_state, clone: 1

    Tdeath = 2
    tG2 = 1
    tstart = 0
    tend = 10
    min_cell_count = 10

    #tau_0 = np.zeros(128)
    tau_0 = np.linspace(0, 2 - 1/128, 128)
    tbirth_0 = np.zeros(128)
    tG1_0 = np.ones(128)
    clone_0 = np.zeros(128)

    simulator = WellMixedSimulator(f, ccm, Tdeath, tG2, tstart, tend,
            min_cell_count=min_cell_count)
    data = simulator.run(tau_0, tbirth_0, tG1_0, clone_0)
    mc_data = WellMixedSimulationData(data)
    for key, value in data.items():
        print('{}: {}'.format(key, value))
    print()

    assert mc_data.get_num_divisions() == 0
    assert mc_data.get_num_transitions() == 0
    assert mc_data.get_num_deaths() == 118
    assert mc_data.get_status() == 3

    print(mc_data)

    # Test 29: exponential cell cycle model with absolute g2 signal and seed
    print("Test 29: exponential cell cycle model with normalised g2 signal and seed")
    seed = 0
    tbirth_0 = np.zeros(4)
    tG1_0 = np.zeros(4)
    t1 = time()
    simulator = WellMixedSimulator(tstart=0, tend=1000, max_cell_count=100,
            min_cell_count=2,
            f=normalised_g2_death_signal, Tdeath=20)
    data = simulator.run(seed=seed, tbirth_0=tbirth_0, tG1_0=tG1_0)
    t2 = time()

    assert data['param']['f'] == normalised_g2_death_signal
    assert data['param']['ccm'] == exponential_ccm
    assert data['param']['Tdeath'](np.array([None])) == 20
    assert data['param']['tG2']((np.array([None]))) == 50
    assert data['param']['tstart'] == 0
    assert data['param']['tend'] == 1000
    assert data['param']['f_args'] == ()
    assert data['param']['ccm_args'] == ()
    assert data['param']['max_cell_count'] == 100
    assert data['param']['min_cell_count'] == 2
    assert data['param']['f_is_stepwise_constant'] == True
    assert np.all(data['init_cond']['tau_0'] == np.zeros(4))
    assert np.all(data['init_cond']['tbirth_0'] == tbirth_0)
    assert np.all(data['init_cond']['tG1_0'] == tG1_0)
    assert np.all(data['init_cond']['clone_0'] == np.arange(4))
    assert data['init_cond']['seed'] == seed

    for key in data:
        print('{}: {}'.format(key, data[key]))
    print()

    mc_data = WellMixedSimulationData(data)

    cell_count = mc_data.get_cell_count()
    total_cell_count = mc_data.get_total_cell_count()
    final_cell_count = cell_count[-1]

    cell_count_alt = np.array([len(cell_indices) for cell_indices in data['cell_indices']])
    total_cell_count_alt = len(data['tbirth'])
    final_cell_count_alt = cell_count_alt[-1]

    t = mc_data.get_t_events()
    t_alt = np.array([t_grid[0] for t_grid in data['t_grid']])

    print("cell_count.shape: {}".format(cell_count.shape))
    print("cell_count_alt.shape: {}".format(cell_count_alt.shape))

    print("total_cell_count: {}".format(total_cell_count))
    print("total_cell_count_alt: {}".format(total_cell_count_alt))

    print("final_cell_count: {}".format(final_cell_count))
    print("final_cell_count_alt: {}".format(final_cell_count_alt))

    print("t.shape: {}".format(t.shape))
    print("t_alt.shape: {}".format(t_alt.shape))

    print("t: cell_count")
    for cur_t, cur_cell_count, in zip(t, cell_count):
        print("{:.2f}: {}".format(
            cur_t,
            cur_cell_count
            ))

    clone_header_str = "t:"
    clones = mc_data.get_unique_clones()
    for clone in clones:
        clone_header_str += " " + str(clone)
    print(clone_header_str)

    clone_cell_count = mc_data.get_clone_cell_count()
    clone_cell_count_0 = mc_data.get_cell_count_for_clone(0)
    clone_cell_count_1 = mc_data.get_cell_count_for_clone(1)
    assert np.all(clone_cell_count[:, 0] == clone_cell_count_0)
    assert np.all(clone_cell_count[:, 1] == clone_cell_count_1)
    for cur_t, cur_clone_cell_count in zip(t, clone_cell_count):
        print("{:.2f}: {}".format(cur_t, cur_clone_cell_count))

    print("Computation time: {}s".format(t2 - t1))

    print("Computation time per final cell: {}".format((t2 - t1) / cell_count[-1]))
    print("Computation time per total cell: {}".format((t2 - t1) / total_cell_count))

    # Print data
    print(mc_data)

    G1_cell_count = mc_data.get_G1_cell_count()
    G2_cell_count = mc_data.get_G2_cell_count()
    combined_cell_count = G1_cell_count + G2_cell_count
    cell_count = mc_data.get_cell_count()
    assert np.all(combined_cell_count == cell_count)
    print(G1_cell_count)
    print(G2_cell_count)
    print(combined_cell_count)
    print(cell_count)

    import matplotlib.pyplot as plt

    plt.figure(figsize=(20,10))
    plt.step(t, cell_count, '-', label='cell count', where='post')
    plt.step(t, G1_cell_count, 'b-', label='G1', where='post')
    plt.step(t, G2_cell_count, 'r-', label='G2', where='post')
    plt.step(t, combined_cell_count, '--', label='combined cell count', where='post')
    plt.legend()
    plt.savefig('cell-count-test.png')
    plt.close()

    t_grid = mc_data.get_t_grid()
    tau_0 = mc_data.get_tau_for_cell_index(0)
    tau_4 = mc_data.get_tau_for_cell_index(4)
    tau_10 = mc_data.get_tau_for_cell_index(10)
    tau_54 = mc_data.get_tau_for_cell_index(54)

    last_tau = mc_data.get_last_tau()
    last_tau_0 = last_tau[0]
    last_tau_4 = last_tau[4]
    last_tau_10 = last_tau[10]
    last_tau_54 = last_tau[54]

    average_f = mc_data.get_average_f()
    average_f_0 = average_f[0]
    average_f_4 = average_f[4]
    average_f_10 = average_f[10]
    average_f_54 = average_f[54]

    time_in_G1 = mc_data.get_time_in_G1()
    time_in_G1_0 = time_in_G1[0]
    time_in_G1_4 = time_in_G1[4]
    time_in_G1_10 = time_in_G1[10]
    time_in_G1_54 = time_in_G1[54]

    tbirth = mc_data.get_tbirth()
    tbirth_0 = tbirth[0]
    tbirth_4 = tbirth[4]
    tbirth_10 = tbirth[10]
    tbirth_54 = tbirth[54]

    print(t_grid[:10])
    print(tau_0[:10])
    print(tau_4[:10])
    print(tau_10[:10])
    print(tau_54[:10])

    plt.figure(figsize=(20,10))
    plt.plot(t_grid, tau_0, 'g', label='0')
    plt.plot(t_grid, tau_4, 'r', label='4')
    plt.plot(t_grid, tau_10, 'c', label='10')
    plt.plot(t_grid, tau_54, 'y', label='54')

    plt.hlines(last_tau_0, t_grid[0], t_grid[-1], linestyles='dashed', label='0', colors='g')
    plt.hlines(last_tau_4, t_grid[0], t_grid[-1], linestyles='dashed', label='4', colors='r')
    plt.hlines(last_tau_10, t_grid[0], t_grid[-1], linestyles='dashed', label='10', colors='c')
    plt.hlines(last_tau_54, t_grid[0], t_grid[-1], linestyles='dashed', label='54', colors='y')

    plt.plot([tbirth_0, tbirth_0 + time_in_G1_0], [0, time_in_G1_0 * average_f_0], 'g-.', label='0')
    plt.plot([tbirth_4, tbirth_4 + time_in_G1_4], [0, time_in_G1_4 * average_f_4], 'r-.', label='4')
    plt.plot([tbirth_10, tbirth_10 + time_in_G1_10], [0, time_in_G1_10 * average_f_10], 'c-.', label='10')
    plt.plot([tbirth_54, tbirth_54 + time_in_G1_54], [0, time_in_G1_54 * average_f_54], 'y-.', label='54')

    plt.legend()
    plt.savefig('tau-test.png')
    plt.close()

    def f(*args):
        signal = normalised_g2_death_signal(*args)

        if len(signal) > 0:
            return signal[0]
        else:
            return 0

    cdf = exponential_cdf

    Nplot = 10**4
    t_plot = np.linspace(t_grid[0], t_grid[-1], Nplot)
    theta = mc_data.get_hindsight_survival_probability(t_plot, f, cdf)
    f_signal = mc_data.get_death_clock_signal(f)
    int_signal = mc_data.get_integrated_death_clock_signal(f)
    print(t_grid)
    print(f_signal)
    print(int_signal)

    fig, ax1 = plt.subplots(figsize=(20,10))
    ax2 = ax1.twinx()

    ax1.step(t_grid, f_signal / f_signal.max(), label='f')
    ax1.step(t, G2_cell_count / (G1_cell_count + G2_cell_count), label='G2 proportion',
            where='post')
    ax1.step(t_plot, theta, label=r'$\theta$')
    ax1.legend()
    ax1.set_ylabel('other')

    ax2.plot(t_grid, int_signal, label='f integrated')
    ax2.legend()
    ax2.set_ylabel('f integrated')

    fig.savefig('f-signal-test.png')
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(20,10))

    mc_data.plot_cell_history(ax, 100)

    fig.savefig('cell-history.png')
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(20,10))

    mc_data.plot_cell_counts(ax)

    fig.savefig('cell-counts.png')
    plt.close(fig)

    # Test 30: two cells with constant base rate and heterotypic cell cycle model
    print("Test 30: heterotypic cell cycle model")
    f = base_rate_death_signal
    ccm = lambda random_state, clone: clone + 1

    Tdeath = 2
    tG2 = 1
    tstart = 0
    tend = 3.5

    tau_0 = [0, 0]
    tbirth_0 = [0, 0]
    tG1_0 = [1, 1]
    clone_0 = [0, 1]

    simulator = WellMixedSimulator(f, ccm, Tdeath, tG2, tstart, tend)
    data = simulator.run(tau_0, tbirth_0, tG1_0, clone_0)
    mc_data = WellMixedSimulationData(data)
    for key, value in data.items():
        print('{}: {}'.format(key, value))
    print()

    assert mc_data.get_num_divisions() == 2
    assert mc_data.get_num_divisions_for_clone(0) == 1
    assert mc_data.get_num_divisions_for_clone(1) == 1
    assert mc_data.get_num_transitions() == 4
    assert mc_data.get_num_transitions_for_clone(0) == 3
    assert mc_data.get_num_transitions_for_clone(1) == 1
    assert mc_data.get_num_deaths() == 0
    assert mc_data.get_num_deaths_for_clone(0) == 0
    assert mc_data.get_num_deaths_for_clone(1) == 0
    assert mc_data.get_status() == 0
    assert np.all(mc_data.get_tbirth() == [0.0, 0.0, 2.0, 2.0, 2.0, 2.0])
    assert np.all(mc_data.get_tG1() == [1.0, 1.0, 1.0, 1.0, 2.0, 2.0])
    assert np.all(mc_data.get_clone() == [0, 1, 0, 0, 1, 1])
    assert np.all(mc_data.get_died() == [False] * 6)
    assert np.all(mc_data.get_tdeath() == [np.inf] * 6)
    assert np.all(mc_data.get_divided() == [True, True, False, False, False, False])
    assert np.all(mc_data.get_tdivision() == [2.0, 2.0, np.inf, np.inf, np.inf, np.inf])
    assert np.all(mc_data.get_transitioned() == [True, True, True, True, False, False])
    assert np.all(mc_data.get_ttransition() == [1.0, 1.0, 3.0, 3.0, np.inf, np.inf])
    assert np.all(mc_data.get_t_last_alive() == [2.0, 2.0, 3.5, 3.5, 3.5, 3.5])
    assert np.all(mc_data.get_max_age() == [2, 2, 1.5, 1.5, 1.5, 1.5])
    assert np.all(mc_data.get_time_in_G1() == [1, 1, 1, 1, 1.5, 1.5])
    assert np.all(mc_data.get_time_in_G2() == [1, 1, 0.5, 0.5, 0.0, 0.0])

    # Test 31: two cells with constant base rate and heterotypic death threshold
    print("Test 31: heterotypic death threshold")
    f = base_rate_death_signal
    ccm = lambda random_state, clone: 1

    Tdeath = lambda clone: clone + 1
    tG2 = 1
    tstart = 0
    tend = 1.5

    tau_0 = [0, 0]
    tbirth_0 = [0, 0]
    tG1_0 = [2, 2]
    clone_0 = [0, 1]

    simulator = WellMixedSimulator(f, ccm, Tdeath, tG2, tstart, tend)
    data = simulator.run(tau_0, tbirth_0, tG1_0, clone_0)
    mc_data = WellMixedSimulationData(data)
    for key, value in data.items():
        print('{}: {}'.format(key, value))
    print()

    assert mc_data.get_num_divisions() == 0
    assert mc_data.get_num_divisions_for_clone(0) == 0
    assert mc_data.get_num_divisions_for_clone(1) == 0
    assert mc_data.get_num_transitions() == 0
    assert mc_data.get_num_transitions_for_clone(0) == 0
    assert mc_data.get_num_transitions_for_clone(1) == 0
    assert mc_data.get_num_deaths() == 1
    assert mc_data.get_num_deaths_for_clone(0) == 1
    assert mc_data.get_num_deaths_for_clone(1) == 0
    assert mc_data.get_status() == 0
    assert np.all(mc_data.get_tbirth() == [0.0, 0.0])
    assert np.all(mc_data.get_tG1() == [2.0, 2.0])
    assert np.all(mc_data.get_clone() == [0, 1])
    assert np.all(mc_data.get_died() == [True, False])
    assert np.all(mc_data.get_tdeath() == [1.0, np.inf])
    assert np.all(mc_data.get_divided() == [False, False])
    assert np.all(mc_data.get_tdivision() == [np.inf, np.inf])
    assert np.all(mc_data.get_transitioned() == [False, False])
    assert np.all(mc_data.get_ttransition() == [np.inf, np.inf])
    assert np.all(mc_data.get_t_last_alive() == [1.0, 1.5])
    assert np.all(mc_data.get_max_age() == [1, 1.5])
    assert np.all(mc_data.get_time_in_G1() == [1, 1.5])
    assert np.all(mc_data.get_time_in_G2() == [0, 0])

    # Test 32: two cells with constant base rate and heterotypic G2 duration
    print("Test 32: heterotypic G2 duration")
    f = base_rate_death_signal
    ccm = lambda random_state, clone: 1

    Tdeath = 2
    tG2 = lambda clone: clone + 1
    tstart = 0
    tend = 2.5

    tau_0 = [0, 0]
    tbirth_0 = [0, 0]
    tG1_0 = [1, 1]
    clone_0 = [0, 1]

    simulator = WellMixedSimulator(f, ccm, Tdeath, tG2, tstart, tend)
    data = simulator.run(tau_0, tbirth_0, tG1_0, clone_0)
    mc_data = WellMixedSimulationData(data)
    for key, value in data.items():
        print('{}: {}'.format(key, value))
    print()

    assert mc_data.get_num_divisions() == 1
    assert mc_data.get_num_divisions_for_clone(0) == 1
    assert mc_data.get_num_divisions_for_clone(1) == 0
    assert mc_data.get_num_transitions() == 2
    assert mc_data.get_num_transitions_for_clone(0) == 1
    assert mc_data.get_num_transitions_for_clone(1) == 1
    assert mc_data.get_num_deaths() == 0
    assert mc_data.get_num_deaths_for_clone(0) == 0
    assert mc_data.get_num_deaths_for_clone(1) == 0
    assert mc_data.get_status() == 0
    assert np.all(mc_data.get_tbirth() == [0.0, 0.0, 2, 2])
    assert np.all(mc_data.get_tG1() == [1, 1, 1, 1])
    assert np.all(mc_data.get_clone() == [0, 1, 0, 0])
    assert np.all(mc_data.get_died() == [False] * 4)
    assert np.all(mc_data.get_tdeath() == [np.inf] * 4)
    assert np.all(mc_data.get_divided() == [True, False, False, False])
    assert np.all(mc_data.get_tdivision() == [2, np.inf, np.inf, np.inf])
    assert np.all(mc_data.get_transitioned() == [True, True, False, False])
    assert np.all(mc_data.get_ttransition() == [1, 1, np.inf, np.inf])
    assert np.all(mc_data.get_t_last_alive() == [2, 2.5, 2.5, 2.5])
    assert np.all(mc_data.get_max_age() == [2, 2.5, 0.5, 0.5])
    assert np.all(mc_data.get_time_in_G1() == [1, 1, 0.5, 0.5])
    assert np.all(mc_data.get_time_in_G2() == [1, 1.5, 0, 0])
    assert np.all(mc_data.get_G1_cell_count() == mc_data.get_G1_cell_count_for_clone(0)
            +  mc_data.get_G1_cell_count_for_clone(1))
    assert np.all(mc_data.get_G2_cell_count() == mc_data.get_G2_cell_count_for_clone(0)
            +  mc_data.get_G2_cell_count_for_clone(1))

    # Test 33: test heterotypic cell cycle models and death signal functions
    print("Test 33: heterotypic cell cycle models and death signal functions")
    # Cell cycle models
    # This assert has a very small of failing even if there is no problem
    assert np.mean([exponential_ccm_heterotypic(clone=1, tG1_param_clone_1=5) for
            i in range(1000)]) < 50

    assert np.mean([uniform_ccm_heterotypic(clone=1, tG1_param_clone_1=100, r_clone_1=1) for
            i in range(1000)]) > 99.5

    # Death signal functions
    # Base rate
    try:
        base_rate_death_signal_heterotypic(
            t = np.array([0]),
            tau = np.array([0]),
            tbirth = np.array([0]),
            tG1 = np.array([0]),
            clone = np.array([2]),
            isinG1 = np.array([False]),
            base_rate_clone_1=5)
    except AssertionError as err:
        pass
    except:
        raise

    assert np.all(base_rate_death_signal_heterotypic(
            t = np.array([0]),
            tau = np.array([0, 0]),
            tbirth = np.array([0, 0]),
            tG1 = np.array([0, 0]),
            clone = np.array([0, 1]),
            isinG1 = np.array([False, False]),
            base_rate_clone_1=5)
            == np.array([1, 5]))

    # Normalised G2 neighbours death signal
    try:
        normalised_g2_death_signal_heterotypic(
            t = np.array([0]),
            tau = np.array([0]),
            tbirth = np.array([0]),
            tG1 = np.array([0]),
            clone = np.array([2]),
            isinG1 = np.array([False]),
            coef_clone_1=5)
    except AssertionError as err:
        pass
    except:
        raise

    assert np.all(normalised_g2_death_signal_heterotypic(
            t = np.array([0]),
            tau = np.array([0, 0]),
            tbirth = np.array([0, 0]),
            tG1 = np.array([0, 0]),
            clone = np.array([0, 1]),
            isinG1 = np.array([True, True]),
            coef_clone_1=5)
            == np.array([0., 0.]))

    assert np.all(normalised_g2_death_signal_heterotypic(
            t = np.array([0]),
            tau = np.array([0, 0]),
            tbirth = np.array([0, 0]),
            tG1 = np.array([0, 0]),
            clone = np.array([0, 1]),
            isinG1 = np.array([True, False]),
            coef_clone_1=5)
            == np.array([1., 5.]))

    assert np.all(normalised_g2_death_signal_heterotypic(
            t = np.array([0]),
            tau = np.array([0, 0]),
            tbirth = np.array([0, 0]),
            tG1 = np.array([0, 0]),
            clone = np.array([0, 1]),
            isinG1 = np.array([False, False]),
            coef_clone_1=5)
            == np.array([1., 5.]))

    # Absolute G2 neighbours death signal
    try:
        g2_death_signal_heterotypic(
            t = np.array([0]),
            tau = np.array([0]),
            tbirth = np.array([0]),
            tG1 = np.array([0]),
            clone = np.array([2]),
            isinG1 = np.array([False]),
            coef_clone_1=5)
    except AssertionError as err:
        pass
    except:
        raise

    assert np.all(g2_death_signal_heterotypic(
            t = np.array([0]),
            tau = np.array([0, 0]),
            tbirth = np.array([0, 0]),
            tG1 = np.array([0, 0]),
            clone = np.array([0, 1]),
            isinG1 = np.array([True, True]),
            coef_clone_1=5)
            == np.array([0., 0.]))

    assert np.all(g2_death_signal_heterotypic(
            t = np.array([0]),
            tau = np.array([0, 0]),
            tbirth = np.array([0, 0]),
            tG1 = np.array([0, 0]),
            clone = np.array([0, 1]),
            isinG1 = np.array([True, False]),
            coef_clone_1=5)
            == np.array([1., 5.]))

    assert np.all(g2_death_signal_heterotypic(
            t = np.array([0]),
            tau = np.array([0, 0]),
            tbirth = np.array([0, 0]),
            tG1 = np.array([0, 0]),
            clone = np.array([0, 1]),
            isinG1 = np.array([False, False]),
            coef_clone_1=5)
            == np.array([2., 10.]))

    # Test 34: two cells with heterotypic base rate
    print("Test 34: heterotypic base rate")
    f = base_rate_death_signal_heterotypic
    f_args = (1, 2)
    ccm = lambda random_state, clone: 1

    Tdeath = 1
    tG2 = 1
    tstart = 0
    tend = 1

    tau_0 = [0, 0]
    tbirth_0 = [0, 0]
    tG1_0 = [2, 2]
    clone_0 = [0, 1]

    simulator = WellMixedSimulator(f, ccm, Tdeath, tG2, tstart, tend, f_args)
    data = simulator.run(tau_0, tbirth_0, tG1_0, clone_0)
    mc_data = WellMixedSimulationData(data)
    for key, value in data.items():
        print('{}: {}'.format(key, value))
    print()

    assert mc_data.get_num_divisions() == 0
    assert mc_data.get_num_divisions_for_clone(0) == 0
    assert mc_data.get_num_divisions_for_clone(1) == 0
    assert mc_data.get_num_transitions() == 0
    assert mc_data.get_num_transitions_for_clone(0) == 0
    assert mc_data.get_num_transitions_for_clone(1) == 0
    assert mc_data.get_num_deaths() == 1
    assert mc_data.get_num_deaths_for_clone(0) == 0
    assert mc_data.get_num_deaths_for_clone(1) == 1
    assert mc_data.get_status() == 0
    assert np.all(mc_data.get_tbirth() == [0.0, 0.0])
    assert np.all(mc_data.get_tG1() == [2.0, 2.0])
    assert np.all(mc_data.get_clone() == [0, 1])
    assert np.all(mc_data.get_died() == [False, True])
    assert np.all(mc_data.get_tdeath() == [np.inf, 0.5])
    assert np.all(mc_data.get_divided() == [False, False])
    assert np.all(mc_data.get_tdivision() == [np.inf, np.inf])
    assert np.all(mc_data.get_transitioned() == [False, False])
    assert np.all(mc_data.get_ttransition() == [np.inf, np.inf])
    assert np.all(mc_data.get_t_last_alive() == [1.0, 0.5])
    assert np.all(mc_data.get_max_age() == [1.0, 0.5])
    assert np.all(mc_data.get_time_in_G1() == [1.0, 0.5])
    assert np.all(mc_data.get_time_in_G2() == [0, 0])
    assert np.all(mc_data.get_G1_cell_count() == mc_data.get_G1_cell_count_for_clone(0)
            +  mc_data.get_G1_cell_count_for_clone(1))
    assert np.all(mc_data.get_G2_cell_count() == mc_data.get_G2_cell_count_for_clone(0)
            +  mc_data.get_G2_cell_count_for_clone(1))

    # Test 35: max cell count for clone 0
    print("Test 35: max cell count for clone 0")
    f = base_rate_death_signal
    ccm = lambda random_state, clone: 1

    Tdeath = 2
    tG2 = 1
    tstart = 0
    tend = 10
    max_cell_count_for_clone = {0 : 2}

    tau_0 = np.zeros(3)
    tbirth_0 = np.zeros(3)
    tG1_0 = np.ones(3)
    clone_0 = np.array([0, 1, 1])

    simulator = WellMixedSimulator(f, ccm, Tdeath, tG2, tstart, tend,
            max_cell_count_for_clone=max_cell_count_for_clone)
    data = simulator.run(tau_0, tbirth_0, tG1_0, clone_0)
    mc_data = WellMixedSimulationData(data)
    for key, value in data.items():
        print('{}: {}'.format(key, value))
    print()

    assert mc_data.get_num_divisions() == 3
    assert mc_data.get_num_divisions_for_clone(0) == 1
    assert mc_data.get_num_divisions_for_clone(1) == 2
    assert mc_data.get_num_transitions() == 3
    assert mc_data.get_num_transitions_for_clone(0) == 1
    assert mc_data.get_num_transitions_for_clone(1) == 2
    assert mc_data.get_num_deaths() == 0
    assert mc_data.get_num_deaths_for_clone(0) == 0
    assert mc_data.get_num_deaths_for_clone(1) == 0
    assert mc_data.get_status() == 4
    assert mc_data.get_status_info() == 0
    assert np.all(mc_data.get_tbirth() == ([0.0] * 3 + [2.0] * 6))
    assert np.all(mc_data.get_tG1() == [1.0] * 9)
    assert np.all(mc_data.get_clone() == [0, 1, 1, 0, 0, 1, 1, 1, 1])
    assert np.all(mc_data.get_died() == [False] * 9)
    assert np.all(mc_data.get_tdeath() == [np.inf] * 9)
    assert np.all(mc_data.get_divided() == [True] * 3 + [False] * 6)
    assert np.all(mc_data.get_tdivision() == [2.0] * 3 + [np.inf] * 6)
    assert np.all(mc_data.get_transitioned() == [True] * 3 + [False] * 6)
    assert np.all(mc_data.get_ttransition() == [1.0] * 3 + [np.inf] * 6)
    assert np.all(mc_data.get_t_last_alive() == [2.0] * 9)
    assert np.all(mc_data.get_max_age() == [2.0] * 3 + [0.0] * 6)
    assert np.all(mc_data.get_time_in_G1() == [1.0] * 3 + [0.0] * 6)
    assert np.all(mc_data.get_time_in_G2() == [1.0] * 3 + [0.0] * 6)
    assert np.all(mc_data.get_G1_cell_count() == mc_data.get_G1_cell_count_for_clone(0)
            +  mc_data.get_G1_cell_count_for_clone(1))
    assert np.all(mc_data.get_G2_cell_count() == mc_data.get_G2_cell_count_for_clone(0)
            +  mc_data.get_G2_cell_count_for_clone(1))

    # Test 36: min cell count for clone 1
    print("Test 35: min cell count for clone 1")
    f = base_rate_death_signal
    ccm = lambda random_state, clone: 1

    Tdeath = 2
    tG2 = 1
    tstart = 0
    tend = 10
    min_cell_count_for_clone = {1 : 1}

    tau_0 = np.array([0, 0, 1.5])
    tbirth_0 = np.zeros(3)
    tG1_0 = np.ones(3)
    clone_0 = np.array([0, 1, 1])

    simulator = WellMixedSimulator(f, ccm, Tdeath, tG2, tstart, tend,
            min_cell_count_for_clone=min_cell_count_for_clone)
    data = simulator.run(tau_0, tbirth_0, tG1_0, clone_0)
    mc_data = WellMixedSimulationData(data)
    for key, value in data.items():
        print('{}: {}'.format(key, value))
    print()

    assert mc_data.get_num_divisions() == 0
    assert mc_data.get_num_divisions_for_clone(0) == 0
    assert mc_data.get_num_divisions_for_clone(1) == 0
    assert mc_data.get_num_transitions() == 0
    assert mc_data.get_num_transitions_for_clone(0) == 0
    assert mc_data.get_num_transitions_for_clone(1) == 0
    assert mc_data.get_num_deaths() == 1
    assert mc_data.get_num_deaths_for_clone(0) == 0
    assert mc_data.get_num_deaths_for_clone(1) == 1
    assert mc_data.get_status() == 5
    assert mc_data.get_status_info() == 1
    assert np.all(mc_data.get_tbirth() == ([0.0] * 3))
    assert np.all(mc_data.get_tG1() == [1.0] * 3)
    assert np.all(mc_data.get_clone() == [0, 1, 1])
    assert np.all(mc_data.get_died() == [False, False, True])
    assert np.all(mc_data.get_tdeath() == [np.inf, np.inf, 0.5])
    assert np.all(mc_data.get_divided() == [False] * 3)
    assert np.all(mc_data.get_tdivision() == [np.inf] * 3)
    assert np.all(mc_data.get_transitioned() == [False] * 3)
    assert np.all(mc_data.get_ttransition() == [np.inf] * 3)
    assert np.all(mc_data.get_t_last_alive() == [0.5] * 3)
    assert np.all(mc_data.get_max_age() == [0.5] * 3)
    assert np.all(mc_data.get_time_in_G1() == [0.5] * 3)
    assert np.all(mc_data.get_time_in_G2() == [0.0] * 3)
    assert np.all(mc_data.get_G1_cell_count() == mc_data.get_G1_cell_count_for_clone(0)
            +  mc_data.get_G1_cell_count_for_clone(1))
    assert np.all(mc_data.get_G2_cell_count() == mc_data.get_G2_cell_count_for_clone(0)
            +  mc_data.get_G2_cell_count_for_clone(1))

    # Test 36: same as test 4, but with apoptosis_at_checkpoint
    print("Test 36: 1 death with apoptosis_at_checkpoint")
    f = base_rate_death_signal
    ccm = lambda random_state, clone: 1

    Tdeath = 2
    tG2 = 1

    tG1_0 = [1]
    clone_0 = [0]

    tstart, tend = 0, 2
    tbirth_0 = [0]
    tau_0 = [1.5]
    simulator = WellMixedSimulator(f, ccm, Tdeath, tG2, tstart, tend,
            apoptosis_at_checkpoint=True)
    data = simulator.run(tau_0, tbirth_0, tG1_0, clone_0)
    mc_data = WellMixedSimulationData(data)
    for key, value in data.items():
        print('{}: {}'.format(key, value))
    print()

    assert mc_data.get_num_divisions() == 0
    assert mc_data.get_num_transitions() == 0
    assert mc_data.get_num_deaths() == 1
    assert mc_data.get_status() == 1
    assert np.all(mc_data.get_tbirth() == [0])
    assert np.all(mc_data.get_tG1() == [1.0])
    assert np.all(mc_data.get_clone() == [0])
    assert np.all(mc_data.get_died() == [True])
    assert np.all(mc_data.get_tdeath() == [1.0] )
    assert np.all(mc_data.get_divided() == [False] )
    assert np.all(mc_data.get_tdivision() == [np.inf])
    assert np.all(mc_data.get_transitioned() == [False])
    assert np.all(mc_data.get_ttransition() == [np.inf])
    assert np.all(mc_data.get_t_last_alive() == [1.0])
    assert np.all(mc_data.get_max_age() == [1.0])
    assert np.all(mc_data.get_time_in_G1() == [1.0])

    # Test 37: same as test 36, but with initial death clock over threshold
    print("Test 37: 1 death with apoptosis_at_checkpoint and initial death clock over threshold")
    f = base_rate_death_signal
    ccm = lambda random_state, clone: 1

    Tdeath = 2
    tG2 = 1

    tG1_0 = [1]
    clone_0 = [0]

    tstart, tend = 0, 2
    tbirth_0 = [0]
    tau_0 = [3]
    simulator = WellMixedSimulator(f, ccm, Tdeath, tG2, tstart, tend,
            apoptosis_at_checkpoint=True)
    data = simulator.run(tau_0, tbirth_0, tG1_0, clone_0)
    mc_data = WellMixedSimulationData(data)
    for key, value in data.items():
        print('{}: {}'.format(key, value))
    print()

    assert mc_data.get_num_divisions() == 0
    assert mc_data.get_num_transitions() == 0
    assert mc_data.get_num_deaths() == 1
    assert mc_data.get_status() == 1
    assert np.all(mc_data.get_tbirth() == [0])
    assert np.all(mc_data.get_tG1() == [1.0])
    assert np.all(mc_data.get_clone() == [0])
    assert np.all(mc_data.get_died() == [True])
    assert np.all(mc_data.get_tdeath() == [1.0] )
    assert np.all(mc_data.get_divided() == [False] )
    assert np.all(mc_data.get_tdivision() == [np.inf])
    assert np.all(mc_data.get_transitioned() == [False])
    assert np.all(mc_data.get_ttransition() == [np.inf])
    assert np.all(mc_data.get_t_last_alive() == [1.0])
    assert np.all(mc_data.get_max_age() == [1.0])
    assert np.all(mc_data.get_time_in_G1() == [1.0])

    # Test 38: same as test 1 but with apoptosis_at_checkpoint
    print("Test 38: no divisions or deaths with apoptosis_at_checkpoint")
    f = base_rate_death_signal
    ccm = lambda random_state, clone: 1

    Tdeath = 2
    tG2 = 1
    tstart = 0
    tend = 1

    tau_0 = [0]
    tbirth_0 = [0]
    tG1_0 = [1]
    clone_0 = [0]

    simulator = WellMixedSimulator(f, ccm, Tdeath, tG2, tstart, tend,
            apoptosis_at_checkpoint=True)
    data = simulator.run(tau_0, tbirth_0, tG1_0, clone_0)
    mc_data = WellMixedSimulationData(data)
    for key, value in data.items():
        print('{}: {}'.format(key, value))
    print()

    assert mc_data.get_num_divisions() == 0
    assert mc_data.get_num_transitions() == 0
    assert mc_data.get_num_deaths() == 0
    assert mc_data.get_status() == 0

    # Test 39: same, but with 1 division
    print("Test 39: 1 division with apoptosis_at_checkpoint")
    tbirth_0 = [-1.5]

    simulator = WellMixedSimulator(f, ccm, Tdeath, tG2, tstart, tend,
            apoptosis_at_checkpoint=True)
    data = simulator.run(tau_0, tbirth_0, tG1_0, clone_0)
    mc_data = WellMixedSimulationData(data)
    for key, value in data.items():
        print('{}: {}'.format(key, value))
    print()

    assert mc_data.get_num_divisions() == 1
    assert mc_data.get_num_transitions() == 0
    assert mc_data.get_num_deaths() == 0
    assert mc_data.get_status() == 0

    # Test 40: same, but with 3 divisions and 2 transitions, also test different
    # clone number
    print("Test 40: 3 divisions and 2 transitions with apoptosis_at_checkpoint")
    tstart, tend = 0, 3
    simulator = WellMixedSimulator(f, ccm, Tdeath, tG2, tstart, tend,
            apoptosis_at_checkpoint=True)
    data = simulator.run(tau_0, tbirth_0, tG1_0, [1])
    mc_data = WellMixedSimulationData(data)
    for key, value in data.items():
        print('{}: {}'.format(key, value))
    print()

    assert mc_data.get_num_divisions() == 3
    assert mc_data.get_num_transitions() == 2
    assert mc_data.get_num_deaths() == 0
    assert mc_data.get_status() == 0
    assert np.all(mc_data.get_tbirth() == [-1.5, 0.5, 0.5, 2.5, 2.5, 2.5, 2.5])
    assert np.all(mc_data.get_tG1() == [1.0] * 7)
    assert np.all(mc_data.get_clone() == [1] * 7)
    assert np.all(mc_data.get_died() == [False] * 7)
    assert np.all(mc_data.get_tdeath() == [np.inf] * 7)
    assert np.all(mc_data.get_divided() == [True, True, True, False, False, False, False])
    assert np.all(mc_data.get_tdivision() == [0.5, 2.5, 2.5, np.inf, np.inf, np.inf, np.inf])
    assert np.all(mc_data.get_transitioned() == [True, True, True, False, False, False, False])
    assert np.all(mc_data.get_ttransition() == [-0.5, 1.5, 1.5, np.inf, np.inf, np.inf, np.inf])
    assert np.all(mc_data.get_t_last_alive() == [0.5, 2.5, 2.5, 3, 3, 3, 3])
    assert np.all(mc_data.get_max_age() == [2, 2, 2, 0.5, 0.5, 0.5, 0.5])
    assert np.all(mc_data.get_time_in_G1() == [1, 1, 1, 0.5, 0.5, 0.5, 0.5,  ])
    assert np.all(mc_data.get_time_in_G2() == [1, 1, 1, 0.0, 0.0, 0.0, 0.0,  ])
    assert np.all(mc_data.get_effective_time_in_G1() == [1, 1, 1])

    # Test 41: 1 death before transition
    print("Test 41: 1 death before transition with apoptosis_at_checkpoint")
    tstart, tend = 0, 1
    tbirth_0 = [-0.5]
    tau_0 = [1.5]
    simulator = WellMixedSimulator(f, ccm, Tdeath, tG2, tstart, tend,
            apoptosis_at_checkpoint=True)
    data = simulator.run(tau_0, tbirth_0, tG1_0, clone_0)
    mc_data = WellMixedSimulationData(data)
    for key, value in data.items():
        print('{}: {}'.format(key, value))
    print()

    assert mc_data.get_num_divisions() == 0
    assert mc_data.get_num_transitions() == 0
    assert mc_data.get_num_deaths() == 1
    assert mc_data.get_status() == 1

    # Test 42: ergodic_rms
    print("Test 42: ergodic_rms")
    Tdeath = 100
    tstart, tend = 0, 1
    tG2 = 1

    ccm = lambda random_state, clone: 1
    f = lambda t, tau, tbirth, tG1, clone, isinG1: t * np.ones(tau.shape)

    tbirth_0 = [0]
    tau_0 = [0]
    tG1_0 = [1]
    clone_0 = [1]

    simulator = WellMixedSimulator(f, ccm, Tdeath, tG2, tstart, tend)
    data = simulator.run(tau_0, tbirth_0, tG1_0, clone_0)
    mc_data = WellMixedSimulationData(data)
    for key, value in data.items():
        print('{}: {}'.format(key, value))
    print()

    assert mc_data.get_num_divisions() == 0
    assert mc_data.get_num_transitions() == 0
    assert mc_data.get_num_deaths() == 0
    assert mc_data.get_status() == 0
    assert mc_data.get_f() == f

    beta = 0.7
    c = .3

    predicted_rms = np.sqrt(((1 - beta)**2 + (1/c - (1 - beta))**2) / 2)
    computed_rms = mc_data.get_ergodic_rms(f, c, beta)

    print('predicted_rms: {}'.format(predicted_rms))
    print('computed_rms: {}'.format(computed_rms))

    assert np.isclose(predicted_rms, computed_rms)
