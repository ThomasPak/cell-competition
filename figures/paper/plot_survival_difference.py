import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.colors as colors

upper_colour = 'tab:grey'
lower_colour = 'tab:cyan'
viability_colour = 'black'
rho_colour = 'black'

lower_style = 'dashed'
upper_style = 'dashed'
rho_style = 'dashed'

B_loser_viability_curve_colour = 'tab:cyan'
A_loser_viability_curve_colour = 'tab:grey'

def plot_heterotypic_survival_difference(ax, df, beta_A, eta_A):
    """
    Plots survival difference against eta_B, beta_B on given ax.

    Parameters
        ax: ax to plot survival difference on
        df: pandas.DataFrame with diff_theta column
        beta_A: value of beta_A to plot predicted border
        eta_A: value of eta_A to plot predicted border
        colourbar
    """

    # Compute average survival difference
    grouped = df.groupby(['eta_B', 'beta_B']).mean()['diff_theta']

    data = grouped.unstack().sort_index(ascending=False)
    data_np = data.to_numpy()

    # Get max_abs_val
    max_abs_val = np.absolute(data_np).max()

    # Plot data
    im = ax.imshow(data_np, cmap='seismic')

    # Format plot
    beta_B_array = data.columns.to_numpy()
    eta_B_array = data.index.to_numpy()

    ax.set_xticks(np.arange(len(beta_B_array)))
    ax.set_yticks(np.arange(len(eta_B_array)))

    ax.set_xticklabels('${:.2f}$'.format(item) for item in beta_B_array)
    ax.set_yticklabels('${:.2f}$'.format(item) for item in eta_B_array)

    ax.set_xlabel(r'$\beta_A$')
    ax.set_ylabel(r'$\eta_A$')

    return im, max_abs_val

def plot_coexisting_simulations(ax, df, diff_theta_threshold):

    # Compute average survival difference
    grouped = df.groupby(['eta_B', 'beta_B']).mean()['diff_theta']

    data = grouped.unstack().sort_index(ascending=False)
    data_np = data.to_numpy()

    for j, row in enumerate(data_np):
        for i, elem in enumerate(row):

            if abs(elem) < diff_theta_threshold:

                ax.plot(i, j, 'ko', mfc='none')

def plot_survival_frequency(ax, df, beta_A=None, eta_A=None, cell_type='B', extent=None):

    # Compute average survival difference
    grouped = df.groupby(['eta_B', 'beta_B']).mean()['theta_{}'.format(cell_type)]

    data = grouped.unstack().sort_index(ascending=False)
    data_np = data.to_numpy()

    # Plot data
    im = ax.imshow(data_np, cmap='seismic', vmin=0, vmax=1, extent=extent)

    # Format plot
    beta_B_array = data.columns.to_numpy()
    eta_B_array = data.index.to_numpy()

    ax.set_xticks(np.arange(len(beta_B_array)))
    ax.set_yticks(np.arange(len(eta_B_array)))

    ax.set_xticklabels('${:.2f}$'.format(item) for item in beta_B_array)
    ax.set_yticklabels('${:.2f}$'.format(item) for item in eta_B_array)

    ax.set_xlabel(r'$\beta_A$')
    ax.set_ylabel(r'$\eta_A$')

    return im

def plot_ergodic_rms(ax, df, cell_type='B', extent=None):

    # Compute average survival difference
    grouped = df.groupby(['eta_B', 'beta_B']).mean()['ergodic_rms_{}'.format(cell_type)]

    data = grouped.unstack().sort_index(ascending=False)
    data_np = data.to_numpy()

    # Plot data
    im = ax.imshow(data_np, cmap='viridis', vmin=0, vmax=1, extent=extent)
    #im = ax.imshow(data_np, cmap='viridis', extent=extent)

    # Format plot
    beta_B_array = data.columns.to_numpy()
    eta_B_array = data.index.to_numpy()

    ax.set_xticks(np.arange(len(beta_B_array)))
    ax.set_yticks(np.arange(len(eta_B_array)))

    ax.set_xticklabels('${:.2f}$'.format(item) for item in beta_B_array)
    ax.set_yticklabels('${:.2f}$'.format(item) for item in eta_B_array)

    ax.set_xlabel(r'$\beta_A$')
    ax.set_ylabel(r'$\eta_A$')

    return im, data_np.min(), data_np.max()

def plot_homotypic_survival_difference(ax, df_competition, df_control, beta_A,
        eta_A, cell_type):

    grouped_control = df_control.groupby(['eta_B', 'beta_B']).mean()['theta_B']

    data_control = grouped_control.unstack().sort_index(ascending=False)
    data_np_control = data_control.to_numpy()

    grouped_competition = df_competition.groupby(['eta_B', 'beta_B']).mean()['theta_' + cell_type]

    data_competition = grouped_competition.unstack().sort_index(ascending=False)
    data_np_competition = data_competition.to_numpy()

    # Calculate homotypic survival probability
    data_np = data_np_competition - data_np_control

    # Get max_abs_val
    max_abs_val = np.absolute(data_np).max()

    # Plot data
    im = ax.imshow(data_np, cmap='seismic')

    # Format plot
    beta_B_array = data_competition.columns.to_numpy()
    eta_B_array = data_competition.index.to_numpy()

    ax.set_xticks(np.arange(len(beta_B_array)))
    ax.set_yticks(np.arange(len(eta_B_array)))

    ax.set_xticklabels('${:.2f}$'.format(item) for item in beta_B_array)
    ax.set_yticklabels('${:.2f}$'.format(item) for item in eta_B_array)

    ax.set_xlabel(r'$\beta_A$')
    ax.set_ylabel(r'$\eta_A$')

    return im, max_abs_val

def plot_neutral_simulations(ax, df_competition, df_control, cell_type,
        diff_theta_threshold):

    grouped_control = df_control.groupby(['eta_B', 'beta_B']).mean()['theta_' + cell_type]

    data_control = grouped_control.unstack().sort_index(ascending=False)
    data_np_control = data_control.to_numpy()

    grouped_competition = df_competition.groupby(['eta_B', 'beta_B']).mean()['theta_' + cell_type]

    data_competition = grouped_competition.unstack().sort_index(ascending=False)
    data_np_competition = data_competition.to_numpy()

    # Calculate homotypic survival probability
    data_np = data_np_competition - data_np_control

    for j, row in enumerate(data_np):
        for i, elem in enumerate(row):

            if abs(elem) < diff_theta_threshold:

                ax.plot(i, j, 'ks', mfc='none')

def plot_point(ax, df, beta, eta, formatstyle='go', beta_B_array=None):

    if beta_B_array is None:
        beta_B_array = df['beta_B'].unique()
    eta_B_array = df['eta_B'].unique()[::-1]

    # Define bounds of box in (beta, eta) coordinates
    beta_diff = beta_B_array[1] - beta_B_array[0]
    eta_diff = eta_B_array[-2] - eta_B_array[-1]

    beta_L = beta_B_array[0] - 0.5 * beta_diff
    beta_R = beta_B_array[-1] + 0.5 * beta_diff

    eta_B = eta_B_array[-1] - 0.5 * eta_diff
    eta_T = eta_B_array[0] + 0.5 * eta_diff

    # beta_B_array and eta_B_array are equally spaced, so we can
    # draw the dividing line
    if np.all(np.isclose(np.ediff1d(eta_B_array), np.ediff1d(eta_B_array)[0])) and \
        np.all(np.isclose(np.ediff1d(beta_B_array), np.ediff1d(beta_B_array)[0])):

        # Plot point beta, eta if within beta_B_array, eta_B_array, respectively
        if (beta_L <= beta) and (beta <= beta_R) and (eta_B <= eta) and (eta <= eta_T):

            # Convert to grid coordinates
            x = np.interp(beta,
                    [beta_L, beta_R],
                    [-0.5, len(beta_B_array) - 0.5])
            y = np.interp(eta,
                    [eta_B, eta_T],
                    [len(eta_B_array) - 0.5, -0.5])

            # Plot point
            ax.plot(x, y, formatstyle)

def plot_neutral_coexistence_point(ax, df, beta_A, eta_A):

    plot_point(ax, df, beta_A, eta_A, 'go')

def plot_coexistence_curve(ax, df, beta_A, eta_A):

    beta_B_array = df['beta_B'].unique()
    eta_B_array = df['eta_B'].unique()[::-1]

    # Define bounds of box in (beta, eta) coordinates
    beta_diff = beta_B_array[1] - beta_B_array[0]
    eta_diff = eta_B_array[-2] - eta_B_array[-1]

    beta_L = beta_B_array[0] - 0.5 * beta_diff
    beta_R = beta_B_array[-1] + 0.5 * beta_diff

    eta_B = eta_B_array[-1] - 0.5 * eta_diff
    eta_T = eta_B_array[0] + 0.5 * eta_diff

    # beta_B_array and eta_B_array are equally spaced, so we can
    # draw the dividing line
    if np.all(np.isclose(np.ediff1d(eta_B_array), np.ediff1d(eta_B_array)[0])) and \
        np.all(np.isclose(np.ediff1d(beta_B_array), np.ediff1d(beta_B_array)[0])):

        # Find candidate points as intersections between box lines
        # and dividing line
        candidate_points = np.array([
                [eta_B * beta_A / eta_A, eta_B],
                [eta_T * beta_A / eta_A, eta_T],
                [beta_L, beta_L * eta_A / beta_A],
                [beta_R, beta_R * eta_A / beta_A]])

        # Filter points if they are within the box
        points = []
        for candidate_point in candidate_points:

            if (beta_L <= candidate_point[0] or np.isclose(beta_L, candidate_point[0])) \
                    and (candidate_point[0] <= beta_R or np.isclose(candidate_point[0], beta_R)) \
                    and (eta_B <= candidate_point[1] or np.isclose(eta_B, candidate_point[1])) \
                    and (candidate_point[1] <= eta_T or np.isclose(candidate_point[1], eta_T)):
                points.append(candidate_point)

        # We need at least two candidate points (with possible duplicates)
        if not len(points) >= 2:
            assert len(points) == 0
            return im

        # Remove duplicates
        final_points = [points[0]]
        for point in points[1:]:
            if (not np.isclose(final_points[0][0], point[0])) and \
                    (not np.isclose(final_points[0][1], point[1])):
                final_points.append(point)
                break

        # Convert to grid coordinates
        x0, x1 = np.interp(np.array([final_points[0][0], final_points[1][0]]),
                [beta_L, beta_R],
                [-0.5, len(beta_B_array) - 0.5])
        y0, y1 = np.interp(np.array([final_points[0][1], final_points[1][1]]),
                [eta_B, eta_T],
                [len(eta_B_array) - 0.5, -0.5])

        # Plot dividing line
        ax.plot([x0, x1], [y0, y1], 'k-')

def plot_B_loser_viability_curve(ax, df, beta_A, eta_A):

    beta_B_array = df['beta_B'].unique()
    eta_B_array = df['eta_B'].unique()[::-1]

    # Define bounds of box in (beta, eta) coordinates
    beta_diff = beta_B_array[1] - beta_B_array[0]
    eta_diff = eta_B_array[-2] - eta_B_array[-1]

    beta_L = beta_B_array[0] - 0.5 * beta_diff
    beta_R = beta_B_array[-1] + 0.5 * beta_diff

    eta_B = eta_B_array[-1] - 0.5 * eta_diff
    eta_T = eta_B_array[0] + 0.5 * eta_diff

    # beta_B_array and eta_B_array are equally spaced, so we can
    # draw the dividing line
    if np.all(np.isclose(np.ediff1d(eta_B_array), np.ediff1d(eta_B_array)[0])) and \
        np.all(np.isclose(np.ediff1d(beta_B_array), np.ediff1d(beta_B_array)[0])):

        # Find candidate points as intersections between box lines
        # and dividing line
        candidate_points = np.array([
                [eta_B / (np.log(2) * (1 - beta_A)), eta_B],
                [eta_T / (np.log(2) * (1 - beta_A)), eta_T],
                [beta_L, beta_L * np.log(2) * (1 - beta_A)],
                [beta_R, beta_R * np.log(2) * (1 - beta_A)]])

        # Filter points if they are within the box
        points = []
        for candidate_point in candidate_points:

            if (beta_L <= candidate_point[0] or np.isclose(beta_L, candidate_point[0])) \
                    and (candidate_point[0] <= beta_R or np.isclose(candidate_point[0], beta_R)) \
                    and (eta_B <= candidate_point[1] or np.isclose(eta_B, candidate_point[1])) \
                    and (candidate_point[1] <= eta_T or np.isclose(candidate_point[1], eta_T)):
                points.append(candidate_point)

        # We need at least two candidate points (with possible duplicates)
        if not len(points) >= 2:
            assert len(points) == 0
            return im

        # Remove duplicates
        final_points = [points[0]]
        for point in points[1:]:
            if (not np.isclose(final_points[0][0], point[0])) and \
                    (not np.isclose(final_points[0][1], point[1])):
                final_points.append(point)
                break

        # Convert to grid coordinates
        x0, x1 = np.interp(np.array([final_points[0][0], final_points[1][0]]),
                [beta_L, beta_R],
                [-0.5, len(beta_B_array) - 0.5])
        y0, y1 = np.interp(np.array([final_points[0][1], final_points[1][1]]),
                [eta_B, eta_T],
                [len(eta_B_array) - 0.5, -0.5])

        # Plot dividing line
        ax.plot([x0, x1], [y0, y1], '-', color=A_loser_viability_curve_colour)

def plot_neutral_competition_curve(ax, df, beta_A, eta_A):

    beta_B_array = df['beta_B'].unique()
    eta_B_array = df['eta_B'].unique()[::-1]

    # Define bounds of box in (beta, eta) coordinates
    beta_diff = beta_B_array[1] - beta_B_array[0]
    eta_diff = eta_B_array[-2] - eta_B_array[-1]

    beta_L = beta_B_array[0] - 0.5 * beta_diff
    beta_R = beta_B_array[-1] + 0.5 * beta_diff

    eta_B = eta_B_array[-1] - 0.5 * eta_diff
    eta_T = eta_B_array[0] + 0.5 * eta_diff

    # beta_B_array and eta_B_array are equally spaced, so we can
    # draw the dividing line
    if np.all(np.isclose(np.ediff1d(eta_B_array), np.ediff1d(eta_B_array)[0])) and \
        np.all(np.isclose(np.ediff1d(beta_B_array), np.ediff1d(beta_B_array)[0])):


        # Plot point neutral competition if beta_A within beta_B_array
        if (beta_L <= beta_A) and (beta_A <= beta_R):

            # Convert to grid coordinates
            x = np.interp(beta_A,
                    [beta_L, beta_R],
                    [-0.5, len(beta_B_array) - 0.5])

            # Plot point
            ax.axvline(x, color='black', linestyle='dashed')

def plot_A_loser_viability_curve(ax, df, beta_A, eta_A):

    beta_B_array = df['beta_B'].unique()
    eta_B_array = df['eta_B'].unique()[::-1]

    # Define bounds of box in (beta, eta) coordinates
    beta_diff = beta_B_array[1] - beta_B_array[0]
    eta_diff = eta_B_array[-2] - eta_B_array[-1]

    beta_L = beta_B_array[0] - 0.5 * beta_diff
    beta_R = beta_B_array[-1] + 0.5 * beta_diff

    eta_B = eta_B_array[-1] - 0.5 * eta_diff
    eta_T = eta_B_array[0] + 0.5 * eta_diff

    beta_plot = 1 - eta_A / (np.log(2) * beta_A)
    eta_min = eta_A / beta_A * (1 - eta_A / (np.log(2) * beta_A))

    # beta_B_array and eta_B_array are equally spaced, so we can
    # draw the dividing line
    if np.all(np.isclose(np.ediff1d(eta_B_array), np.ediff1d(eta_B_array)[0])) and \
        np.all(np.isclose(np.ediff1d(beta_B_array), np.ediff1d(beta_B_array)[0])):

        # Plot point A loser viability curve if beta_B within beta_B_array
        if (beta_L <= beta_plot) and (beta_plot <= beta_R):

            # Convert to grid coordinates
            x = np.interp(beta_plot,
                    [beta_L, beta_R],
                    [-0.5, len(beta_B_array) - 0.5])

            ymin, ymax = np.interp(np.array([eta_min, eta_T]),
                    [eta_B, eta_T],
                    [len(eta_B_array) - 0.5, -0.5])

            # Plot line
            ax.vlines(x, ymin, ymax, color=B_loser_viability_curve_colour, linestyle='solid')

            ymin, ymax = np.interp(np.array([eta_B, eta_min]),
                    [eta_B, eta_T],
                    [len(eta_B_array) - 0.5, -0.5])

            # Plot point
            ax.vlines(x, ymin, ymax, color=B_loser_viability_curve_colour, linestyle='dotted')

def plot_B_winner_viability_curve(ax, df, beta_A, eta_A, N=100):

    beta_B_array = df['beta_B'].unique()
    eta_B_array = df['eta_B'].unique()[::-1]

    # Define bounds of box in (beta, eta) coordinates
    beta_diff = beta_B_array[1] - beta_B_array[0]
    eta_diff = eta_B_array[-2] - eta_B_array[-1]

    beta_L = beta_B_array[0] - 0.5 * beta_diff
    beta_R = beta_B_array[-1] + 0.5 * beta_diff

    eta_B = eta_B_array[-1] - 0.5 * eta_diff
    eta_T = eta_B_array[0] + 0.5 * eta_diff

    beta_plot_array1 = np.linspace(beta_L, 1 - eta_A / (np.log(2) * beta_A), N)
    eta_plot_array1 = np.log(2) * beta_plot_array1 * (1 - beta_plot_array1)

    beta_plot_array2 = np.linspace(1 - eta_A / (np.log(2) * beta_A), beta_R, N)
    eta_plot_array2 = np.log(2) * beta_plot_array2 * (1 - beta_plot_array2)

    # beta_B_array and eta_B_array are equally spaced, so we can
    # draw the dividing line
    if np.all(np.isclose(np.ediff1d(eta_B_array), np.ediff1d(eta_B_array)[0])) and \
        np.all(np.isclose(np.ediff1d(beta_B_array), np.ediff1d(beta_B_array)[0])):

        # Convert to grid coordinates
        x = np.interp(beta_plot_array1,
                [beta_L, beta_R],
                [-0.5, len(beta_B_array) - 0.5])

        y = np.interp(eta_plot_array1,
                [eta_B, eta_T],
                [len(eta_B_array) - 0.5, -0.5])

        # Plot point
        ax.plot(x, y, color='tab:blue', linestyle='solid')

        # Convert to grid coordinates
        x = np.interp(beta_plot_array2,
                [beta_L, beta_R],
                [-0.5, len(beta_B_array) - 0.5])

        y = np.interp(eta_plot_array2,
                [eta_B, eta_T],
                [len(eta_B_array) - 0.5, -0.5])

        # Plot point
        ax.plot(x, y, color='tab:blue', linestyle='dotted')

def plot_homotypic_viability_curve(ax, df, N=100, coef=np.log(2), offset=0,
        beta_B_array=None, color='black'):

    if beta_B_array is None:
        beta_B_array = df['beta_B'].unique()
    eta_B_array = df['eta_B'].unique()[::-1]

    # Define bounds of box in (beta, eta) coordinates
    beta_diff = beta_B_array[1] - beta_B_array[0]
    eta_diff = eta_B_array[-2] - eta_B_array[-1]

    beta_L = beta_B_array[0] - 0.5 * beta_diff
    beta_R = beta_B_array[-1] + 0.5 * beta_diff

    eta_B = eta_B_array[-1] - 0.5 * eta_diff
    eta_T = eta_B_array[0] + 0.5 * eta_diff

    beta_plot_array = np.linspace(beta_L, beta_R, N)
    eta_plot_array = coef * beta_plot_array * (1 - beta_plot_array)

    #beta_plot_array2 = np.linspace(1 - eta_A / (np.log(2) * beta_A), beta_R, N)
    #eta_plot_array2 = np.log(2) * beta_plot_array2 * (1 - beta_plot_array2)

    # beta_B_array and eta_B_array are equally spaced, so we can
    # draw the dividing line
    if np.all(np.isclose(np.ediff1d(eta_B_array), np.ediff1d(eta_B_array)[0])) and \
        np.all(np.isclose(np.ediff1d(beta_B_array), np.ediff1d(beta_B_array)[0])):

        # Convert to grid coordinates
        x = np.interp(beta_plot_array,
                [beta_L, beta_R],
                [-0.5, len(beta_B_array) - 0.5])

        y = np.interp(eta_plot_array,
                [eta_B, eta_T],
                [len(eta_B_array) - 0.5, -0.5])

        # Plot point
        ax.plot(x + offset, y, color=color, linestyle='solid')

def plot_lower_curve(ax, df, rho, N=100, offset=0, beta_B_array=None):

    if beta_B_array is None:
        beta_B_array = df['beta_B'].unique()
    eta_B_array = df['eta_B'].unique()[::-1]

    # Define bounds of box in (beta, eta) coordinates
    beta_diff = beta_B_array[1] - beta_B_array[0]
    eta_diff = eta_B_array[-2] - eta_B_array[-1]

    beta_L = beta_B_array[0] - 0.5 * beta_diff
    beta_R = beta_B_array[-1] + 0.5 * beta_diff

    eta_B = eta_B_array[-1] - 0.5 * eta_diff
    eta_T = eta_B_array[0] + 0.5 * eta_diff

    # Compute intersection of lower curve with bottom of box
    D = (1 + rho)**2 - 4 * (eta_B + rho)
    if D < 0:
        return

    beta_min = ((1 + rho) - np.sqrt(D)) / 2

    beta_plot_array = np.linspace(beta_min, beta_R, N)
    eta_plot_array = (beta_plot_array - rho) * (1 - beta_plot_array)

    # beta_B_array and eta_B_array are equally spaced, so we can
    # draw the dividing line
    if np.all(np.isclose(np.ediff1d(eta_B_array), np.ediff1d(eta_B_array)[0])) and \
        np.all(np.isclose(np.ediff1d(beta_B_array), np.ediff1d(beta_B_array)[0])):

        # Convert to grid coordinates
        x = np.interp(beta_plot_array,
                [beta_L, beta_R],
                [-0.5, len(beta_B_array) - 0.5])

        y = np.interp(eta_plot_array,
                [eta_B, eta_T],
                [len(eta_B_array) - 0.5, -0.5])

        # Plot point
        ax.plot(x + offset, y, color=lower_colour, linestyle=lower_style)

def plot_upper_curve(ax, df, rho, N=100, offset=0, beta_B_array=None):

    if beta_B_array is None:
        beta_B_array = df['beta_B'].unique()
    eta_B_array = df['eta_B'].unique()[::-1]

    # Define bounds of box in (beta, eta) coordinates
    beta_diff = beta_B_array[1] - beta_B_array[0]
    eta_diff = eta_B_array[-2] - eta_B_array[-1]

    beta_L = beta_B_array[0] - 0.5 * beta_diff
    beta_R = beta_B_array[-1] + 0.5 * beta_diff

    eta_B = eta_B_array[-1] - 0.5 * eta_diff
    eta_T = eta_B_array[0] + 0.5 * eta_diff

    beta_plot_array = np.linspace(beta_L, beta_R, N)
    eta_plot_array = (beta_plot_array + rho) * (1 - beta_plot_array)

    # beta_B_array and eta_B_array are equally spaced, so we can
    # draw the dividing line
    if np.all(np.isclose(np.ediff1d(eta_B_array), np.ediff1d(eta_B_array)[0])) and \
        np.all(np.isclose(np.ediff1d(beta_B_array), np.ediff1d(beta_B_array)[0])):

        # Convert to grid coordinates
        x = np.interp(beta_plot_array,
                [beta_L, beta_R],
                [-0.5, len(beta_B_array) - 0.5])

        y = np.interp(eta_plot_array,
                [eta_B, eta_T],
                [len(eta_B_array) - 0.5, -0.5])

        # Plot point
        ax.plot(x + offset, y, color=upper_colour, linestyle=upper_style)

def plot_homotypic_well_mixed_extinction_frequency(ax, df, include_sim_end=False, extent=None):
    """
    Plots extinction frequency against eta_B, beta_B on given ax.

    Parameters
        ax: ax to plot survival difference on
        df: pandas.DataFrame with diff_theta column
        beta_A: value of beta_A to plot predicted border
        eta_A: value of eta_A to plot predicted border
        colourbar
    """

    # Get unique status
    unique_status = df['status'].unique()

    # Compute average survival difference
    grouped = df.groupby(['eta_B', 'beta_B'])['status'].value_counts()

    data_grid = grouped.unstack().fillna(0)

    status_count_data = { status : data_grid[status].unstack().sort_index(ascending=False).to_numpy()
            for status in unique_status}

    # 0 end simulation
    # 1 zero cell count
    # 2 max cell count
    # 3 min cell count
    # 4 clone-specific max cell count
    # 5 clone-specific min cell count

    assert not 1 in status_count_data
    assert not 4 in status_count_data
    assert not 5 in status_count_data

    if not 0 in unique_status:
        status_count_data[0] = status_count_data[unique_status[0]] * 0
    if not 2 in unique_status:
        status_count_data[2] = status_count_data[unique_status[0]] * 0
    if not 3 in unique_status:
        status_count_data[3] = status_count_data[unique_status[0]] * 0

    #print('Number of simulations that reached end of simulations: {}'.format(np.sum(data_np0)))

    total = status_count_data[2] + status_count_data[3]

    if include_sim_end:
        total += status_count_data[0]

    data_np = status_count_data[0] * 0
    mask = total > 0
    data_np[mask] = status_count_data[3][mask] / total[mask]
    
    # Plot data
    im = ax.imshow(data_np, cmap='Greys', extent=extent)

    # Format plot
    grouped = df.groupby(['eta_B', 'beta_B']).mean()['theta_B']
    data = grouped.unstack().sort_index(ascending=False)
    beta_B_array = data.columns.to_numpy()
    eta_B_array = data.index.to_numpy()

    ax.set_xticks(np.arange(len(beta_B_array)))
    ax.set_yticks(np.arange(len(eta_B_array)))

    ax.set_xticklabels('${:.2f}$'.format(item) for item in beta_B_array)
    ax.set_yticklabels('${:.2f}$'.format(item) for item in eta_B_array)

    ax.set_xlabel(r'$\beta_A$')
    ax.set_ylabel(r'$\eta_A$')

    return im

def plot_homotypic_vertex_extinction_frequency(ax, df, include_sim_end=False,
        extent=None, end_time=10000, threshold=100):
    """
    Plots extinction frequency against eta_B, beta_B on given ax.

    Parameters
        ax: ax to plot survival difference on
        df: pandas.DataFrame with diff_theta column
        beta_A: value of beta_A to plot predicted border
        eta_A: value of eta_A to plot predicted border
        colourbar
    """

    # copy dataframe
    df = pd.DataFrame(df)

    # 0 end simulation
    # 1 zero cell count
    # 2 max cell count
    # 3 min cell count
    # 4 clone-specific max cell count
    # 5 clone-specific min cell count

    # Failed simulations get status -1
    df['status'] = -1

    df.loc[(df['State'] == 'COMPLETED') & (df['final_timestep'] == end_time), 'status'] = 0
    df.loc[(df['State'] == 'COMPLETED') & (df['final_timestep'] < end_time) & (df['final_cell_count'] <= threshold), 'status'] = 3
    df.loc[(df['State'] == 'COMPLETED') & (df['final_timestep'] < end_time) & (df['final_cell_count'] > threshold), 'status'] = 2

    return plot_homotypic_well_mixed_extinction_frequency(ax, df, include_sim_end, extent)

def plot_rho_curve(ax, df, rho, beta_B_array=None):

    if beta_B_array is None:
        beta_B_array = df['beta_B'].unique()
    eta_B_array = df['eta_B'].unique()[::-1]

    # Define bounds of box in (beta, eta) coordinates
    beta_diff = beta_B_array[1] - beta_B_array[0]
    eta_diff = eta_B_array[-2] - eta_B_array[-1]

    beta_L = beta_B_array[0] - 0.5 * beta_diff
    beta_R = beta_B_array[-1] + 0.5 * beta_diff

    eta_B = eta_B_array[-1] - 0.5 * eta_diff
    eta_T = eta_B_array[0] + 0.5 * eta_diff

    # beta_B_array and eta_B_array are equally spaced, so we can
    # draw the dividing line
    if np.all(np.isclose(np.ediff1d(eta_B_array), np.ediff1d(eta_B_array)[0])) and \
        np.all(np.isclose(np.ediff1d(beta_B_array), np.ediff1d(beta_B_array)[0])):


        # Convert to grid coordinates
        x = np.interp(rho,
                [beta_L, beta_R],
                [-0.5, len(beta_B_array) - 0.5])

        # Plot point
        ax.axvline(x, color=rho_colour, linestyle=rho_style)

if __name__ == '__main__':

    import pandas as pd
    import matplotlib as mpl

    # Plot customisation
    mpl.rcParams['lines.linewidth'] = 3
    mpl.rcParams['lines.markersize'] = 10
    mpl.rcParams['font.size'] = 10
    plt.rc('text', usetex=True)

    df = pd.read_csv('heterotypic-model-2-exponential-data-18-06-2021.csv')

    # Compute survival differential
    df['theta_A'] = df['num_divisions_A'] / (df['num_divisions_A'] + df['num_deaths_A'])
    df['theta_B'] = df['num_divisions_B'] / (df['num_divisions_B'] + df['num_deaths_B'])

    df['diff_theta'] = df['theta_B'] - df['theta_A']

    # Get unique eta_A and beta_A values
    unique_eta_As = df['eta_A'].unique()
    unique_beta_As = df['beta_A'].unique()

    # Get bounds of survival difference
    min_diff_theta = df['diff_theta'].min()
    max_diff_theta = df['diff_theta'].max()

    max_abs_diff_theta = max(abs(min_diff_theta), abs(max_diff_theta))

    # Create subfigures
    fig, axes = plt.subplots(len(unique_eta_As), len(unique_beta_As),
            sharex=True, sharey=True)

    # Invert vertical ordering of axes
    if isinstance(axes, np.ndarray):
        axes = axes[::-1, :]
    else:
        axes = np.array([[axes]])

    # Loop over eta_A and beta_A
    for i, eta_A in enumerate(unique_eta_As):
        for j, beta_A in enumerate(unique_beta_As):

            # Select current axes
            ax = axes[i,j]

            # Select eta_A and beta_A
            cur_df = df[(df['eta_A'] == eta_A) & (df['beta_A'] == beta_A)]

            # Set title
            ax.set_title(r'$\eta_B = {}, \beta_B = {}$'.format(eta_A, beta_A))

            im = plot_heterotypic_survival_difference(ax, cur_df, beta_A, eta_A, max_abs_diff_theta)

    fig.subplots_adjust(right=0.8)
    cbar_ax = fig.add_axes([0.85, 0.15, 0.05, 0.7])
    fig.colorbar(im, cax=cbar_ax)

    # Formatting
    fig.set_figheight(10)
    fig.set_figwidth(10)

    fig.suptitle('This is a test')

    ax = axes[0,0]

    yticks_old = ax.get_yticks()
    yticklabels_old = ax.get_yticklabels()

    ax.set_yticks(yticks_old[::2])
    ax.set_yticklabels(yticklabels_old[::2])

    xticks_old = ax.get_xticks()
    xticklabels_old = ax.get_xticklabels()

    ax.set_xticks(xticks_old[::6])
    ax.set_xticklabels(xticklabels_old[::6])

    plt.show()
