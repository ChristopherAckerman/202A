import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
matplotlib.rcParams['text.usetex'] = True


def simulate_log_linear_economy(
                                Y_ss=1.15188802402873,
                                rho=0.95,
                                beta=0.952380952380952,
                                delta=0.082,
                                mu=0,
                                sigma=0.002,
                                n_obs=1100,
                                initial_k_tilde=0,
                                initial_c_tilde=0,
                                initial_z_tilde=0,
                                my_seed=100
                            ):
    """
    Takes in parameter values.
    Returns dataframe with a simulated log-linear economy.
    All x_t are defined as log deviations from the steady state.
    """
    # initialize values for c, k, and z
    # lists are vectors with one observation for each period.
    c = [initial_c_tilde]
    k = [initial_k_tilde]
    z = [initial_z_tilde]
    c_tilde = initial_c_tilde
    k_tilde = initial_k_tilde
    z_tilde = initial_z_tilde
    t = 0
    # set seed; draw random shocks
    np.random.seed(my_seed)
    epsilon = np.random.normal(mu, sigma, n_obs)
    while t < n_obs:
        # use the analytic formulas to calculate next period values
        A = (1 + gamma) ** t
        c_t = (beta * theta * K_ss ** (theta - 1)) * (1 / A) * \
                ((theta - 1) * \
                (((theta * K_ss ** (theta - 1) + 1 - delta)/A) * k_tilde + \
                (K_ss ** (theta - 1))/A * z_tilde - \
                C_ss/(A * K_ss) * c_tilde)) + \
                rho * z_tilde + \
                c_tilde

        k_t = ((theta * K_ss ** (theta - 1) + 1 - delta)/A) * k_tilde + \
            ((K_ss ** (theta - 1))/A) * z_tilde - \
            (C_ss/(A * K_ss)) * c_tilde

        z_t = rho * z_tilde + epsilon[t]

        # store the new values in each variable's list
        c.append(c_t)
        k.append(k_t)
        z.append(z_t)

        # update the old values with the new values
        k_tilde = k_t
        c_tilde = c_t
        z_tilde = z_t

        # advance time period
        t += 1

    # send the variables to a dataframe
    economy = pd.DataFrame(index=range(n_obs + 1))
    economy['c'] = c
    economy['k'] = k
    economy['z'] = z
    return economy


def graph_log_linear_economy(
        log_linear_economy,
        filename='log-linear-simulations.pdf'
                            ):
    fig, ax = plt.subplots(1, 1)
    ax.plot(log_linear_economy.index, log_linear_economy['c'],
            label=r'$\tilde{c}$')
    ax.plot(log_linear_economy.index, log_linear_economy['k'],
            label=r'$\tilde{k}$')
    ax.plot(log_linear_economy.index, log_linear_economy['z'],
            label=r'$\tilde{z}$')
    ax.legend(frameon=False)
    if filename:
        plt.savefig(filename)


def remove_log_linearization(log_linear_economy):
    my_economy = pd.DataFrame(index=log_linear_economy.index)
    my_economy['A'] = [(1 + gamma)**t for t in my_economy.index]
    my_economy['C'] = my_economy['A'] * C_ss * np.exp(log_linear_economy['c'])
    my_economy['K'] = my_economy['A'] * K_ss * np.exp(log_linear_economy['k'])
    # add Z; multiply in front of Y
    my_economy['Z'] = np.exp(log_linear_economy['z'])
    my_economy['Y'] = my_economy['Z'] * my_economy['A'] ** (1 - theta)\
        * my_economy['K'] ** theta
    my_economy['I'] = my_economy['Y'] - my_economy['C']
    return my_economy


def graph_my_economy(my_economy, filename='my-economy-simulations.pdf'):
    fig, ax = plt.subplots(1, 1)
    ax.plot(my_economy.index, my_economy['C'],
           label=r'$C$')
    ax.plot(my_economy.index, my_economy['K'],
           label=r'$K$')
    ax.plot(my_economy.index, my_economy['Y'],
           label=r'$Y$')
    ax.plot(my_economy.index, my_economy['I'],
           label=r'$I$')
    ax.legend(frameon=False)
    if filename:
        plt.savefig(filename)


if __name__ == '__main__':
    # These parameter values used in multiple functions
    C_ss=0.988469588504255
    K_ss=1.60214152474980
    theta=0.3
    gamma=0.02
    log_linear_economy = simulate_log_linear_economy()
    graph_log_linear_economy(log_linear_economy)
    my_economy = remove_log_linearization(log_linear_economy)
    graph_my_economy(my_economy)
    my_economy.to_csv('simulated_economy.csv')
    log_linear_economy.to_csv('simulated_log_linear_economy.csv')
