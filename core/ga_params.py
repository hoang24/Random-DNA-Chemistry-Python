import numpy as np

'''
    Genetic algorithm script to find the fittest random DNA chemistry (lowest NRMSE over 10 experiements)
    Chromosome is a list encoding
    'n': number of single strands [5, 10)
    'p': ratio of upper to lower strands [0.5, 1)
    'y': ratio of upper strands with complements [0, 1)
    'a_in': ratio of influces to the overall number of strands [0, 1)
    'a_out': ratio of outfluces to the overall number of strands [0, 1)
    'theta': (positive) normal distribution of rate constants
        'mean': [0.05, 0.2)
        'variance': [0, 0.02)
    'theta_in': normal distribution of influx rate constants
        'mean': [0, 0.0006)
        'variance': 0.0000001
    'theta_out': normal distribution of efflux rate constants
        'mean': [0, 0.0006)
        'variance': 0.0000001
    'phi': (positive) normal distribution of partial double strands per upper strand
        'mean': [0, 4)
        'variance': [0, 0.5)

    [n, p, y, a_in, a_out, theta_mean, theta_var, theta_in, theta_out, phi_mean, phi_var]
    11-gene chromosome
'''

def get_input_params():
    n = int(np.round(np.random.uniform(low=5, high=10)))
    p = np.random.uniform(low=0.5, high=1)
    y = np.random.uniform(low=0, high=1)
    a_in = 2/n # hamming distance case
    # a_in = np.random.uniform(low=0, high=1) # regular case
    a_out = np.random.uniform(low=0, high=1)
    theta_mean = np.random.uniform(low=0.05, high=0.2)
    theta_var = np.random.uniform(low=0, high=0.02)
    theta_in = np.random.uniform(low=0, high=0.0006)
    theta_out = np.random.uniform(low=0, high=0.0006)
    phi_mean = np.random.uniform(low=0, high=4)
    phi_var = np.random.uniform(low=0, high=0.5)

    input_params_list = [n, p, y, a_in, a_out, theta_mean, theta_var, 
                         theta_in, theta_out, phi_mean, phi_var]

    input_params = {
        'n': n, # Number of single strands [5, 10)
        'p': p, # ratio of upper to lower strands [0.5, 1)
        'y': y, # ratio of upper strands with complements [0, 1)
        'a_in': a_in, # ratio of influces to the overall number of strands [0, 1)
        'a_out': a_out, # ratio of outfluces to the overall number of strands [0, 1)
        'theta': { # (positive) normal distribution of rate constants
            'mean': theta_mean, # [0.05, 0.2)
            'variance': theta_var # [0, 0.02)
        },
        'theta_in': { # normal distribution of influx rate constants
            'mean': theta_in, # [0, 0.0006)
            'variance': 0.0000001
        },
        'theta_out': { # normal distribution of efflux rate constants
            'mean': theta_out, # [0, 0.0006)
            'variance': 0.0000001
        },
        'phi': { # (positive) normal distribution of partial double strands per upper strand
            'mean': phi_mean, # [0, 4)
            'variance': phi_var # [0, 0.5)
        },
    }
    return input_params, input_params_list
