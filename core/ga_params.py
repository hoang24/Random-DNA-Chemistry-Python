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

def initialize_params_list():
    n = int(np.random.choice(a=range(5, 10+1))) # 0
    p = np.random.uniform(low=0.5, high=1) # 1
    y = np.random.uniform(low=0, high=1) # 2
    a_in = 2/n # hamming distance case
    # a_in = np.random.uniform(low=0, high=1) # regular case
    a_out = np.random.uniform(low=0, high=1) # 3
    theta_mean = np.random.uniform(low=0.05, high=0.2) # 4
    theta_var = np.random.uniform(low=0, high=0.02) # 5
    theta_in = np.random.uniform(low=0, high=0.0006) # 6
    theta_out = np.random.uniform(low=0, high=0.0006) # 7

    nL = int(round(n / (1 + p)))
    nU = n-nL
    min_single = np.min([nL, nU])
    nF = int(round(y * min_single))

    phi_mean = np.random.uniform(low=0, high=nL-nF) # 8
    phi_var = np.random.uniform(low=0, high=0.5) # 9

    # 10 params to optimize
    input_params_list = [n, p, y, a_out, theta_mean, theta_var, 
                         theta_in, theta_out, phi_mean, phi_var]


    return input_params_list

def get_params_dict(input_params_list):
    input_params_dict = {
        'n': input_params_list[0], # Number of single strands [5, 10)
        'p': input_params_list[1], # ratio of upper to lower strands [0.5, 1)
        'y': input_params_list[2], # ratio of upper strands with complements [0, 1)
        'a_in': 2/input_params_list[0], # ratio of influces to the overall number of strands [0, 1)
        'a_out': input_params_list[3], # ratio of outfluces to the overall number of strands [0, 1)
        'theta': { # (positive) normal distribution of rate constants
            'mean': input_params_list[4], # [0.05, 0.2)
            'variance': input_params_list[5] # [0, 0.02)
        },
        'theta_in': { # normal distribution of influx rate constants
            'mean': input_params_list[6], # [0, 0.0006)
            'variance': 0.0000001
        },
        'theta_out': { # normal distribution of efflux rate constants
            'mean': input_params_list[7], # [0, 0.0006)
            'variance': 0.0000001
        },
        'phi': { # (positive) normal distribution of partial double strands per upper strand
            'mean': input_params_list[8], # [0, 4)
            'variance': input_params_list[9] # [0, 0.5)
        },
    }
    return input_params_dict