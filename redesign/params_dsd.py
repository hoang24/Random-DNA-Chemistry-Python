'''
    User defined parameters for random DNA strand displacement circuit.
    For testing Visual DSD program
'''

input_params = {
    'n': 9, # Number of single strands [5, 10)
    'p': 0.846, # ratio of upper strand per lower strand [0.5, 1)
    's': 2, # number of switching domains
    'a_in': 1, # number of influx species
    'd': { # uniform distribution of initial species count
        'ub': 1000, # upper bound
        'lb': 0 # lower bound
    },
    'theta': { # (positive) normal distribution of rate constants
        'mean': 0.148, # [0.05, 0.2)
        'variance': 0.0000000 # [0, 0.02)
    },
    'theta_in': { # normal distribution of influx rate constants
        'mean': 0.0001, # Sm_base [0, 0.0006)
        'variance': 0.0000000
    },
    'theta_out': { # normal distribution of efflux rate constants
        'mean': 0.000152, # [0, 0.0006)
        'variance': 0.0000000
    },
}

time_params = {
    't_start': 0.00, # start simulation time
    't_end': 1.00, # end simulation time
    't_perturb': 0.01, # time to start perturbing the chemistry
    't_hold': 0.10, # input hold time (time delay between each perturbation)
}

exp = 'exp4'