'''
    User defined parameters for random DNA strand displacement circuit.
    For testing Visual DSD program
'''

input_params = {
    # 'V': 1e-6, # Volume of the chamber containing the chemistry
    'n': 9, # Number of single strands [5, 10)
    'p': 0.846, # ratio of upper to lower strands [0.5, 1)
    'y': 0.214, # ratio of upper strands with complements [0, 1)
    'a_in': 0.222, # ratio of influces to the overall number of strands [0, 1)
    'a_out': 0.0193, # ratio of outfluces to the overall number of strands [0, 1)
    'theta': { # (positive) normal distribution of rate constants
        'mean': 0.148, # [0.05, 0.2)
        'variance': 0.00530 # [0, 0.02)
    },
    'theta_in': { # normal distribution of influx rate constants
        'mean': 0.0001, # Sm_base [0, 0.0006)
        'variance': 0.0000001
    },
    'theta_out': { # normal distribution of efflux rate constants
        'mean': 0.000152, # [0, 0.0006)
        'variance': 0.0000001
    },
    'phi': { # (positive) normal distribution of partial double strands per upper strand
        'mean': 2.48, # [0, 4) # k = 3 --> |phi-3| < 0.5
        'variance': 0.136 # [0, 0.5)
    },
}

time_params = {
    't_start': 0.00, # start simulation time
    't_end': 1.00, # end simulation time
    't_perturb': 0.01, # time to start perturbing the chemistry
    't_hold': 0.10, # input hold time (time delay between each perturbation)
}