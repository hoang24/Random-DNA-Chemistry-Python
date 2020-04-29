'''
    Optimized parameters from GA
    Use for short-term and long-term memory tasks
'''
n = 9
p = 0.846
y = 0.214
a_in = 1/9; # 0.222
a_out = 0.613
theta_mean = 0.148
theta_var = 0.0185
theta_in = 0.000315
theta_out = 0.000228
phi_mean = 2.57
phi_var = 0.0973

input_params1 = {
    # 'V': 1e-6, # Volume of the chamber containing the chemistry
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
        'mean': 0.0001, # Sm_base [0, 0.0006)
        'variance': 0.0000001
    },
    'theta_out': { # normal distribution of efflux rate constants
        'mean': theta_out, # [0, 0.0006)
        'variance': 0.0000001
    },
    'phi': { # (positive) normal distribution of partial double strands per upper strand
        'mean': phi_mean, # [0, 4) # k = 3 --> |phi-3| < 0.5
        'variance': phi_var # [0, 0.5)
    },
}

time_params1 = {
    't_start': 0.00, # start simulation time
    't_end': 1.00, # end simulation time
    't_perturb': 0.01, # time to start perturbing the chemistry
    't_hold': 0.10, # input hold time (time delay between each perturbation)
}


input_params2 = {
    # 'V': 1e-6, # Volume of the chamber containing the chemistry
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
        'mean': 0.0002, # Sm_base [0, 0.0006)
        'variance': 0.0000001
    },
    'theta_out': { # normal distribution of efflux rate constants
        'mean': theta_out, # [0, 0.0006)
        'variance': 0.0000001
    },
    'phi': { # (positive) normal distribution of partial double strands per upper strand
        'mean': phi_mean, # [0, 4) # k = 3 --> |phi-3| < 0.5
        'variance': phi_var # [0, 0.5)
    },
}

time_params2 = {
    't_start': 0.00, # start simulation time
    't_end': 1.00, # end simulation time
    't_perturb': 0.01, # time to start perturbing the chemistry
    't_hold': 0.20, # input hold time (time delay between each perturbation)
}


input_params3 = {
    # 'V': 1e-6, # Volume of the chamber containing the chemistry
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
        'mean': 0.0003, # Sm_base [0, 0.0006)
        'variance': 0.0000001
    },
    'theta_out': { # normal distribution of efflux rate constants
        'mean': theta_out, # [0, 0.0006)
        'variance': 0.0000001
    },
    'phi': { # (positive) normal distribution of partial double strands per upper strand
        'mean': phi_mean, # [0, 4) # k = 3 --> |phi-3| < 0.5
        'variance': phi_var # [0, 0.5)
    },
}

time_params3 = {
    't_start': 0.00, # start simulation time
    't_end': 1.00, # end simulation time
    't_perturb': 0.01, # time to start perturbing the chemistry
    't_hold': 0.30, # input hold time (time delay between each perturbation)
}


input_params4 = {
    # 'V': 1e-6, # Volume of the chamber containing the chemistry
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
        'mean': 0.0004, # Sm_base [0, 0.0006)
        'variance': 0.0000001
    },
    'theta_out': { # normal distribution of efflux rate constants
        'mean': theta_out, # [0, 0.0006)
        'variance': 0.0000001
    },
    'phi': { # (positive) normal distribution of partial double strands per upper strand
        'mean': phi_mean, # [0, 4) # k = 3 --> |phi-3| < 0.5
        'variance': phi_var # [0, 0.5)
    },
}

time_params4 = {
    't_start': 0.00, # start simulation time
    't_end': 1.00, # end simulation time
    't_perturb': 0.01, # time to start perturbing the chemistry
    't_hold': 0.40, # input hold time (time delay between each perturbation)
}


input_params5 = {
    # 'V': 1e-6, # Volume of the chamber containing the chemistry
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
        'mean': 0.0005, # Sm_base [0, 0.0006)
        'variance': 0.0000001
    },
    'theta_out': { # normal distribution of efflux rate constants
        'mean': theta_out, # [0, 0.0006)
        'variance': 0.0000001
    },
    'phi': { # (positive) normal distribution of partial double strands per upper strand
        'mean': phi_mean, # [0, 4) # k = 3 --> |phi-3| < 0.5
        'variance': phi_var # [0, 0.5)
    },
}

time_params5 = {
    't_start': 0.00, # start simulation time
    't_end': 1.00, # end simulation time
    't_perturb': 0.01, # time to start perturbing the chemistry
    't_hold': 0.50, # input hold time (time delay between each perturbation)
}


input_params6 = {
    # 'V': 1e-6, # Volume of the chamber containing the chemistry
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
        'mean': 0.0006, # Sm_base [0, 0.0006)
        'variance': 0.0000001
    },
    'theta_out': { # normal distribution of efflux rate constants
        'mean': theta_out, # [0, 0.0006)
        'variance': 0.0000001
    },
    'phi': { # (positive) normal distribution of partial double strands per upper strand
        'mean': phi_mean, # [0, 4) # k = 3 --> |phi-3| < 0.5
        'variance': phi_var # [0, 0.5)
    },
}

time_params6 = {
    't_start': 0.00, # start simulation time
    't_end': 1.00, # end simulation time
    't_perturb': 0.01, # time to start perturbing the chemistry
    't_hold': 0.60, # input hold time (time delay between each perturbation)
}

