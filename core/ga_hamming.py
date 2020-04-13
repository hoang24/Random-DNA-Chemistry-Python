from params import time_params1 # hold time = 0.1
from ga_params import get_input_params
from hamming_dist import hamming_dist
import numpy as np
import colorama
from termcolor import colored


colorama.init()

"""
Genetic algorithm parameters:
    Mating pool size (number of parents mating)
    Population size (solution per population)
"""
mate_pool_size = 4
sol_per_pop = 8
population_list = []
population_dict = []
for sol in range(sol_per_pop):
    input_params, input_params_list = get_input_params()
    population_dict.append(input_params)
    population_list.append(input_params_list)
population_list = np.array(population_list)
pop_size = population_list.shape

def eval_hamming(num_exp, num_epoch, input_params):
    print(colored('baseIn = {}; thold = {}'
        .format(input_params['theta_in']['mean'], time_params1['t_hold']), 'white', 'on_red'))
    NRMSE_means = []
    NRMSE_stds = []
    for exp in range(num_exp):
        print(colored('Exp #{}'.format(exp+1), 'white', 'on_green'))
        NRMSE, __ = hamming_dist(input_params=input_params, time_params=time_params1, num_epoch=num_epoch)
        NRMSE_means.append(np.mean(NRMSE))
        NRMSE_stds.append(np.std(NRMSE))
    return NRMSE_means, NRMSE_stds

num_gen = 10
for gen in range(num_gen):
    print(colored('Generation: {}'.format(gen), 'white', 'on_blue'))
    NRMSE_means_per_sol = []
    for sol in population_dict:
        NRMSE_means, __ = eval_hamming(num_exp=1, num_epoch=1, input_params=sol)
        NRMSE_means_per_sol.append(NRMSE_means)
    print(NRMSE_means_per_sol)


