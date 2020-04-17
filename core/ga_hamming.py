from params import time_params1 # hold time = 0.1
from ga_params import initialize_params_list, get_params_dict
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
sol_per_pop = 32
mate_pool_size = sol_per_pop / 2
if (mate_pool_size % 2) != 0:
    mate_pool_size = np.floor(mate_pool_size + 1)
mate_pool_size = int(mate_pool_size)

population_list = []
population_dict = []
for sol in range(sol_per_pop):
    input_params_list = initialize_params_list()
    input_params_dict = get_params_dict(input_params_list)
    population_list.append(input_params_list)
    population_dict.append(input_params_dict)
population_list = np.array(population_list)
pop_size = population_list.shape # shape[0]: sol_per_pop (8), shape[1]: number of params (10)


def eval_hamming(num_exp, num_epoch, input_params, time_params):
    '''
        Method to evaluate Hamming distance for GA
        Args:
            num_exp (int): number of experiments
            num_epoch (int): number of epochs
            input_params (dict): dictionary of input parameters
        Returns:
            NRMSE_mean (numpy.float64): mean of NRMSE over number of experiments
            NRMSE_std (numpy.float64): standard deviation of NRMSE over number of experiments
    '''

    NRMSE_per_exp = []
    for exp in range(num_exp):
        print(colored('Exp #{}'.format(exp+1), 'white', 'on_green'))
        NRMSE, __ = hamming_dist(input_params=input_params, time_params=time_params, num_epoch=num_epoch)
        NRMSE_per_exp.append(NRMSE)
    NRMSE_mean = np.mean(NRMSE_per_exp)
    NRMSE_std = np.std(NRMSE_per_exp)
        
    return NRMSE_mean, NRMSE_std


def get_error_fitness(population_dict):
    '''
        Method to get the fitness in terms of NRMSE from Hamming distance task
        Args:
            population_dict (dict): dictionary contains chemistry input parameters
        Returns:
            NRMSE_means_per_sol (list of numpy.float64): list of NRMSE, len=sol_per_pop
    '''

    NRMSE_means_per_sol = []
    for sol_idx, sol in enumerate(population_dict):
        print(colored('Chromosome #{}: baseIn = {}, thold = {}'
            .format(sol_idx, sol['theta_in']['mean'], time_params1['t_hold']), 'white', 'on_red'))
        NRMSE_means, __ = eval_hamming(num_exp=1, num_epoch=1, input_params=sol, time_params=time_params1)
        NRMSE_means_per_sol.append(NRMSE_means)
    return NRMSE_means_per_sol

def select_mating_pool(population, error, num_parents):
    '''
        Method to select parents for mating based on the low error (high fitness)
        Args:
            population (numpy.ndarray): array of parameters to optimize
                population.shape = (number of solutions, number of params)
            error (list of numpy.float64): list of NRMSE for each solution in a population, len=sol_per_pop
            num_parents (int): number of parents - size of mate pool
        Returns:
            parents (numpy.ndarray): array of parents 
                parents.shape = (number of parents, number of params)
            best_result (numpy.ndarray): array of the best params in the current generation
                best_result.shape = (1, number of params)
    '''

    for e_idx in range(len(error)): # make nan NRMSE to inf so not get selected
        if np.isnan(error[e_idx]):
            error[e_idx] = np.inf
    parents = np.empty((num_parents, population.shape[1]))
    for parent_num in range(num_parents):
        min_error_idx = np.where(error == np.min(error))
        min_error_idx = min_error_idx[0][0]
        if parent_num is 0:
            best_result = population[min_error_idx, :]
        parents[parent_num, :] = population[min_error_idx, :]
        error[min_error_idx] = np.inf
    return parents, best_result


def crossover(parents, offspring_size):
    '''
        Method to create offspring from parents through crossover
        Args:
            parents (numpy.ndarray): array of parents 
                parents.shape = (number of parents, number of params)
            offspring_size (tuple): (number of offsprings, number of params)
        Returns
            offsprings (numpy.ndarray): array of crossovered offsprings
                offsprings.shape = (number of offsprings, number of params)
    '''

    offsprings = np.empty(offspring_size)
    # The point at which crossover takes place between two parents. Usually it is at the center.
    crossover_point = np.uint8(offspring_size[1]/2)

    for k in range(offspring_size[0]):
        # Index of the first parent to mate.
        parent1_idx = k % parents.shape[0]
        # Index of the second parent to mate.
        parent2_idx = (k+1) % parents.shape[0]
        # The new offspring will have its first half of its genes taken from the first parent.
        offsprings[k, 0:crossover_point] = parents[parent1_idx, 0:crossover_point]
        # The new offspring will have its second half of its genes taken from the second parent.
        offsprings[k, crossover_point:] = parents[parent2_idx, crossover_point:]
    return offsprings


def mutation(offsprings):
    '''
        Method to create a crossovered offsprings using mutation
        Args:
            offsprings (numpy.ndarray): array of crossovered offsprings
                offsprings.shape = (number of offsprings, number of params)
        Returns:
            offsprings (numpy.ndarray): array of mutated offsprings
                offsprings.shape = (number of offsprings, number of params)
    '''

    # Mutation changes a single gene in each offspring randomly.
    for idx in range(offsprings.shape[0]):
        # The random value to be added to the gene.
        rand_param_idx = np.random.choice(a=offsprings.shape[1])
        if rand_param_idx is 0:
            rand_val = int(np.random.choice(a=range(5, 10))) # 0
        elif rand_param_idx is 1:
            rand_val = np.random.uniform(low=0.5, high=1) # 1
        elif rand_param_idx is 2:
            rand_val = np.random.uniform(low=0, high=1) # 2
        elif rand_param_idx is 3:
            rand_val = np.random.uniform(low=0, high=1) # 3
        elif rand_param_idx is 4:
            rand_val = np.random.uniform(low=0.05, high=0.2) # 4
        elif rand_param_idx is 5:
            rand_val = np.random.uniform(low=0, high=0.02) # 5
        elif rand_param_idx is 6:
            rand_val = np.random.uniform(low=0, high=0.0006) # 6
        elif rand_param_idx is 7:
            rand_val = np.random.uniform(low=0, high=0.0006) # 7
        elif rand_param_idx is 8:
            rand_val = np.random.uniform(low=0, high=4) # 8
        elif rand_param_idx is 9:
            rand_val = np.random.uniform(low=0, high=0.5) # 9
        else:
            raise BaseException
        offsprings[idx, rand_param_idx] = rand_val
    return offsprings


num_gen = 5
for gen in range(num_gen): # 
    print(colored('Generation: {}'.format(gen), 'white', 'on_blue'))
    
    # Calculate fitness: lower NRMSE higher fitness
    NRMSE_means_per_sol = get_error_fitness(population_dict)

    # Select best parents for mating
    parents, best_result = select_mating_pool(population=population_list, error=NRMSE_means_per_sol, num_parents=mate_pool_size)
    population_list[0:parents.shape[0], :] = parents

    # Generate next generation using crossover and mutation
    p_c = 1.0 # crossover probability
    p_m = 0.5 # mutation probability
    if np.random.uniform(0,1) < p_c: # crossover
        offspring_size = (pop_size[0] - parents.shape[0], population_list.shape[1])
        offsprings = crossover(parents=parents, offspring_size=offspring_size)    
        if np.random.uniform(0,1) < p_m/p_c: # mutation
            offsprings = mutation(offsprings=offsprings)
        population_list[parents.shape[0]:, :] = offsprings

    print("Best result: ", best_result)

    # update population dictionary
    for pop_idx in range(population_list.shape[0]):
        population_dict[pop_idx] = get_params_dict(population_list[pop_idx])

# Getting the best solution after iterating finishing all generations.
# Calculate fitness for final generation
NRMSE_means_per_sol = get_error_fitness(population_dict)
# Then return the index of that solution corresponding to the best fitness.
best_match_idx = np.where(NRMSE_means_per_sol == np.min(NRMSE_means_per_sol))
best_match_idx = best_match_idx[0][0]

print(colored('Best solution: ', 'red', 'on_white'))
print('n = ', population_list[best_match_idx, 0])
print('p = ', population_list[best_match_idx, 1])
print('y = ', population_list[best_match_idx, 2])
print('a_in = ', 2/population_list[best_match_idx, 0])
print('a_out = ', population_list[best_match_idx, 3])
print('theta[mean] = ', population_list[best_match_idx, 4])
print('theta[variance] = ', population_list[best_match_idx, 5])
print('theta_in = ', population_list[best_match_idx, 6])
print('theta_out = ', population_list[best_match_idx, 7])
print('phi[mean] = ', population_list[best_match_idx, 8])
print('phi[variance] = ', population_list[best_match_idx, 9])

print("Best fitness (lowest NRMSE for Hamming distance task): ", NRMSE_means_per_sol[best_match_idx])