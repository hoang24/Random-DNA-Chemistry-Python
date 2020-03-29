from perturb_chem import RandomDNAChemPerturbationGillespy2
import numpy as np
import matplotlib.pyplot as plt
from collections import OrderedDict


def simulate_chem (num_trajectories, num_time_element, randomDNAChem):
    '''
        Method to simulate the Gillespy2 model of the perturbed Random DNA Chemistry
        Args: 
            num_trajectories (int): number of Gillespy2 trajectories
            num_time_element (int): number of time element in each Gillespy2 period
        Returns:
            gillespy2_results (numpy.ndarray): results from Gillespy2 stochastic simulation in form
                [period][trajectory]['time'/'species_name']
    '''

    gillespy2_results = []

    # Creating the Gillespy2 chemistry model in the non-perturb period
    gillespy2_model = RandomDNAChemPerturbationGillespy2(non_gillespy2_chem=randomDNAChem,
                                                         rate_in_timeIndex=0, # index 0 for time t=0
                                                         period_start=randomDNAChem.time_params['time_array'][0],
                                                         period_end=randomDNAChem.time_params['time_array'][1],
                                                         numel=num_time_element,
                                                         previous_gillespy2_result=None)
    gillespy2_result = gillespy2_model.run(number_of_trajectories=num_trajectories)
    gillespy2_results.append(gillespy2_result)

    # Creating the Gillespy2 chemistry models in the perturb period
    for time_index in range(1, len(randomDNAChem.time_params['time_array']) - 1):
        time_offset = randomDNAChem.time_params['t_perturb'] + randomDNAChem.time_params['t_hold'] * (time_index - 1)
        trajectories = []
        for index in range(num_trajectories):
            previous_trajectory = gillespy2_result[index] # take the gillespy2_result from the previous period
            # Creating one Gillespy2 chemistry model in one trajectory in one perturb period
            gillespy2_model = RandomDNAChemPerturbationGillespy2(non_gillespy2_chem=randomDNAChem,
                                                                 rate_in_timeIndex=time_index, # index 1 for time t=1 and so on
                                                                 period_start=randomDNAChem.time_params['time_array'][time_index] - time_offset,
                                                                 period_end=randomDNAChem.time_params['time_array'][time_index + 1] - time_offset,
                                                                 numel=num_time_element,
                                                                 previous_gillespy2_result=previous_trajectory)
            trajectory = gillespy2_model.run(number_of_trajectories=1) # run the single trajectory
            trajectories.append(trajectory)
        gillespy2_result = trajectories # gillespy2_result of the current period
        gillespy2_results.append(gillespy2_result)

    return gillespy2_results

def create_time_concentration_lookup (traj_in_use, randomDNAChem, gillespy2_results):
    '''
        Method to create lookup dictionary for a particular Gillespy2 trajectory
        Args:
            traj_in_use (int): index number of the trajectory to create lookup dict for
            randomDNAChem (class): random DNA chemistry class
            gillespy2_results (numpy.ndarray): results from Gillespy2 stochastic simulation in form
                [period][trajectory]['time'/'species_name']
        Returns:
            time_lookup (list): time array from 0 to end of simulation time
            concentration_lookup (dict): dictionary of concentration over time for each species 
                {'species': [concentration (list)]}
    '''

    time_lookup = list(gillespy2_results[0][traj_in_use]['time']) # time vector at non-perturbed period
    for time_index in range(1, randomDNAChem.time_params['num_perturb'] + 1): # time vector at perturbed period
        time_offset = randomDNAChem.time_params['t_perturb'] + randomDNAChem.time_params['t_hold'] * (time_index - 1)
        time_per_period = np.delete(gillespy2_results[time_index][traj_in_use]['time'], 0) # delete the 0th element since repeat
        time_lookup += list(time_per_period + time_offset)

    concentration_lookup = {}
    for species_name in randomDNAChem.species_lookup['S']:
        concentrations = list(gillespy2_results[0][traj_in_use][species_name]) # concentration vector at non-perturbed period
        for time_index in range(1, randomDNAChem.time_params['num_perturb'] + 1): # concentration vector at perturbed period
            concentration_per_period = np.delete(gillespy2_results[time_index][traj_in_use][species_name], 0) # delete the 0th element since repeat
            concentrations += list(concentration_per_period)
        concentration_lookup.update({'{}'.format(species_name): concentrations})

    return time_lookup, concentration_lookup

def plot_concentration (time_lookup, concentration_lookup):
    '''
        Method to plot concentration over time for a particular Gillespy2 trajectory
        Args:
            time_lookup (list): time array from 0 to end of simulation time
            concentration_lookup (dict): dictionary of concentration over time for each species 
                {'species': [concentration (list)]}
    '''

    color_array = ['#000000', '#0000FF', '#00FF00', '#00FFFF', '#000080',
                   '#008000', '#008080', '#800000', '#800080', '#808000',
                   '#808080', '#C0C0C0', '#FF0000', '#FF00FF', '#FFFF00',
                   '#8B0000', '#006400', '#BDB76B', '#008B8B', '#191970']

    plt.figure(figsize = (18,10))
    plt.title('Concentration of all species')
    plt.xlabel('time')
    plt.ylabel('concentration')

    for species_index, (species, concentration) in enumerate(concentration_lookup.items()):
        plt.plot(time_lookup, concentration, color=color_array[species_index], label=species)

    handles, labels = plt.gca().get_legend_handles_labels()
    by_label = OrderedDict(zip(labels, handles))
    plt.legend(by_label.values(), by_label.keys(), loc='best')

    # plot_name = 'random_dna_chem_plot_test'
    try:
        plot_name
    except NameError:
        plt.show()
    else:
        plt.savefig('plots/' + plot_name + '.eps')

def create_influx_lookup (randomDNAChem):
    '''
        Method to make influx rate array the same length as time array and create a lookup for influx rate of all species
        Args:
            randomDNAChem (class): random DNA chemistry class
            time_lookup (list): time array from 0 to end of simulation time

        Returns:
            influx_lookup (dict): dictionary for influx rate of each species with length of time_lookup list
                {'species': [influx_rate (list)]}
    '''

    influx_lookup = randomDNAChem.rateConst_lookup['rate_IN'].copy()
    for r_in, rate_in in randomDNAChem.rateConst_lookup['rate_IN'].items():
        influx_rate_per_reaction = []
        for ir in rate_in:
            influx_rate_per_reaction += [ir] * (num_time_element-1)
        influx_rate_per_reaction.append(influx_rate_per_reaction[-1])
        influx_lookup.update({'{}'.format(r_in): influx_rate_per_reaction})

    return influx_lookup

def plot_influx (time_lookup, influx_lookup):
    '''
        Method to plot influx rate over time
        Args:
            time_lookup (list): time array from 0 to end of simulation time
            concentration_lookup (dict): dictionary of concentration over time for each species 
                {'species': [concentration (list)]}
    '''

    color_array = ['#000000', '#0000FF', '#00FF00', '#00FFFF', '#000080',
                   '#008000', '#008080', '#800000', '#800080', '#808000',
                   '#808080', '#C0C0C0', '#FF0000', '#FF00FF', '#FFFF00',
                   '#8B0000', '#006400', '#BDB76B', '#008B8B', '#191970']
    plt.figure(figsize = (18,10))
    plt.title('Plot of influx rates')
    plt.xlabel('time')
    plt.ylabel('influx rate')
    for reaction_index, (reaction, influx_rate) in enumerate(influx_lookup.items()):
        plt.plot(time_lookup, influx_rate, color=color_array[reaction_index], label=reaction)
    handles, labels = plt.gca().get_legend_handles_labels()
    by_label = OrderedDict(zip(labels, handles))
    plt.legend(by_label.values(), by_label.keys(), loc='best')

    # plot_name = 'influx_rate_plot_test'
    try:
        plot_name
    except NameError:
        plt.show()
    else:
        plt.savefig('plots/' + plot_name + '.eps')

