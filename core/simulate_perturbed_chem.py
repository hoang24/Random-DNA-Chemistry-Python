from params import input_params, time_params
from random_dna_chem import RandomDNAStrandDisplacementCircuit
from perturb_chem import RandomDNAChemPerturbationGillespy2
from collections import OrderedDict
import matplotlib.pyplot as plt
import sys, getopt


randomDNAChem = RandomDNAStrandDisplacementCircuit(input_params=input_params, 
                                                   time_params=time_params)

color_array = ['#000000', '#0000FF', '#00FF00', '#00FFFF', '#000080',
               '#008000', '#008080', '#800000', '#800080', '#808000',
               '#808080', '#C0C0C0', '#FF0000', '#FF00FF', '#FFFF00',
               '#8B0000', '#006400', '#BDB76B', '#008B8B', '#191970']

plt.figure(figsize = (18,10))
plt.title('Plot of all molecules over a simulation time')
plt.xlabel('Time (s)')
plt.ylabel('Number of molecules')

gillespy2_results = []
num_trajectories = 5

# Creating the Gillespy2 chemistry model in the non-perturb period
gillespy2_model = RandomDNAChemPerturbationGillespy2(non_gillespy2_chem=randomDNAChem,
                                                     rate_in_timeIndex=0, # index 0 for time t=0
                                                     period_start=randomDNAChem.time_params['time_array'][0],
                                                     period_end=randomDNAChem.time_params['time_array'][1],
                                                     previous_gillespy2_result=None)

# Result of stochastic Gillespie simulation for non-perturbation period
gillespy2_result = gillespy2_model.run(number_of_trajectories=num_trajectories)

# Result of stochastic Gillespie simulation for entire simulation time
gillespy2_results.append(gillespy2_result)

# Plot non-perturb period
for index in range(num_trajectories):
    trajectory = gillespy2_result[index]
    for species_index, species in enumerate(randomDNAChem.species_lookup['S']):
        species_plot = plt.plot(trajectory['time'],
                                trajectory['{}'.format(species)],
                                color=color_array[species_index],
                                label=species)

# Creating the Gillespy2 chemistry models in the perturb period
for time_index in range(1, len(randomDNAChem.time_params['time_array']) - 1):
    # Calculate the Gillespy2 time offset
    time_offset = randomDNAChem.time_params['t_perturb'] + randomDNAChem.time_params['t_hold'] * (time_index - 1)

    trajectories = []
    for index in range(num_trajectories):
        previous_trajectory = gillespy2_result[index] # take the gillespy2_result from the previous period

        # Creating one Gillespy2 chemistry model in one trajectory in one perturb period
        gillespy2_model = RandomDNAChemPerturbationGillespy2(non_gillespy2_chem=randomDNAChem,
                                                             rate_in_timeIndex=time_index, # index 1 for time t=1 and so on
                                                             period_start=randomDNAChem.time_params['time_array'][time_index] - time_offset,
                                                             period_end=randomDNAChem.time_params['time_array'][time_index + 1] - time_offset,
                                                             previous_gillespy2_result=previous_trajectory)
        trajectory = gillespy2_model.run(number_of_trajectories=1) # run the single trajectory
        trajectories.append(trajectory)
    gillespy2_result = trajectories # gillespy2_result of the current period

    # Result of stochastic Gillespie simulation for entire simulation time
    gillespy2_results.append(gillespy2_result)

    # Plot each of the perturb period
    for index in range(num_trajectories):
        trajectory = gillespy2_result[index]
        for species_index, species in enumerate(randomDNAChem.species_lookup['S']):
            species_plot = plt.plot(trajectory['time'] + time_offset,
                                    trajectory['{}'.format(species)],
                                    color=color_array[species_index],
                                    label=species)

handles, labels = plt.gca().get_legend_handles_labels()
by_label = OrderedDict(zip(labels, handles))
plt.legend(by_label.values(), by_label.keys(), loc='best')
# plot_name = 'random_dna_chem_with_perturb_junk'
try:
    plot_name
except NameError:
    plt.show()
else:
    plt.savefig('plots/' + plot_name + '.eps')

# print(gillespy2_results)