import matplotlib.pyplot as plt
from random_dna_chem import RandomDNAStrandDisplacementCircuit
from construct_chem import RandomDNAChemConstructionGillespy2
from params import input_params, time_params
import sys, getopt
from collections import OrderedDict


def main(argv):
    '''
        Create the random DNA chemistry by calling the RandomDNAStrandDisplacementCircuit class.
        Construc the gillespy2 random DNA chemistry by calling the RandomDNAChemConstructionGillespy2 class.
        Plot the species population vs. time.
        Call this method from Random-DNA-Chemistry-Python/
            python core/simulate_chem.py -h (display help message)
            python core/simulate_chem.py -t <number_of_trajectories> (set the number of trajectories and show the plot)
            python core/simulate_chem.py -t <number_of_trajectories> -p <plot_name> (set the number of trajectories and save the plot)
    '''
    try:
        opts, args = getopt.getopt(argv,"hp:t:",["trajectories="])
    except getopt.GetoptError:
        print('Incorrect syntax. Usage: python core/simulate_chem.py -t <number_of_trajectories> -p <plot_name>')
        sys.exit(2)

    if len(argv) <= 1:
        print('Incorrect syntax. Usage: python core/simulate_chem.py -t <number_of_trajectories> -p <plot_name>')
        sys.exit(2)

    for opt, arg in opts:
        if opt in ('-h', '--help'):
            print('Usage: python core/simulate_chem.py -t <number_of_trajectories> -p <plot_name>')
            sys.exit()
        elif opt in ('-t', '--trajectories'):
            num_trajectories = int(arg)
        elif opt in ('-p', '--plot_name'):
            plot_name = str(arg)

    randomDNAChem = RandomDNAStrandDisplacementCircuit(input_params=input_params, 
                                                       time_params=time_params)

    gillespy2_model = RandomDNAChemConstructionGillespy2(non_gillespy2_chem=randomDNAChem)

    gillespy2_results = gillespy2_model.run(number_of_trajectories=num_trajectories)

    color_array = ['#000000', '#0000FF', '#00FF00', '#00FFFF', '#000080', 
                   '#008000', '#008080', '#800000', '#800080', '#808000', 
                   '#808080', '#C0C0C0', '#FF0000', '#FF00FF', '#FFFF00',
                   '#8B0000', '#006400', '#BDB76B', '#008B8B', '#191970']

    plt.figure(figsize = (18,10))
    plt.title('Plot of all molecules over a simulation time')
    plt.xlabel('Time (s)')
    plt.ylabel('Number of molecules')

    for index in range(num_trajectories):
        trajectory = gillespy2_results[index]
        for species_index, species in enumerate(randomDNAChem.species_lookup['S']):
            species_plot = plt.plot(trajectory['time'], 
                           trajectory['{}'.format(species)], 
                           color=color_array[species_index],
                           label=species)

    handles, labels = plt.gca().get_legend_handles_labels()
    by_label = OrderedDict(zip(labels, handles))
    plt.legend(by_label.values(), by_label.keys(), loc='best')
    try:
        plot_name
    except NameError:
        plt.show()
    else:
        plt.savefig('plots/' + plot_name + '.png')

if __name__ == '__main__':
    main(sys.argv[1:])
