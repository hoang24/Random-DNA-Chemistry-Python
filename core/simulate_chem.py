import matplotlib.pyplot as plt
from random_dna_chem import RandomDNAStrandDisplacementCircuit
from construct_chem import RandomDNAChemConstructionGillespy2
from params import input_params, time_params
import sys, getopt


def main(argv):
    try:
        opts, args = getopt.getopt(argv,"h:t:",["trajectories="])
    except getopt.GetoptError:
        print('Incorrect syntax. Usage: test.py -t <number_of_trajectories>')
        sys.exit(2)

    if len(argv) <= 1:
        print('Incorrect syntax. Usage: test.py -t <number_of_trajectories>')
        sys.exit(2)

    for opt, arg in opts:
        if opt in ('-h', '--help'):
            print('Usage: test.py -t <number_of_trajectories>')
            sys.exit()
        elif opt in ('-t', '--trajectories'):
            num_trajectories = int(arg)

    randomDNAChem = RandomDNAStrandDisplacementCircuit(input_params=input_params, 
                                                       time_params=time_params)

    gillespy2_model = RandomDNAChemConstructionGillespy2(non_gillespy2_chem=randomDNAChem)

    gillespy2_results = gillespy2_model.run(number_of_trajectories=num_trajectories)
    
    for index in range(num_trajectories):
        trajectory = gillespy2_results[index]
        plotcolor = 0x000000
        color_array = ['#000000', '#0000FF', '#00FF00', '#00FFFF', '#000080', 
                       '#008000', '#008080', '#800000', '#800080', '#808000', 
                       '#808080', '#C0C0C0', '#FF0000', '#FF00FF', '#FFFF00',
                       '#8B0000', '#006400', '#BDB76B', '#008B8B', '#191970']
        for species_index, species in enumerate(randomDNAChem.species_lookup['S']):
            species_plot = plt.plot(trajectory['time'], 
                           trajectory['{}'.format(species)], 
                           color=color_array[species_index])

    plt.xlabel('Time (s)')
    plt.ylabel('Number of molecules')
    plt.show()

if __name__ == '__main__':
    main(sys.argv[1:])
