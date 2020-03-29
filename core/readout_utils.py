from params import input_params, time_params
from random_dna_chem import RandomDNAStrandDisplacementCircuit
from simulate_perturbed_chem import simulate_chem, create_time_concentration_lookup
from readout_layer import ReadOutLayer, train_readout, test_readout
import copy


def load_chem_data():
    '''
        Load chemistry data
        Returns:
            time_lookup (list of list): one for testset, one for trainset
            concentration_lookup (list of dict): one for testset, one for trainset
    '''
    print("Generating data ...")
    randomDNAChem = RandomDNAStrandDisplacementCircuit(input_params=input_params, 
                                                       time_params=time_params)

    gillespy2_results = simulate_chem(num_trajectories=2, num_time_element=1001, randomDNAChem=randomDNAChem)

    time_lookup = [[], []]
    concentration_lookup = [{}, {}]
    for traj in range(2):
        time_lookup[traj], concentration_lookup[traj] = create_time_concentration_lookup(traj_in_use=traj,
                                                                                         randomDNAChem=randomDNAChem,
                                                                                         gillespy2_results=gillespy2_results)
    return time_lookup, concentration_lookup

def create_trainset(concentration_lookup):
    trainset = copy.deepcopy(concentration_lookup[0]) # trainset is a scaled concentration lookup
    for species_name, concentration in trainset.items():
        scale_factor = max(concentration) / 1 # scale the concentration value between 0 and 1
        for i in range(len(concentration)):
            trainset[species_name][i] = concentration[i] / scale_factor
    return trainset

def create_testset(concentration_lookup):
    testset = copy.deepcopy(concentration_lookup[1]) # testset is a scaled concentration lookup
    for species_name, concentration in testset.items():
        scale_factor = max(concentration) / 1 # scale the concentration value between 0 and 1
        for i in range(len(concentration)):
            testset[species_name][i] = concentration[i] / scale_factor
    return testset