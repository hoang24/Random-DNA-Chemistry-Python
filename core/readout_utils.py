from params import input_params, time_params
from random_dna_chem import RandomDNAStrandDisplacementCircuit
from simulate_perturbed_chem import simulate_chem, create_time_concentration_lookup, create_influx_lookup
from readout_layer import ReadOutLayer, train_readout, test_readout
import copy


def load_chem_data():
    '''
        Load chemistry data and return 2 copies of time array and 2 copies of concentration dictionary
        Returns:
            time_lookup (list): time array
            concentration_lookup (list of dict): one for testset, one for trainset
            randomDNAChem (class): original random DNA chemistry
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
    if time_lookup[0] == time_lookup[1]:
        time_lookup = time_lookup[0]
    else:
        raise BaseException
    return time_lookup, concentration_lookup, randomDNAChem

def create_trainset(concentration_lookup):
    '''
        Creating the trainset
        Args:
            concentration_lookup (list of dict)
        Returns:
            trainset (dict)
    '''

    trainset = copy.deepcopy(concentration_lookup[0]) # trainset is a scaled concentration lookup
    for species_name, concentration in trainset.items():
        scale_factor = max(concentration) / 1 # scale the concentration value between 0 and 1
        for i in range(len(concentration)):
            trainset[species_name][i] = concentration[i] / scale_factor
    return trainset

def create_testset(concentration_lookup):
    '''
        Creating the trainset
        Args:
            concentration_lookup (list of dict)
        Returns:
            testset (dict)
    '''

    testset = copy.deepcopy(concentration_lookup[1]) # testset is a scaled concentration lookup
    for species_name, concentration in testset.items():
        scale_factor = max(concentration) / 1 # scale the concentration value between 0 and 1
        for i in range(len(concentration)):
            testset[species_name][i] = concentration[i] / scale_factor
    return testset

# def create_ST_target(randomDNAChem):
#     influx_lookup = create_influx_lookup(randomDNAChem)
#     ST_lookup = influx_lookup.copy()
#     for r_in, rate_in in ST_lookup.items():
#         ST_target_per_reaction = []
#         for ir_index in range(2, len(time_lookup1) + 1): # from index 2 to index end+1
#             ST_target_per_reaction.append(rate_in[ir_index - 1] + 2*rate_in[ir_index - 2])
#         ST_lookup.update({'{}'.format(r_in): ST_target_per_reaction})

#     return ST_lookup

# def create_LT_target(randomDNAChem):
#     influx_lookup = create_influx_lookup(randomDNAChem)
#     t_hold = randomDNAChem.time_params['t_hold']
#     t_hold_index = time_lookup1.index(t_hold)
#     t_hold_32_index = time_lookup1.index(t_hold*(3/2))
#     LT_lookup = influx_lookup.copy() # lookup dict for short term memory task target for all reactions
#     for r_in, rate_in in LT_lookup.items():
#         LT_target_per_reaction = []
#         for ir_index in range(t_hold_32_index, len(time_lookup1) + t_hold_index): # from index hold_time to index end + index (3/4)*hold_time
#             LT_target_per_reaction.append(rate_in[ir_index - t_hold_index] + (1/2)*rate_in[ir_index - t_hold_32_index])
#         LT_lookup.update({'{}'.format(r_in): LT_target_per_reaction})

#     return LT_lookup
