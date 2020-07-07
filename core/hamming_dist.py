from simulate_perturbed_chem import plot_concentration, create_influx_lookup, plot_influx
from readout_utils import load_chem_data, create_trainset, create_testset, analyze_error, plot_error
from readout_layer import ReadOutLayer, train_readout, test_readout
import torch
import numpy as np
import matplotlib.pyplot as plt
from collections import OrderedDict


def hamming_dist(input_params, time_params, num_epoch, plot_chem=False, error_plot=False, plot_input=False):
    '''
        Method to calculate the Hamming distance between 2 input bitstream
        Args: 
            input_params (dict): dictionary of input parameters 
            time_params (dict): dictionary of time parameters
            num_epoch (int): number of epoch
            plot_chem (bool): whether to generate the concentration plots
            error_plot (bool): whether to generate the error plots
        Returns:
            NRMSE_per_epoch[-1] (numpy.float64): NRMSE of the final epoch 
            fitness_per_epoch[-1] (numpy.float64): fitness of the final epoch 
    '''

    # Load data from chemistry
    time_lookup, concentration_lookup, randomDNAChem = load_chem_data(input_params, time_params)


    # Plot chemistry
    if plot_chem:
        plot_concentration(time_lookup=time_lookup, concentration_lookup=concentration_lookup[0])
        plot_concentration(time_lookup=time_lookup, concentration_lookup=concentration_lookup[1])

    # Create trainset and testset
    trainset = create_trainset(concentration_lookup=concentration_lookup)
    testset = create_testset(concentration_lookup=concentration_lookup)
    if len(trainset) != len(testset):
        raise BaseException


    # Create readout layer
    readout = ReadOutLayer(numIn=randomDNAChem.species_lookup['nS'])
    device = torch.device("cuda:0" if torch.cuda.is_available() else "cpu")


    # Create lookup dictionary for TARGET of Hamming distance task
    influx_lookup = create_influx_lookup(randomDNAChem=randomDNAChem, num_time_element=1001)

    if plot_input:
        plot_influx(time_lookup, influx_lookup)

    for reaction, influx in influx_lookup.items():
        scale_factor = max(influx) / 1 # scale the influx value between 0 and 1
        for i in range(len(influx)):
            influx_lookup[reaction][i] = influx[i] / scale_factor
    influx_list = []
    for val in influx_lookup.values():
        influx_list.append(val)


    def hamming_distance(s1, s2):
        '''
            Return the Hamming distance between equal-length sequences.
            Returns:
                diff_list (list): list of the bitwise differentce between 2 sequence
                hamm_dist (int): hamming distance (sum of diff_list)
        '''
        if len(s1) != len(s2):
            raise ValueError("Undefined for sequences of unequal length.")
        else:
            diff_list = [el1 != el2 for el1, el2 in zip(s1, s2)]
            hamm_dist = sum(diff_list)/len(s1)
        return diff_list, hamm_dist


    # Training
    print('Training model: ')
    train_target, hamm_dist = hamming_distance(influx_list[0], influx_list[1])
    losses, outputs = train_readout(readout=readout, trainset=trainset, target=train_target, epochs=num_epoch, device=device)


    # Performance analysis
    RMSE_per_epoch, NRMSE_per_epoch, fitness_per_epoch, avgLoss_per_epoch = analyze_error(losses=losses, num_epoch=num_epoch)
    if error_plot:
        plot_error(type_per_epoch=RMSE_per_epoch, num_epoch=num_epoch, plot_type='RMSE')
        plot_error(type_per_epoch=NRMSE_per_epoch, num_epoch=num_epoch, plot_type='NRMSE')
        plot_error(type_per_epoch=fitness_per_epoch, num_epoch=num_epoch, plot_type='fitness')
        plot_error(type_per_epoch=avgLoss_per_epoch, num_epoch=num_epoch, plot_type='avgLoss')


    # Testing
    print('Testing model: ')
    test_target, hamm_dist = hamming_distance(influx_list[0], influx_list[1])
    final_accuracy = test_readout(readout=readout, testset=testset, target=test_target, device=device)


    # Return the NRMSE and fitness of the last epoch
    return NRMSE_per_epoch[-1], fitness_per_epoch[-1]