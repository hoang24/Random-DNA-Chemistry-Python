from simulate_perturbed_chem import plot_concentration, create_influx_lookup
from readout_utils import load_chem_data, create_trainset, create_testset, analyze_error, plot_error
from readout_layer import ReadOutLayer, train_readout, test_readout
import torch
import numpy as np
import matplotlib.pyplot as plt
from collections import OrderedDict


def short_term_memory(input_params, time_params, num_epoch, plot_chem=False, error_plot=False, plot_target=False, plot_output=False):
    # Load data from chemistry
    time_lookup, concentration_lookup, randomDNAChem = load_chem_data(input_params, time_params)


    # Plot chemistry
    if plot_chem:
        plot_concentration(time_lookup, concentration_lookup[0])
        plot_concentration(time_lookup, concentration_lookup[1])


    # Create trainset and testset
    trainset = create_trainset(concentration_lookup=concentration_lookup)
    testset = create_testset(concentration_lookup=concentration_lookup)
    if len(trainset) != len(testset):
        raise BaseException


    # Create readout layer
    readout = ReadOutLayer(numIn=randomDNAChem.species_lookup['nS'])
    device = torch.device("cuda:0" if torch.cuda.is_available() else "cpu")


    # Create lookup dictionary for TARGET of short-term memory task
    influx_lookup = create_influx_lookup(randomDNAChem=randomDNAChem, num_time_element=1001)
    ST_lookup = influx_lookup.copy()
    for r_in, rate_in in ST_lookup.items():
        ST_target_per_reaction = []
        for ir_index in range(2, len(time_lookup) + 1): # from index 2 to index end+1
            ST_target_per_reaction.append(rate_in[ir_index - 1] + 2*rate_in[ir_index - 2])
        ST_lookup.update({'{}'.format(r_in): ST_target_per_reaction})

    for reaction, influx in ST_lookup.items():
        scale_factor = max(influx) / 1 # scale the influx value between 0 and 1
        for i in range(len(influx)):
            ST_lookup[reaction][i] = influx[i] / scale_factor
    ST_target_list = []
    for val in ST_lookup.values():
        ST_target_list.append(val)


    # Training
    print('Training model: ')
    for species, concentration in trainset.items():
        trainset.update({'{}'.format(species): concentration[2:]})
    train_target = ST_target_list[0][:-1]
    losses, outputs = train_readout(readout=readout, trainset=trainset, target=train_target, epochs=num_epoch, device=device)


    # Target plots
    if plot_target:
        color_array = ['#000000', '#0000FF', '#00FF00', '#00FFFF', '#000080',
                       '#008000', '#008080', '#800000', '#800080', '#808000',
                       '#808080', '#C0C0C0', '#FF0000', '#FF00FF', '#FFFF00',
                       '#8B0000', '#006400', '#BDB76B', '#008B8B', '#191970']
        plt.figure(figsize = (18,10))
        plt.title('Plot of Short Term Memory Task target')
        plt.xlabel('time')
        plt.ylabel('ST memory target')
        for reaction_index, (reaction, ST_target) in enumerate(ST_lookup.items()):
            plt.plot(time_lookup[2:], ST_target[:-1], color=color_array[reaction_index], label=reaction)
        handles, labels = plt.gca().get_legend_handles_labels()
        by_label = OrderedDict(zip(labels, handles))
        plt.legend(by_label.values(), by_label.keys(), loc='best')
        plt.show()


    # Output and Target plot
    if plot_output:
        plt.rc('font', family='serif')
        plt.rc('xtick', labelsize='x-large')
        plt.rc('ytick', labelsize='x-large')
        plt.figure(figsize = (9,5))
        # plt.title('Short Term Memory Task')
        plt.xlabel('time (s)', fontsize='x-large')
        plt.ylabel('Target and Output', fontsize='x-large')
        plt.plot(time_lookup[2:], train_target, color='red')
        plt.plot(time_lookup[2:], outputs, color='black')
        plt.legend(('Target', 'Output'), fontsize='x-large')
        plt.tight_layout()
        plt.show()


    # Performance analysis
    RMSE_per_epoch, NRMSE_per_epoch, fitness_per_epoch, avgLoss_per_epoch = analyze_error(losses=losses, num_epoch=num_epoch)
    if error_plot:
        plot_error(type_per_epoch=RMSE_per_epoch, num_epoch=num_epoch, plot_type='RMSE')
        plot_error(type_per_epoch=NRMSE_per_epoch, num_epoch=num_epoch, plot_type='NRMSE')
        plot_error(type_per_epoch=fitness_per_epoch, num_epoch=num_epoch, plot_type='fitness')
        plot_error(type_per_epoch=avgLoss_per_epoch, num_epoch=num_epoch, plot_type='avgLoss')


    # Testing
    print('Testing model: ')
    for species, concentration in testset.items():
        testset.update({'{}'.format(species): concentration[2:]})
    test_target = ST_target_list[0][:-1]
    final_accuracy = test_readout(readout=readout, testset=testset, target=test_target, device=device)


    # Return the NRMSE and fitness of the last epoch
    return NRMSE_per_epoch[-1], fitness_per_epoch[-1]