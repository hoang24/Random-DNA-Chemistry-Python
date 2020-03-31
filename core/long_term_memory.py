from simulate_perturbed_chem import plot_concentration, create_influx_lookup
from readout_utils import load_chem_data, create_trainset, create_testset, analyze_error, plot_error
from readout_layer import ReadOutLayer, train_readout, test_readout
import torch
import numpy as np
import matplotlib.pyplot as plt
from collections import OrderedDict


# Load data from chemistry
time_lookup, concentration_lookup, randomDNAChem = load_chem_data()


# Plot chemistry
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


# Create lookup dictionary for TARGET of long-term memory task
t_hold = randomDNAChem.time_params['t_hold']
t_hold_index = time_lookup.index(t_hold)
t_hold_32_index = time_lookup.index(t_hold*(3/2))
influx_lookup = create_influx_lookup(randomDNAChem=randomDNAChem, num_time_element=1001)
LT_lookup = influx_lookup.copy() # lookup dict for short term memory task target for all reactions
for r_in, rate_in in LT_lookup.items():
    LT_target_per_reaction = []
    for ir_index in range(t_hold_32_index, len(time_lookup) + t_hold_index): # from index hold_time to index end + index (3/4)*hold_time
        LT_target_per_reaction.append(rate_in[ir_index - t_hold_index] + (1/2)*rate_in[ir_index - t_hold_32_index])
    LT_lookup.update({'{}'.format(r_in): LT_target_per_reaction})

color_array = ['#000000', '#0000FF', '#00FF00', '#00FFFF', '#000080',
               '#008000', '#008080', '#800000', '#800080', '#808000',
               '#808080', '#C0C0C0', '#FF0000', '#FF00FF', '#FFFF00',
               '#8B0000', '#006400', '#BDB76B', '#008B8B', '#191970']
plt.figure(figsize = (18,10))
plt.title('Plot of Long Term Memory Task target')
plt.xlabel('time')
plt.ylabel('LT memory target')
for reaction_index, (reaction, LT_target) in enumerate(LT_lookup.items()):
    plt.plot(time_lookup[t_hold_32_index:], LT_target[:-t_hold_index], color=color_array[reaction_index], label=reaction)
handles, labels = plt.gca().get_legend_handles_labels()
by_label = OrderedDict(zip(labels, handles))
plt.legend(by_label.values(), by_label.keys(), loc='best')
plt.show()

for reaction, influx in LT_lookup.items():
    scale_factor = max(influx) / 1 # scale the influx value between 0 and 1
    for i in range(len(influx)):
        LT_lookup[reaction][i] = influx[i] / scale_factor


# Training
print('Training model: ')
for species, concentration in trainset.items():
    trainset.update({'{}'.format(species): concentration[t_hold_32_index:]})
train_target = LT_lookup['0 --> U0'][:-t_hold_index]
num_epoch = 5
losses = train_readout(readout=readout, trainset=trainset, target=train_target, epochs=num_epoch, device=device)


# Performance analysis
RMSE, NRMSE, fitness, avgLoss = analyze_error(losses=losses, num_epoch=num_epoch)
plot_error(type_per_epoch=RMSE, num_epoch=num_epoch, plot_type='RMSE')
plot_error(type_per_epoch=NRMSE, num_epoch=num_epoch, plot_type='NRMSE')
plot_error(type_per_epoch=fitness, num_epoch=num_epoch, plot_type='fitness')
plot_error(type_per_epoch=avgLoss, num_epoch=num_epoch, plot_type='avgLoss')


# Testing
print('Testing model: ')
for species, concentration in testset.items():
    testset.update({'{}'.format(species): concentration[t_hold_32_index:]})
test_target = LT_lookup['0 --> U0'][:-t_hold_index]
final_accuracy = test_readout(readout=readout, testset=testset, target=test_target, device=device)
