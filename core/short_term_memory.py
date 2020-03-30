from simulate_perturbed_chem import plot_concentration, create_influx_lookup
from readout_utils import load_chem_data, create_trainset, create_testset
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


# Create lookup dictionary for TARGET of short-term memory task
influx_lookup = create_influx_lookup(randomDNAChem=randomDNAChem, num_time_element=1001)
ST_lookup = influx_lookup.copy()
for r_in, rate_in in ST_lookup.items():
    ST_target_per_reaction = []
    for ir_index in range(2, len(time_lookup) + 1): # from index 2 to index end+1
        ST_target_per_reaction.append(rate_in[ir_index - 1] + 2*rate_in[ir_index - 2])
    ST_lookup.update({'{}'.format(r_in): ST_target_per_reaction})

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


# Training
print('Training model: ')
for species, concentration in trainset.items():
    trainset.update({'{}'.format(species): concentration[2:]})
train_target = ST_lookup['0 --> U0'][:-1]
num_epoch = 5
losses = train_readout(readout=readout, trainset=trainset, target=train_target, epochs=num_epoch, device=device)
avg_losses_per_epoch = []
for i in range(num_epoch): # list index 0 (epoch 1) to list index 4 (epoch 5)
    avg_losses_per_epoch.append(np.mean(losses[i]))
plt.figure(figsize = (18,10))
plt.title('Plot of losses vs. epochs')
plt.xlabel('Epoch')
plt.ylabel('Loss')
plt.plot(range(1, num_epoch+1), avg_losses_per_epoch)
# plt.savefig('losses_vs_epochs' + '.eps')
plt.show()


# Testing
print('Testing model: ')
for species, concentration in testset.items():
    testset.update({'{}'.format(species): concentration[2:]})
test_target = ST_lookup['0 --> U0'][:-1]
final_accuracy = test_readout(readout=readout, testset=testset, target=test_target, device=device)

