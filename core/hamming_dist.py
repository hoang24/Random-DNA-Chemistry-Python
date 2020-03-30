from simulate_perturbed_chem import plot_concentration, create_influx_lookup, plot_influx
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


# Create lookup dictionary for TARGET of Hamming distance task
influx_lookup = create_influx_lookup(randomDNAChem=randomDNAChem, num_time_element=1001)
plot_influx(time_lookup, influx_lookup)
for reaction, influx in influx_lookup.items():
    scale_factor = max(influx) / 1 # scale the influx value between 0 and 1
    for i in range(len(influx)):
        influx_lookup[reaction][i] = influx[i] / scale_factor


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
train_target, hamm_dist = hamming_distance(influx_lookup['0 --> U0'], influx_lookup['0 --> L0'])
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
test_target, hamm_dist = hamming_distance(influx_lookup['0 --> U0'], influx_lookup['0 --> L0'])
final_accuracy = test_readout(readout=readout, testset=testset, target=test_target, device=device)

