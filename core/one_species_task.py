from simulate_perturbed_chem import plot_concentration
from readout_utils import load_chem_data, create_trainset, create_testset, analyze_error, plot_error
from readout_layer import ReadOutLayer, train_readout, test_readout
import torch
import numpy as np
import matplotlib.pyplot as plt


# Load data from chemistry
time_lookup, concentration_lookup, randomDNAChem = load_chem_data()


# Plot chemistry
plot_concentration(time_lookup, concentration_lookup[0])
plot_concentration(time_lookup, concentration_lookup[1])


# Create trainset and testset
trainset = create_trainset(concentration_lookup=concentration_lookup[0])
testset = create_testset(concentration_lookup=concentration_lookup[1])
if len(trainset) != len(testset):
    raise BaseException


# Create readout layer
readout = ReadOutLayer(numIn=randomDNAChem.species_lookup['nS'])
device = torch.device("cuda:0" if torch.cuda.is_available() else "cpu")


# Training
print('Training model: ')
num_epoch = 5
losses = train_readout(readout=readout, trainset=trainset, target=trainset['U0'], epochs=num_epoch, device=device)


# Performance analysis
RMSE, NRMSE, fitness, avgLoss = analyze_error(losses=losses, num_epoch=num_epoch)
plot_error(type_per_epoch=RMSE, num_epoch=num_epoch, plot_type='RMSE')
plot_error(type_per_epoch=NRMSE, num_epoch=num_epoch, plot_type='NRMSE')
plot_error(type_per_epoch=fitness, num_epoch=num_epoch, plot_type='fitness')
plot_error(type_per_epoch=avgLoss, num_epoch=num_epoch, plot_type='avgLoss')


# Testing
print('Testing model: ')
final_accuracy = test_readout(readout=readout, testset=testset, target=testset['U0'], device=device)
