from simulate_perturbed_chem import plot_concentration
from readout_utils import load_chem_data, create_trainset, create_testset
from readout_layer import ReadOutLayer, train_readout, test_readout
import torch
import numpy as np
import matplotlib.pyplot as plt


# Load data from chemistry
time_lookup, concentration_lookup = load_chem_data()

# Plot chemistry
plot_concentration(time_lookup[0], concentration_lookup[0])
plot_concentration(time_lookup[1], concentration_lookup[1])

# Create trainset and testset
trainset = create_trainset(concentration_lookup=concentration_lookup)
testset = create_testset(concentration_lookup=concentration_lookup)
if len(trainset) != len(testset):
    raise BaseException

# Create readout layer
readout = ReadOutLayer(numIn=len(trainset))
device = torch.device("cuda:0" if torch.cuda.is_available() else "cpu")

# Training
print('Training model: ')
num_epoch = 5
losses = train_readout(readout=readout, trainset=trainset, target=trainset['U0'], epochs=num_epoch, device=device)

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
final_accuracy = test_readout(readout=readout, testset=testset, target=testset['U0'], device=device)
