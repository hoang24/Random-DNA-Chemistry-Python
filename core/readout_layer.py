from params import input_params, time_params
from random_dna_chem import RandomDNAStrandDisplacementCircuit
from perturb_chem import RandomDNAChemPerturbationGillespy2

import torch
from torch import nn
from torch.nn import functional
import torch.optim as optim


class ReadOutLayer(nn.Module):
    '''
        Readout layer of Random DNA Strand Displacement Circuit Reservoir Computer
        Single Perceptron that takes concentrations of species inside the reservoir as inputs
        Tasks: Hamming distance learning, Long-term memory task, Short-term memory task
    '''

    def __init__(self, numIn):
        super(ReadOutLayer, self).__init__()
        self.output = nn.Linear(in_features=numIn, out_features=1, bias=True)

    def forward(self, x):
        return self.output(x)


# Train the readout layer
def train_readout(readout, results, epochs, device):

    # Define criterion and optimizer
    criterion = nn.MSELoss()
    optimizer = optim.SGD(readout.parameters(), lr=0.001)

    # Initialize a list of losses and a running loss
    losses = []
    running_loss = 0

    # Go through each epoch
    for epoch in range(epochs):
        # Go through each timestep
        for t in randomDNAChem.time_params.time_array:

            # Build the input and target
            x = torch.Tensor()
            y = torch.Tensor()

            # Zero the parameter gradients
            optimizer.zerograd()

            # Make a forward pass, calculate loss, and backpropogate
            output = readout(x)
            loss = criterion(output, y)
            loss.backward()
            optimizer.step()

            # Stats
            losses.append(loss.item())
            running_loss += loss.item()
            if t % 10000 == 0:
                print("\tepoch {}, inst {:<5}\trunning loss: {}".format(epoch, i, running_loss))
                running_loss = 0

    # Return list of losses
    return losses


# Test the readout layer
def test_readout(readout, results, device):

    # Initialize tracking variables
    total = 0
    correct = 0

    # Go through each timestep
    for t in randomDNAChem.time_params.time_array:

        # Build the input and target
        x = torch.Tensor()
        y = torch.Tensor()

        # Make a forward pass, check for accuracy
        output = readout(x)
        total += 1
        if abs(output.item() - y.item()) < 0.01:
            correct += 1

        # Stats
        if t % 10000 == 0:
            print("\tinst {:<5}\tcurrent accuracy: {:.3f}%".format(t, (correct / total) * 100))

    # Print final result
    print("\tFinal accuracy: {:.3f}%".format((correct / total) * 100))


# Load data
print("\nGenerating data...")
device = torch.device("cuda:0" if torch.cuda.is_available() else "cpu")

randomDNAChem = RandomDNAStrandDisplacementCircuit(input_params=input_params, 
                                                   time_params=time_params)
numIn = randomDNAChem.species_lookup['nS']



# Training
print('Training model: ')
readout = ReadOutLayer(numIn=numIn)
losses = train_readout(readout=readout, results=train_results, epochs=1, device=device)


# Testing
print('Testing model: ')