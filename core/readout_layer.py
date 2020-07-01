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


def train_readout(readout, trainset, target, epochs, device):
    '''
        Method to train the readout layer
        Args:
            readout (class): readout layer
            trainset (dict): training set (concentration normalized from 0 to 1) over time for each species 
                {'species': [scaledConcentration (list)]}
            target (list): list of target for the readout over time 
            epochs (int): number of epochs to train the network
            device (torch.device): device to train on
        Returns:
            losses (list): list of losses, containing list of losses per epoch
                losses[epochs][loss (list - time array length)]
    '''

    # Define criterion and optimizer
    criterion = nn.MSELoss()
    optimizer = optim.SGD(readout.parameters(), lr=0.001)

    # Initialize a list of losses and a running loss
    losses = []
    # running_loss = 0

    x_matrix = [] # matrix of all concentration vectors
    for concentration in trainset.values():
        x_matrix.append(concentration)

    # Go through each epoch
    for epoch in range(epochs):
        losses_per_epoch = []
        # Go through each timestep
        for i in range(len(target)):
            # Build the vertical slice of input and target at time i
            y_hat = [target[i]] # target at a particular time i
            x_matrix_i = [] # vertical vector contains each species concentration at time i
            for species_index in range(len(trainset)):
                x_matrix_i.append(x_matrix[species_index][i])
            x = torch.Tensor(x_matrix_i)
            y = torch.Tensor(y_hat)

            # Zero the parameter gradients
            optimizer.zero_grad()

            # Make a forward pass, calculate loss, and backpropogate
            output = readout(x)
            loss = criterion(output, y)
            loss.backward()
            optimizer.step()

            # Stats
            losses_per_epoch.append(loss.item())
            # running_loss += loss.item()
            # if i % 1000 == 0: # calculate cumulative loss over 1000 timestep
                # print("\tepoch {}, inst {:<4}\trunning loss: {}".format(epoch, i, running_loss))
                # running_loss = 0

        # Print weights
        print("Training weights: ")
        for weight in readout.parameters():
            print(weight.data)

        losses.append(losses_per_epoch)

    return losses


def test_readout(readout, testset, target, device):
    '''
        Method to test the readout layer
        Args:
            readout (class): readout layer
            testset (dict): testing set (concentration normalized from 0 to 1) over time for each species 
                {'species': [scaledConcentration (list)]}
            target (list): list of target for the readout over time
            device (torch.device): device to train on
        Returns:
            final_accuracy (float): final accuracy of the readout in percentage (%)
    '''

    # Initialize tracking variables
    total = 0
    correct = 0

    x_matrix = [] # matrix of all concentration vectors
    for concentration in testset.values():
        x_matrix.append(concentration)

    # Go through each timestep
    for i in range(len(target)):
        # Build the vertical slice of input and target at time i
        y_hat = [target[i]] # target at a particular time i 
        x_matrix_i = [] # vertical vector contains each species concentration at time i
        for species_index in range(len(testset)):
            x_matrix_i.append(x_matrix[species_index][i]) 
        x = torch.Tensor(x_matrix_i)
        y = torch.Tensor(y_hat)

        # Make a forward pass, check for accuracy
        output = readout(x)
        total += 1
        if abs(output.item() - y.item()) < 0.01:
            correct += 1

        accuracy = (correct / total) * 100

        # Stats
        # if i % 1000 == 0: # calculate cumulative accuracy over the entire simulation time up to each 1000 timestep
            # print("\tinst {:<4}\tcurrent accuracy: {:.3f}%".format(i, accuracy))

    # Print final result
    final_accuracy = accuracy
    print("\tFinal accuracy: {:.3f}%".format(final_accuracy))

    return final_accuracy
