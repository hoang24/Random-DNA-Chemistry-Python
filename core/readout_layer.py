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

    def __init__(self):
        super(ReadOutLayer, self).__init__()
        # self.convolution = nn.Conv2d(in_channels=1,
        #                              out_channels=1,
        #                              kernel_size=1)
        randomDNAChem = RandomDNAStrandDisplacementCircuit(input_params=input_params, 
                                                           time_params=time_params)
        numIn = randomDNAChem.species_lookup['nS']
        print(numIn)
        self.fc = nn.Linear(in_features=numIn,
                            out_features=1,
                            bias=True)

    def forward(self, x):
        # activation_func_output = functional.relu(input=self.convolution(x))
        # x = functional.max_pool2d(kernel_size=activation_func_output, stride=(1, 1))
        x = functional.relu(self.fc(x))
        return F.log_softmax(x, dim=1)

readout = ReadOutLayer()
print(readout)

optimizer = optim.Adam(readout.parameters(), lr=0.001)

# Training the network over epochs
# EPOCHS = 3
# for epoch in range(EPOCHS): # iterate over multiple epochs
#     for data in trainset:
#         # data is a batch of featuresets and labels
#         X, y = data
#         readout.zerograd() # zero the gradient everytime before putting data thru neural net, gradient contains loss, optimizer use gradient
#         output = readout(X.view(-1, 17)) # view() reshape image, -1 batch size, 28*28 is the dimension of the image
#         loss = functional.nll_loss(output, y) # calculate loss
#         loss.backward() # backpropagation
#         optimizer.step() # adjust the weight