import copy
import numpy as np
import matplotlib.pyplot as plt


def create_trainset(concentration_lookup):
    '''
        Creating the trainset
        Args:
            concentration_lookup (dict of float)
        Returns:
            trainset (dict of float)
    '''

    trainset = copy.deepcopy(concentration_lookup) # trainset is a scaled concentration lookup
    for species_name, concentration in trainset.items():
        scale_factor = max(concentration) / 1 # scale the concentration value between 0 and 1
        for i in range(len(concentration)):
            trainset[species_name][i] = concentration[i] / scale_factor
    return trainset

def create_testset(concentration_lookup):
    '''
        Creating the trainset
        Args:
            concentration_lookup (dict of float)
        Returns:
            testset (dict of float)
    '''

    testset = copy.deepcopy(concentration_lookup) # testset is a scaled concentration lookup
    for species_name, concentration in testset.items():
        scale_factor = max(concentration) / 1 # scale the concentration value between 0 and 1
        for i in range(len(concentration)):
            testset[species_name][i] = concentration[i] / scale_factor
    return testset

def analyze_error(losses, num_epoch):
    '''
        Function to do performance evaluation from losses per epoch
        Args:
            losses (list of list): list containing losses per epoch 
            num_epoch (int): number of epochs
        Returns:
            RMSE_per_epoch (list of numpy.float64): list containing RMSE per epoch
            NRMSE_per_epoch (list of numpy.float64): list containing NRMSE per epoch
            fitness_per_epoch (list of numpy.float64): list containing fitness per epoch
            avgLoss_per_epoch (list of numpy.float64): list containing average loss per epoch
    '''

    losses = np.array(losses)
    RMSE_per_epoch = []
    NRMSE_per_epoch = []
    fitness_per_epoch = []
    avgLoss_per_epoch = []
    for epoch in range(num_epoch): # list index 0 (epoch 1) to list index n-1 (epoch n)
        RMSE = np.sqrt(sum(losses[epoch]**2) / len(losses[epoch]))
        RMSE_per_epoch.append(RMSE)
        NRMSE = RMSE / (max(losses[epoch]) - min(losses[epoch]))
        NRMSE_per_epoch.append(NRMSE)
        fitness = 1 / RMSE
        fitness_per_epoch.append(fitness)
        avgLoss = np.mean(losses[epoch])
        avgLoss_per_epoch.append(avgLoss)
    return RMSE_per_epoch, NRMSE_per_epoch, fitness_per_epoch, avgLoss_per_epoch

def plot_error(type_per_epoch, num_epoch, plot_type):
    '''
        Function to compute and plot the average plots per epoch
        Args:
            type_per_epoch (list): list of type of error to be plotted
            num_epoch (int): number of epochs
            plot_type (str): type of plots ('RMSE', 'NRMSE', 'fitness', 'avgLoss')
    '''

    plt.figure(figsize = (18,10))
    plt.xlabel('Epoch')
    plt.plot(range(1, num_epoch+1), type_per_epoch) # plot from epoch 1 to epoch n

    if plot_type is 'RMSE':
        plt.title('RMSE vs. epochs')
        plt.ylabel('RMSE')
        # plt.savefig('RMSE_vs_epochs' + '.eps')
    elif plot_type is 'NRMSE':
        plt.title('NRMSE vs. epochs')
        plt.ylabel('NRMSE')
        # plt.savefig('NRMSE_vs_epochs' + '.eps')
    elif plot_type is 'fitness':
        plt.title('fitness vs. epochs')
        plt.ylabel('fitness')
        # plt.savefig('fitness_vs_epochs' + '.eps')
        plt.show()
    elif plot_type is 'avgLoss':
        plt.title('Average losses vs. epochs')
        plt.ylabel('avgLoss')
        # plt.savefig('avgLoss_vs_epochs' + '.eps')
    else:
        raise ValueError("Undefined type of error to be plot.")
    plt.show()
