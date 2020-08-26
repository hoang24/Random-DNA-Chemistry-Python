import readout_utils as utils
import readout_layer as rlayer
import torch
import numpy as np
import pandas as pd


class ShortTermMemoryTask():

    readout = ''
    device = ''
    target = []

    def __init__(self, directory, result_name, num_epoch):
        '''
            Short Term Memory Task class
            Args:
                directory (str): directory where DSD and Python results are
                result_name (str): name of DSD or Python results (pyResult or dsdResult)
            Attributes:
                time_lookup (list of float): time array
                concentration_lookup (dict of list of float): concentration arrays for each species
                readout (class): readout layer of the RC
                device (str): device to be run on, cuda or cpu
                target (list of float): target for short-term memory task
                trainset (dict of list of float): training set for short-term memory task
        '''
        # Load results
        self.directory = directory
        self.result_name = result_name
        self.num_epoch = num_epoch
        self.time_lookup, self.concentration_lookup = self.load_result()

        # Create readout
        self.readout, self.device = self.create_readout(num_species=len(self.concentration_lookup))

        # Create target
        self.target = self.create_target()

        # Create trainset and training
        self.trainset = utils.create_trainset(concentration_lookup=self.concentration_lookup)
        self.NRMSE = self.train()

    def load_result(self):
        df_result = pd.read_csv(f'visualDSD/{self.directory}/{self.result_name}.csv', engine='python')
        dict_result = df_result.to_dict('list')
        try:
            time_lookup = np.array(dict_result['Time '])
            concentration_lookup = {key: dict_result[key] for key in dict_result if key != 'Time '}
        except KeyError:
            time_lookup = np.array(dict_result['Time'])
            concentration_lookup = {key: dict_result[key] for key in dict_result if key != 'Time'}

        return time_lookup, concentration_lookup

    def create_readout(self, num_species):
        readout = rlayer.ReadOutLayer(numIn=num_species)
        device = torch.device("cuda:0" if torch.cuda.is_available() else "cpu")
        return readout, device

    def create_target(self):
        # Load influx and create look up dictionary for influx with the same length as the time array
        df_influx = pd.read_csv(f'visualDSD/{self.directory}/influx.csv', engine='python')
        original_influx_lookup = df_influx.to_dict('list')
        influx_lookup = original_influx_lookup.copy()
        for r_in, rate_in in original_influx_lookup.items():
            influx_rate_per_reaction = []
            for ir in rate_in:
                influx_rate_per_reaction += [ir] * (len(self.time_lookup)-1)
            influx_rate_per_reaction.append(influx_rate_per_reaction[-1])
            influx_lookup.update({f'{r_in}': influx_rate_per_reaction})

        # Make a look up dictionary for target of the Short Term Memory task
        target_lookup = influx_lookup.copy()
        for r_in, rate_in in target_lookup.items():
            target_per_reaction = []
            for ir_index in range(2, len(self.time_lookup) + 1): # from index 2 to index end+1
                target_per_reaction.append(rate_in[ir_index - 1] + 2*rate_in[ir_index - 2])
            target_lookup.update({f'{r_in}': target_per_reaction})

        # Convert dict to list for the targets
        for reaction, influx in target_lookup.items():
            scale_factor = max(influx) / 1 # scale the influx value between 0 and 1
            for i in range(len(influx)):
                target_lookup[reaction][i] = influx[i] / scale_factor
        target = []
        for val in target_lookup.values():
            target.append(val)
        return target

    def train(self):
        for species, concentration in self.trainset.items():
            self.trainset.update({f'{species}': concentration[2:]})
        losses, outputs = rlayer.train_readout(readout=self.readout, 
                                        trainset=self.trainset, 
                                        target=self.target[0][:-1], 
                                        epochs=self.num_epoch, 
                                        device=self.device)

        # Performance analysis
        _, NRMSE_per_epoch, _, _ = utils.analyze_error(losses=losses, num_epoch=self.num_epoch)
        return NRMSE_per_epoch[-1]

if __name__ == '__main__':
    directory = 'exp1'
    epochs = 10

    # Short-term Memory class for DSD and Python results
    short_term_python = ShortTermMemoryTask(directory=directory, result_name='pyResult', num_epoch=epochs)
    short_term_dsd = ShortTermMemoryTask(directory=directory, result_name='dsdResult', num_epoch=epochs)

    # Get short-term memory task errors for DSD and Python results
    print(f'NRMSE of DSD results: {short_term_dsd.NRMSE} at length {len(short_term_dsd.time_lookup)}')
    print(f'NRMSE of Python results: {short_term_python.NRMSE} at length {len(short_term_python.time_lookup)}')
