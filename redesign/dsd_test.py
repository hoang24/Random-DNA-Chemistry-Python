import pandas as pd
import os
import pickle
import numpy as np
from random_dsd_chem import RandomDSDChemistry
from dsd_rendering import DSDRenderProgram
from params_dsd import input_params, time_params, exp


"""
    Perform the following steps in a loop until finish
    1. (script) Generate visual DSD script based on random params in params_dsd.py
    2. (user) Paste in Visual DSD, simulate, and download SimulationResult ("{index}").csv
    3. (user) Copy SimulationResult ("{index}").csv to visualDSD/"tau_{t_hold}"/"exp{index}"/ directory.')
    4. (script) Idle until the SimulationResult ("{index}").csv is found.
    5. (script) Using the number of species at the end of the previous SimulationResult ("{index}").csv to generate the next visual DSD script (step 1)
"""

def generate_DSD(directory, chemistry, initial, final, past_result, run_index):
    '''
        directory (str): name of the directory to save .txt and .csv files in
        chemistry (class): random DNA strand displacement circuit class
        initial (float): start time of the simulation period
        final (float): end time of the simulation period
        past_result (str): name of the previous simulation result CSV file generated in Visual DSD
        run_index (int): index of the run in an experiment (number of time it has been run)
    '''
    if run_index == 0:
        initial_conditions = chemistry.species_count_lookup
        past_result = None
    else:
        # Get new initial conditions (final conditions of the previous period)
        df_dsdResult = pd.read_csv(f'visualDSD/{directory}/{past_result}.csv', engine='python')
        dsdResult = df_dsdResult.to_dict('list')
        for key, value in dsdResult.items():
            dsdResult.update({key: [value[-1]]}) # take only the last value
        try:
            del dsdResult['Time ']
        except KeyError:
            del dsdResult['Time']

        # Create initial_conditions with same format as species_count_lookup of chemistry
        conU = {}
        conL = {}
        conF = {}
        conPS = {}
        conPW = {}
        for key, value in dsdResult.items():
            if key in chemistry.species_lookup['U']:
                conU.update({key: value})
            elif key in chemistry.species_lookup['L']:
                conL.update({key: value})
            elif key in chemistry.species_lookup['F']:
                conF.update({key: value})
            elif key in chemistry.species_lookup['PS']:
                conPS.update({key: value})
            elif key in chemistry.species_lookup['PW']:
                conPW.update({key: value})
            else:
                raise Exception('Unknown species.')
        initial_conditions = {
            'conU': conU,
            'conL': conL,
            'conF': conF,
            'conPS': conPS,
            'conPW': conPW
        }

    # Render a DSD file for a particular period
    filename = f'visualDSD/{directory}/rendered_{initial}s-{final}s.txt'
    DSDRenderProgram(chemistry=chemistry,
                     initial_conditions=initial_conditions, 
                     influx_index=run_index,
                     dsd_filename=filename,
                     initial=initial, final=final)

    return filename


if __name__ == '__main__':
    directory = os.path.join('tau_{}/'.format(time_params['t_hold']), exp)
    try:
        os.makedirs(f'visualDSD/{directory}/')
    except FileExistsError:
        pass

    # Generate chemistry and save to a pickle file
    rand_dsd_chem = RandomDSDChemistry(input_params=input_params, time_params=time_params)
    with open(f'visualDSD/{directory}/chemistry.pickle', 'wb') as f:
        pickle.dump(rand_dsd_chem, f)

    # Generate influx data
    df_influx = pd.DataFrame(rand_dsd_chem.rateConst_lookup['rate_IN'])
    df_influx.to_csv(f'visualDSD/{directory}/influx.csv', index=False)

    # Loop through the time_array and generate the Visual DSD script and save DSD results
    time_array = rand_dsd_chem.time_params['time_array']
    for index, (time1, time2) in enumerate(zip(time_array[:-1], time_array[1:])):
        newDSDFile = generate_DSD(directory=directory, chemistry=rand_dsd_chem, 
                                  initial=time1, final=time2, 
                                  past_result=f'SimulationResult ({index-1})', run_index=index)
        print(f'"{newDSDFile}" file generated.')
        print(f'Step 1: Paste in Visual DSD, simulate, and download SimulationResult ({index}).csv')
        print(f'Step 2: Copy to visualDSD/{directory}/ directory.')
        first_while = True # print the waiting message once
        while not os.path.exists(f'visualDSD/{directory}/SimulationResult ({index}).csv'):
            if first_while:
                print(f'Waiting for "visualDSD/{directory}/SimulationResult ({index}).csv" file')
            first_while = False
        print(f'"visualDSD/{directory}/SimulationResult ({index}).csv" file found.\n')
