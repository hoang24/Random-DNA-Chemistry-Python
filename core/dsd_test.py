from python_to_dsd import PythonToDSD
import pandas as pd
import os
from random_dna_chem import RandomDNAStrandDisplacementCircuit
from simulate_perturbed_chem import simulate_chem, create_time_concentration_lookup
from params_dsd import input_params, time_params, exp
import pickle
import numpy as np


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
        initial_conditions = chemistry.concentration_lookup
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

        # Create initial_conditions with same format as concentration_lookup of chemistry
        conU = {}
        conL = {}
        conF = {}
        conP = {}
        for key, value in dsdResult.items():
            if key in chemistry.species_lookup['U']:
                conU.update({key: value})
            elif key in chemistry.species_lookup['L']:
                conL.update({key: value})
            elif key in chemistry.species_lookup['F']:
                conF.update({key: value})
            elif key in chemistry.species_lookup['P']:
                conP.update({key: value})
            else:
                raise Exception('Unknown species.')
        initial_conditions = {
            'conU': conU,
            'conL': conL,
            'conF': conF,
            'conP': conP
        }

    # Render a DSD file for a particular period
    filename = f'visualDSD/{directory}/rendered_{initial}s-{final}s.txt'
    PythonToDSD(chemistry=chemistry,
                initial_conditions=initial_conditions, 
                influx_index=run_index,
                initial=initial, final=final, points=1001, filename=filename)

    return filename

def load_chem_data(chemistry, num_trajectories):
    '''
        Load chemistry data and return 2 copies of time array and 2 copies of concentration dictionary
        Returns:
            time_lookup (list of float): time array
            concentration_lookup (list of dict of float): one for testset, one for trainset
            randomDNAChem (class): original random DNA chemistry
    '''
    gillespy2_results = simulate_chem(num_trajectories=num_trajectories, num_time_element=1001, randomDNAChem=chemistry)

    time_list = []
    concentration_list = []
    for traj in range(num_trajectories):
        time, concentrations = create_time_concentration_lookup(traj_in_use=traj,
                                                               randomDNAChem=randomDNAChem,
                                                               gillespy2_results=gillespy2_results)
        time = np.array(time)
        time_list.append(time)
        for key in concentrations.keys():
            concentrations[key] = np.array(concentrations[key])
        concentration_list.append(concentrations)

    # Averaging over num_trajectories
    time_avg = sum(time_list) / num_trajectories
    concentration_avg = {}
    for key in concentration_list[0].keys():
        concentration_avg.update({key: []})

    for key in concentration_avg.keys():
        concentration_values = []
        for traj in range(num_trajectories):
            concentration_values.append(concentration_list[traj][key])
        concentration_avg[key] = sum(concentration_values) / num_trajectories

    return time_avg, concentration_avg


if __name__ == '__main__':
    directory = os.path.join('tau_{}/'.format(time_params['t_hold']), exp)
    try:
        os.makedirs(f'visualDSD/{directory}/')
    except FileExistsError:
        pass

    # Generate chemistry and save to a pickle file
    randomDNAChem = RandomDNAStrandDisplacementCircuit(input_params=input_params, time_params=time_params)
    with open(f'visualDSD/{directory}/chemistry.pickle', 'wb') as f:
        pickle.dump(randomDNAChem, f)

    # Save Python results
    time_lookup, concentration_lookup = load_chem_data(chemistry=randomDNAChem, num_trajectories=10)
    pyResult = {**{'Time': time_lookup}, **concentration_lookup} # append time list and concentration list together
    df_pyResult = pd.DataFrame(pyResult)
    df_pyResult.to_csv(f'visualDSD/{directory}/pyResult.csv', index=False)
    df_influx = pd.DataFrame(randomDNAChem.rateConst_lookup['rate_IN'])
    df_influx.to_csv(f'visualDSD/{directory}/influx.csv', index=False)

    # Loop through the time_array and generate the Visual DSD script and save DSD results
    time_array = randomDNAChem.time_params['time_array']
    for index, (time1, time2) in enumerate(zip(time_array[:-1], time_array[1:])):
        newDSDFile = generate_DSD(directory=directory, chemistry=randomDNAChem, 
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
