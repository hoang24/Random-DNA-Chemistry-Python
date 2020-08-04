from python_to_dsd import PythonToDSD
import pandas as pd
import os
from readout_utils import load_chem_data
from params_dsd import input_params, time_params


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
        df_dsdResults = pd.read_csv(f'visualDSD/{directory}/{past_result}.csv', engine='python')
        dsdResults = df_dsdResults.to_dict('list')
        for key, value in dsdResults.items():
            dsdResults.update({key: [value[-1]]}) # take only the last value
        try:
            del dsdResults['Time ']
        except KeyError:
            del dsdResults['Time']

        # Create initial_conditions with same format as concentration_lookup of chemistry
        conU = {}
        conL = {}
        conF = {}
        conP = {}
        for key, value in dsdResults.items():
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


if __name__ == '__main__':

    # Save Python results
    time_lookup, concentration_lookup, randomDNAChem = load_chem_data(input_params, time_params)
    pyResult = {**{'Time': time_lookup}, **concentration_lookup[0]} # append time list and concentration list together
    df_pyResult = pd.DataFrame(pyResult)
    df_pyResult.to_csv(f'visualDSD/{directory}/pyResult.csv', index=False)
    df_influx = pd.DataFrame(chemistry.rateConst_lookup['rate_IN'])
    df_influx.to_csv(f'visualDSD/{directory}/influx.csv', index=False)

    # Loop through the time_array and generate the Visual DSD script and save DSD results
    time_array = randomDNAChem.time_params['time_array']
    for index, (time1, time2) in enumerate(zip(time_array[:-1], time_array[1:])):
        newDSDFile = generate_DSD(directory='exp1', chemistry=randomDNAChem, 
                                  initial=time1, final=time2, 
                                  past_result=f'SimulationResult ({index-1})', run_index=index)

        print(f'"{newDSDFile}" file generated.')
        print(f'Step 1: Paste in Visual DSD, simulate, and download SimulationResult ({index}).csv')
        print(f'Step 2: Copy to visualDSD/{directory}/ directory.')

        first_while = True # print the waiting message once
        while not os.path.exists(f'visualDSD/{directory}/SimulationResult ({index}).csv'):
            if first_run:
                print(f'Waiting for "visualDSD/{directory}/SimulationResult ({index}).csv" file')
            first_while = False
        print(f'"visualDSD/{directory}/SimulationResult ({index}).csv" file found.\n')
