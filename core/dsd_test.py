from simulate_perturbed_chem import plot_concentration, create_influx_lookup, plot_influx
from readout_utils import load_chem_data, create_trainset, create_testset, analyze_error, plot_error
from params_dsd import input_params, time_params
from python_to_dsd import PythonToDSD
from string import Template
import pandas as pd
from datetime import datetime


timestamp = datetime.now().strftime('%m-%d-%Y-%Hh-%Mm-%Ss')

# Load data from chemistry (Python)
time_lookup, concentration_lookup, randomDNAChem = load_chem_data(input_params, time_params)
SimulationResult = {**{'Time': time_lookup}, **concentration_lookup[0]}

df_concentration = pd.DataFrame(SimulationResult)
df_concentration.to_csv(f'visualDSD/concentration_{timestamp}.csv', index=False)

# Render a DSD file for the first period
PythonToDSD(chemistry=randomDNAChem, initial=0, final=0.01, points=1001, filename=f'visualDSD/rendered_{timestamp}.txt')

# Print out the influx rate dictionary
df_influx = pd.DataFrame(randomDNAChem.rateConst_lookup['rate_IN'])
df_influx.to_csv(f'visualDSD/influx_{timestamp}.csv', index=False)
