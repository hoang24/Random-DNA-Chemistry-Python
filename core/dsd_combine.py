import numpy as np
import pandas as pd
import pickle


directory = 'exp2'
with open(f'visualDSD/{directory}/chemistry.pickle','rb') as f:
    randDNAChem = pickle.load(f)

df = pd.read_csv(f'visualDSD/{directory}/SimulationResult (0).csv', engine='python')
dsd_results = df.to_dict('list')

for index in range(1, randDNAChem.time_params['num_perturb'] + 1):
    df = pd.read_csv(f'visualDSD/{directory}/SimulationResult ({index}).csv', engine='python')
    results = df.to_dict('list')
    for key in results.keys():
        dsd_results[key] = dsd_results[key] + results[key][1:]

df_dsdResult = pd.DataFrame(dsd_results)
df_dsdResult.to_csv(f'visualDSD/{directory}/dsdResult.csv', index=False)
