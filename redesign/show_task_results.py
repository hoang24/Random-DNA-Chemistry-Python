import numpy as np
from params_dsd import input_params, time_params
import os


"""
    Script to summarize the NRMSE of short- and long-term memory task for Python and DSD model
    Then calculate the means and stds for both tasks of both models
"""

tau_dir = 'tau_{}'.format(time_params['t_hold'])

dsd_NRMSE_short = []
dsd_NRMSE_long = []
for i in range(1, 11):
    with open(f'visualDSD/{tau_dir}/exp{i}/short_NRMSE.txt', 'r') as f:
        info = f.read().split(',')
        dsd_NRMSE_short.append(float(info[0]))
    with open(f'visualDSD/{tau_dir}/exp{i}/long_NRMSE.txt', 'r') as f:
        info = f.read().split(',')
        dsd_NRMSE_long.append(float(info[0]))


print(f'DSD Model: NRMSE for Short-term Memory Task = {dsd_NRMSE_short}')
print(f'DSD Model: NRMSE for Long-term Memory Task = {dsd_NRMSE_long}')

# Short-term Memory Task
dsd_mean_short = np.mean(dsd_NRMSE_short)
dsd_std_short = np.std(dsd_NRMSE_short)
print('Short-term Memory Task:')
print(f'DSD Model: {dsd_mean_short} +/- {dsd_std_short}')

# Long-term Memory Task
dsd_mean_long = np.mean(dsd_NRMSE_long)
dsd_std_long = np.std(dsd_NRMSE_long)
print('Long-term Memory Task:')
print(f'DSD Model: {dsd_mean_long} +/- {dsd_std_long}')

with open(f'visualDSD/{tau_dir}/task_results.txt', 'w') as f:
    f.write(f'DSD Model: NRMSE for Short-term Memory Task = {dsd_NRMSE_short}\n')
    f.write(f'DSD Model: NRMSE for Long-term Memory Task = {dsd_NRMSE_long}\n')
    f.write('\n')
    f.write('Short-term Memory Task:\n')
    f.write(f'DSD Model: {dsd_mean_short} +/- {dsd_std_short}\n')
    f.write('\n')
    f.write('Long-term Memory Task:\n')
    f.write(f'DSD Model: {dsd_mean_long} +/- {dsd_std_long}\n')
    f.write('\n')
