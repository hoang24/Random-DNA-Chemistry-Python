import numpy as np
from params_dsd import input_params, time_params


"""
    Script to summarize the NRMSE of short- and long-term memory task for Python and DSD model
    Then calculate the means and stds for both tasks of both models
"""

tau_dir = f'tau_{time_params['t_hold']}'

dsd_NRMSE_short = []
py_NRMSE_short = []
dsd_NRMSE_long = []
py_NRMSE_long = []
for i in range(1, 11):
    with open(f'visualDSD/{tau_dir}/exp{i}/short_NRMSE.txt', 'r') as f:
        info = f.read().split(',')
        dsd_NRMSE_short.append(float(info[0]))
        py_NRMSE_short.append(float(info[1]))
    with open(f'visualDSD/{tau_dir}/exp{i}/long_NRMSE.txt', 'r') as f:
        info = f.read().split(',')
        dsd_NRMSE_long.append(float(info[0]))
        py_NRMSE_long.append(float(info[1]))


print(f'DSD Model: NRMSE for Short-term Memory Task = {dsd_NRMSE_short}')
print(f'Abstract Model: NRMSE for Short-term Memory Task = {py_NRMSE_short}')
print(f'DSD Model: NRMSE for Long-term Memory Task = {dsd_NRMSE_long}')
print(f'Abstract Model: NRMSE for Long-term Memory Task = {py_NRMSE_long}')

# Short-term Memory Task
dsd_mean_short = np.mean(dsd_NRMSE_short)
dsd_std_short = np.std(dsd_NRMSE_short)
py_mean_short = np.mean(py_NRMSE_short)
py_std_short = np.std(py_NRMSE_short)
print('Short-term Memory Task:')
print(f'Abstract Model: {py_mean_short} +/- {py_std_short}')
print(f'DSD Model: {dsd_mean_short} +/- {dsd_std_short}')

# Long-term Memory Task
dsd_mean_long = np.mean(dsd_NRMSE_long)
dsd_std_long = np.std(dsd_NRMSE_long)
py_mean_long = np.mean(py_NRMSE_long)
py_std_long = np.std(py_NRMSE_long)
print('Long-term Memory Task:')
print(f'Abstract Model: {py_mean_long} +/- {py_std_long}')
print(f'DSD Model: {dsd_mean_long} +/- {dsd_std_long}')
