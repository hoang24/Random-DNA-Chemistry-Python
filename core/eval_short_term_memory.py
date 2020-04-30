from time_series_params import input_params1, time_params1
from time_series_params import input_params2, time_params2
from time_series_params import input_params3, time_params3
from time_series_params import input_params4, time_params4
from time_series_params import input_params5, time_params5
from time_series_params import input_params6, time_params6
from short_term_memory import short_term_memory
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import colorama
from termcolor import colored


colorama.init()
input_params = [input_params1, input_params2, input_params3, input_params4, input_params5, input_params6]
time_params = [time_params1, time_params2, time_params3, time_params4, time_params5, time_params6]

NRMSE_means = []
NRMSE_stds = []
fitness_means = []
fitness_stds = []
base_influx_list = []
hold_time_list = []

num_exp = 20
num_epoch = 10
for input_param in input_params:
    for time_param in time_params:
        print(colored('baseIn = {}; thold = {}'.format(input_param['theta_in']['mean'], time_param['t_hold']), 'white', 'on_red'))
        NRMSE_per_exp = []
        fitness_per_exp = []
        for i in range(num_exp):
            print(colored('Exp #{}'.format(i+1), 'white', 'on_green'))
            NRMSE, fitness = short_term_memory(input_params=input_param, time_params=time_param, num_epoch=num_epoch)
            NRMSE_per_exp.append(NRMSE)            
            fitness_per_exp.append(fitness)
        NRMSE_means.append(np.mean(NRMSE_per_exp))
        NRMSE_stds.append(np.std(NRMSE_per_exp))
        fitness_means.append(np.mean(fitness_per_exp))
        fitness_stds.append(np.std(fitness_per_exp))
        base_influx_list.append(input_param['theta_in']['mean'])
        hold_time_list.append(time_param['t_hold'])

df_NRMSE = pd.DataFrame({'input range (species/sec)': base_influx_list,                                                                                  
                         'tau (sec)': hold_time_list,                                                                                  
                         'NRMSE_means': NRMSE_means,
                         'NRMSE_stds': NRMSE_stds})
df_NRMSE.to_csv('dat/ST_NRMSE_{}ep_{}ex.csv'.format(num_epoch, num_exp), index=False)

df_fitness = pd.DataFrame({'input range (species/sec)': base_influx_list,                                                                                  
                           'tau (sec)': hold_time_list,                                                                                  
                           'fitness_means': fitness_means,
                           'fitness_stds': fitness_stds})
df_fitness.to_csv('dat/ST_fitness_{}ep_{}ex.csv'.format(num_epoch, num_exp), index=False)
