from params import input_params5, time_params1
from hamming_dist import hamming_dist
from short_term_memory import short_term_memory
from long_term_memory import long_term_memory
import numpy as np
import matplotlib.pyplot as plt

NRMSE, fitness = hamming_dist(input_params=input_params5, time_params=time_params1, num_epoch=10, plot_chem=True, error_plot=True)
NRMSE, fitness = short_term_memory(input_params=input_params5, time_params=time_params1, num_epoch=50, plot_chem=True, error_plot=True)
NRMSE, fitness = long_term_memory(input_params=input_params5, time_params=time_params1, num_epoch=50, plot_chem=True, error_plot=True)

# print('NRMSE:', NRMSE)
# print('fitness:', fitness)