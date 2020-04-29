from params import input_params5, time_params1
from hamming_dist import hamming_dist
import numpy as np
import matplotlib.pyplot as plt

NRMSE, fitness = hamming_dist(input_params=input_params5, time_params=time_params1, num_epoch=1, plot_chem=True)

print('NRMSE:', NRMSE)
print('fitness:', fitness)