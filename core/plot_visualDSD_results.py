import pandas as pd
import matplotlib.pyplot as plt
from collections import OrderedDict


plt.rc('font', family='serif')
plt.rc('xtick', labelsize='large')
plt.rc('ytick', labelsize='large')
fig = plt.figure(figsize=(4, 3))
ax = fig.add_subplot(1, 1, 1)

df_NRMSE = pd.read_csv('visualDSD/dat_dsd_tau_sweep.csv')

df_NRMSE_mean = df_NRMSE.pivot(index='tau (sec)',columns='Time-series Task',values='NRMSE_means')
df_NRMSE_stds = df_NRMSE.pivot(index='tau (sec)',columns='Time-series Task',values='NRMSE_stds')

ax_NRMSE = df_NRMSE_mean.plot(kind='bar', yerr=df_NRMSE_stds, capsize=3, grid=True, rot=0)
ax_NRMSE.grid(linestyle=':')
ax_NRMSE.legend(loc='best', handlelength=2.0, framealpha=0.5, fontsize='large')

ax_NRMSE.set_xlabel(r'Input Hold Time $\tau$ (seconds)', fontfamily='serif', fontsize='x-large')
ax_NRMSE.set_ylabel('NRMSE', fontsize='x-large')

plt.tight_layout()
plt.savefig('visualDSD/plot_dsd_tau_sweep.png')
plt.savefig('visualDSD/plot_dsd_tau_sweep.eps')