import pandas as pd
import matplotlib.pyplot as plt
from collections import OrderedDict

# model = ['Goudarzi, 2013', 'Yahiro, 2018', 'Our model']
# ST_NRMSE = [0.24, 0.17, 0.055]
# LT_NRMSE = [0.12, 0.18, 0.05]

# CRN_model = ['Goudarzi, 2013', 'Goudarzi, 2013', 'Yahiro, 2018', 'Yahiro, 2018', 'Our model', 'Our model']
# Task = ['Short-term Memory Task',
#         'Long-term Memory Task',
#         'Short-term Memory Task',
#         'Long-term Memory Task',
#         'Short-term Memory Task',
#         'Long-term Memory Task']
# NRMSE_mean = [0.24, 0.12, 0.17, 0.18, 0.055, 0.05]
# NRMSE_stds = [0.05, 0.025, 0.025, 0.026, 0.025, 0.0125]

# df_NRMSE = pd.DataFrame({
#         'CRN Model': CRN_model,
#         'Task': Task,
#         'NRMSE_mean': NRMSE_mean,
#         'NRMSE_stds': NRMSE_stds
#     })


plt.rc('font', family='serif')
plt.rc('xtick', labelsize='large')
plt.rc('ytick', labelsize='large')
fig = plt.figure(figsize=(4, 3))
ax = fig.add_subplot(1, 1, 1)

df_NRMSE = pd.read_csv('dat/model_comparison.csv')
df_NRMSE_mean = df_NRMSE.pivot(index='CRN Model',columns='Task',values='NRMSE_mean')
df_NRMSE_stds = df_NRMSE.pivot(index='CRN Model',columns='Task',values='NRMSE_stds')

ax_NRMSE = df_NRMSE_mean.plot(kind='bar', yerr=df_NRMSE_stds, capsize=3, grid=True, rot=0)
ax_NRMSE.grid(linestyle=':')
ax_NRMSE.legend(loc='best', handlelength=2.0, framealpha=0.5, fontsize='x-large')
# plt.title('Model Comparison')
ax_NRMSE.set_xlabel('CRN Model', fontfamily='serif', fontsize='xx-large')
ax_NRMSE.set_ylabel('NRMSE', fontsize='xx-large')
plt.tight_layout()
plt.savefig('plots/model_comparison_2020.png')
plt.savefig('plots/model_comparison_2020.eps')
