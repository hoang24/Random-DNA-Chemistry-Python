import pandas as pd
import matplotlib.pyplot as plt


num_epoch = 10
num_exp = 20

plt.rc('font', family='serif')
plt.rc('xtick', labelsize='x-large')
plt.rc('ytick', labelsize='x-large')
fig = plt.figure(figsize=(4, 3))
ax = fig.add_subplot(1, 1, 1)

df_NRMSE = pd.read_csv('dat/LT_NRMSE_{}ep_{}ex.csv'.format(num_epoch, num_exp))
df_NRMSE_mean = df_NRMSE.pivot(index='input range (species/sec)',columns='tau (sec)',values='NRMSE_means')
df_NRMSE_stds = df_NRMSE.pivot(index='input range (species/sec)',columns='tau (sec)',values='NRMSE_stds')
ax_NRMSE = df_NRMSE_mean.plot(kind='bar', yerr=df_NRMSE_stds.values, grid=True, rot=0, error_kw=dict(lw=1, capsize=2, capthick=1))
ax_NRMSE.grid(linestyle=':')
ax_NRMSE.legend(loc='upper left', bbox_to_anchor=(0.12, 1), handlelength=1.0, framealpha=0.5, title=r'$\tau$ (s)', fontsize='medium')
# plt.title('Long-term Memory Task')
ax_NRMSE.set_xlabel(r'Base Influx Rate $\theta_{in}$ (species/sec)', fontfamily='serif', fontsize='x-large')
ax_NRMSE.set_ylabel('NRMSE', fontsize='x-large')
plt.tight_layout()
plt.savefig('plots/LT_NRMSE_{}ep_{}ex.png'.format(num_epoch, num_exp))
plt.savefig('plots/LT_NRMSE_{}ep_{}ex.eps'.format(num_epoch, num_exp))


# df_fitness = pd.read_csv('dat/LT_fitness_{}ep_{}ex.csv'.format(num_epoch, num_exp))
# df_fitness_mean = df_fitness.pivot(index='input range (species/sec)',columns='tau (sec)',values='fitness_means')
# df_fitness_stds = df_fitness.pivot(index='input range (species/sec)',columns='tau (sec)',values='fitness_stds')
# ax_fitness = df_fitness_mean.plot(kind='bar', yerr=df_fitness_stds.values, capsize=3, grid=True)
# ax_fitness.grid(linestyle=':')
# ax_fitness.legend(fancybox=True, framealpha=0.5)
# plt.title('Long-term Memory Task')
# plt.ylabel('fitness')
# plt.tight_layout()
# plt.savefig('plots/LT_fitness_{}ep_{}ex.png'.format(num_epoch, num_exp))
# plt.savefig('plots/LT_fitness_{}ep_{}ex.eps'.format(num_epoch, num_exp))
