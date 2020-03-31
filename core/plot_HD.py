import pandas as pd
import matplotlib.pyplot as plt


num_epoch = 10
num_exp = 10


df_NRMSE = pd.read_csv('dat/HD_NRMSE_{}ep_{}ex.csv'.format(num_epoch, num_exp))
df_NRMSE_mean = df_NRMSE.pivot(index='input range (species/sec)',columns='tau (sec)',values='NRMSE_means')
df_NRMSE_stds = df_NRMSE.pivot(index='input range (species/sec)',columns='tau (sec)',values='NRMSE_stds')
ax_NRMSE = df_NRMSE_mean.plot(kind='bar', yerr=df_NRMSE_stds.values, capsize=3, grid=True)
ax_NRMSE.grid(linestyle=':')
ax_NRMSE.legend(fancybox=True, framealpha=0.5)
plt.title('Hamming Distance Task')
plt.ylabel('NRMSE')
plt.tight_layout()
plt.savefig('plots/HD_NRMSE_{}ep_{}ex.png'.format(num_epoch, num_exp))
plt.savefig('plots/HD_NRMSE_{}ep_{}ex.eps'.format(num_epoch, num_exp))


df_fitness = pd.read_csv('dat/HD_fitness_{}ep_{}ex.csv'.format(num_epoch, num_exp))
df_fitness_mean = df_fitness.pivot(index='input range (species/sec)',columns='tau (sec)',values='fitness_means')
df_fitness_stds = df_fitness.pivot(index='input range (species/sec)',columns='tau (sec)',values='fitness_stds')
ax_fitness = df_fitness_mean.plot(kind='bar', yerr=df_fitness_stds.values, capsize=3, grid=True)
ax_fitness.grid(linestyle=':')
ax_fitness.legend(fancybox=True, framealpha=0.5)
plt.title('Hamming Distance Task')
plt.ylabel('fitness')
plt.tight_layout()
plt.savefig('plots/HD_fitness_{}ep_{}ex.png'.format(num_epoch, num_exp))
plt.savefig('plots/HD_fitness_{}ep_{}ex.eps'.format(num_epoch, num_exp))
