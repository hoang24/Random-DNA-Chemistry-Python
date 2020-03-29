# Long-term memory task
t_hold = randomDNAChem.time_params['t_hold']
t_hold_index = time_lookup1.index(t_hold)
t_hold_32_index = time_lookup1.index(t_hold*(3/2))
LT_lookup = influx_lookup.copy() # lookup dict for short term memory task target for all reactions
for r_in, rate_in in LT_lookup.items():
    LT_target_per_reaction = []
    for ir_index in range(t_hold_32_index, len(time_lookup1) + t_hold_index): # from index hold_time to index end + index (3/4)*hold_time
        LT_target_per_reaction.append(rate_in[ir_index - t_hold_index] + (1/2)*rate_in[ir_index - t_hold_32_index])
    LT_lookup.update({'{}'.format(r_in): LT_target_per_reaction})

plt.figure(figsize = (18,10))
plt.title('Plot of Long Term Memory Task target')
plt.xlabel('time')
plt.ylabel('LT memory target')
for reaction_index, (reaction, LT_target) in enumerate(LT_lookup.items()):
    plt.plot(time_lookup1[t_hold_32_index:], LT_target[:-t_hold_index], color=color_array[reaction_index], label=reaction)
handles, labels = plt.gca().get_legend_handles_labels()
by_label = OrderedDict(zip(labels, handles))
plt.legend(by_label.values(), by_label.keys(), loc='best')
plt.show()
