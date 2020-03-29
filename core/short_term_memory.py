# Short-term memory task
ST_lookup = influx_lookup.copy() # lookup dict for short term memory task target for all reactions
for r_in, rate_in in ST_lookup.items():
    ST_target_per_reaction = []
    for ir_index in range(2, len(time_lookup1) + 1): # from index 2 to index end+1
        ST_target_per_reaction.append(rate_in[ir_index - 1] + 2*rate_in[ir_index - 2])
    ST_lookup.update({'{}'.format(r_in): ST_target_per_reaction})

plt.figure(figsize = (18,10))
plt.title('Plot of Short Term Memory Task target')
plt.xlabel('time')
plt.ylabel('ST memory target')
for reaction_index, (reaction, ST_target) in enumerate(ST_lookup.items()):
    plt.plot(time_lookup1[2:], ST_target[:-1], color=color_array[reaction_index], label=reaction)
handles, labels = plt.gca().get_legend_handles_labels()
by_label = OrderedDict(zip(labels, handles))
plt.legend(by_label.values(), by_label.keys(), loc='best')
plt.show()