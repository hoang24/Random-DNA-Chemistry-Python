import gillespy2
import numpy


class RandomDNAChemPerturbationGillespy2(gillespy2.Model):
    '''
        Class to construct the random DNA strand displacement circuit using gillespy2
        Attributes:
            randomDNAChem (class): random DNA chemistry class
            rate_in_timeIndex (int): index of time i.e. (0: 0sec), (1: 0.01sec), (2: 0.21sec), (3:0.41sec), (4:0.61sec), (5:0.81sec)
            period_start (float): start of an simulation period i.e. 0
            period_end (float): end of an simulation period i.e. 0.01
            previous_gillespy2_result (obj): gillespy2 result after run from the previous period
    '''

    def __init__(self, non_gillespy2_chem, rate_in_timeIndex, period_start, period_end, numel, previous_gillespy2_result):
        '''
            Init method to call the parent class (gillespy2.Model) and the class methods
            Args:
                non_gillespy2_chem (class): (random DNA chem centric) class that is not using gillespy2
                rate_in_timeIndex (int): index of time i.e. (0: 0sec), (1: 0.01sec), (2: 0.21sec), (3:0.41sec), (4:0.61sec), (5:0.81sec)
                period_start (float): start of an simulation period i.e. 0
                period_end (float): end of an simulation period i.e. 0.01
                previous_gillespy2_result (obj): gillespy2 result after run from the previous period
        '''

        super().__init__(self)
        self.randomDNAChem = non_gillespy2_chem
        self.rate_in_timeIndex = rate_in_timeIndex
        self.period_start = period_start
        self.period_end = period_end
        self.numel = numel
        self.previous_gillespy2_result = previous_gillespy2_result
        self.convert_rates()
        self.convert_species()
        self.convert_reactions()
        self.set_timespan()

    def convert_rates(self):
        '''
            Method to convert rate constants from non-gillespy2 Random DNA Chemistry to the gillespy2 Random DNA Chemistry
        '''

        gillespy2_rates = []

        for index_bind, (r_bind, rate_bind) in enumerate(self.randomDNAChem.rateConst_lookup['rate_BIND'].items()):
            name_bind = 'bind{}'.format(index_bind)
            gillespy2_kbind = gillespy2.Parameter(name=name_bind, expression=rate_bind[0])
            gillespy2_rates.append(gillespy2_kbind)

        for index_displace, (r_displace, rate_displace) in enumerate(self.randomDNAChem.rateConst_lookup['rate_DISPLACE'].items()):
            name_displace = 'displace{}'.format(index_displace)
            gillespy2_kdisplace = gillespy2.Parameter(name=name_displace, expression=rate_displace[0])
            gillespy2_rates.append(gillespy2_kdisplace)

        for index_in, (r_in, rate_in) in enumerate(self.randomDNAChem.rateConst_lookup['rate_IN'].items()):
            name_in = 'influx{}'.format(index_in)
            gillespy2_kin = gillespy2.Parameter(name=name_in, expression=rate_in[self.rate_in_timeIndex])
            gillespy2_rates.append(gillespy2_kin)

        for index_out, (r_out, rate_out) in enumerate(self.randomDNAChem.rateConst_lookup['rate_OUT'].items()):
            name_out = 'efflux{}'.format(index_out)
            gillespy2_kout = gillespy2.Parameter(name=name_out, expression=rate_out[0])
            gillespy2_rates.append(gillespy2_kout)

        self.add_parameter(gillespy2_rates)

    def convert_species(self):
        '''
            Method to convert species from non-gillespy2 Random DNA Chemistry to the gillespy2 Random DNA Chemistry
        '''
        gillespy2_species = []

        if self.previous_gillespy2_result == None:
            for u, ucon in self.randomDNAChem.concentration_lookup['conU'].items():
                gillespy2_u = gillespy2.Species(name=u, initial_value=ucon[0])
                gillespy2_species.append(gillespy2_u)

            for l, lcon in self.randomDNAChem.concentration_lookup['conL'].items():
                gillespy2_l = gillespy2.Species(name=l, initial_value=lcon[0])
                gillespy2_species.append(gillespy2_l)

            for f, fcon in self.randomDNAChem.concentration_lookup['conF'].items():
                gillespy2_f = gillespy2.Species(name=f, initial_value=fcon[0])
                gillespy2_species.append(gillespy2_f)

            for p, pcon in self.randomDNAChem.concentration_lookup['conP'].items():
                gillespy2_p = gillespy2.Species(name=p, initial_value=pcon[0])
                gillespy2_species.append(gillespy2_p)
        else:
            for u, ucon in self.randomDNAChem.concentration_lookup['conU'].items():
                gillespy2_u = gillespy2.Species(name=u, initial_value=int(self.previous_gillespy2_result[u][-1]))
                gillespy2_species.append(gillespy2_u)

            for l, lcon in self.randomDNAChem.concentration_lookup['conL'].items():
                gillespy2_l = gillespy2.Species(name=l, initial_value=int(self.previous_gillespy2_result[l][-1]))
                gillespy2_species.append(gillespy2_l)

            for f, fcon in self.randomDNAChem.concentration_lookup['conF'].items():
                gillespy2_f = gillespy2.Species(name=f, initial_value=int(self.previous_gillespy2_result[f][-1]))
                gillespy2_species.append(gillespy2_f)

            for p, pcon in self.randomDNAChem.concentration_lookup['conP'].items():
                gillespy2_p = gillespy2.Species(name=p, initial_value=int(self.previous_gillespy2_result[p][-1]))
                gillespy2_species.append(gillespy2_p)

        self.add_species(gillespy2_species)

    def convert_reactions(self):
        '''
            Method to convert reactions from non-gillespy2 Random DNA Chemistry to the gillespy2 Random DNA Chemistry
        '''
        gillespy2_reactions = []

        for i_bind, r_bind in enumerate(self.randomDNAChem.reaction_lookup['R_BIND']):
            for species1 in self.listOfSpecies:
                if species1 == r_bind[0:2]:
                    bind_reactant1 = species1
            for species2 in self.listOfSpecies:
                if species2 == r_bind[5:7]:
                    bind_reactant2 = species2
            for species3 in self.listOfSpecies:
                if species3 == r_bind[12:17]:
                    bind_product = species3

            gillespy2_bind = gillespy2.Reaction(name=r_bind, 
                                                rate=list(self.listOfParameters.values())[i_bind],
                                                reactants={bind_reactant1:1, bind_reactant2:1}, 
                                                products={bind_product:1})
            gillespy2_reactions.append(gillespy2_bind)

        for i_displace, r_displace in enumerate(self.randomDNAChem.reaction_lookup['R_DISPLACE']):
            for species1 in self.listOfSpecies:
                if species1 == r_displace[0:2]:
                    displace_reactant1 = species1
            for species2 in self.listOfSpecies:
                if species2 == r_displace[5:9]:
                    displace_reactant2 = species2
            for species3 in self.listOfSpecies:
                if species3 == r_displace[14:18]:
                    displace_product1 = species3
            for species3 in self.listOfSpecies:
                if species3 == r_displace[21:23]:
                    displace_product2 = species3

            gillespy2_displace = gillespy2.Reaction(name=r_displace, 
                                                    rate=list(self.listOfParameters.values())[i_displace],
                                                    reactants={displace_reactant1:1, displace_reactant2:1}, 
                                                    products={displace_product1:1, displace_product2:1})
            gillespy2_reactions.append(gillespy2_displace)

        for i_in, r_in in enumerate(self.randomDNAChem.reaction_lookup['R_IN']):
            if len(r_in) == 10:
                for species1 in self.listOfSpecies:
                    if species1 == r_in[6:10]:
                        in_product = species1
            elif len(r_in) == 8:
                for species1 in self.listOfSpecies:
                    if species1 == r_in[6:8]:
                        in_product = species1
            else:
                raise BaseException

            gillespy2_in = gillespy2.Reaction(name=r_in, 
                                              rate=list(self.listOfParameters.values())[i_in],
                                              reactants={}, 
                                              products={in_product:1})
            gillespy2_reactions.append(gillespy2_in)

        for i_out, r_out in enumerate(self.randomDNAChem.reaction_lookup['R_OUT']):
            if len(r_out) == 10:
                for species1 in self.listOfSpecies:
                    if species1 == r_out[0:4]:
                        out_reactant = species1
            elif len(r_out) == 8:
                for species1 in self.listOfSpecies:
                    if species1 == r_out[0:2]:
                        out_reactant = species1
            else:
                raise BaseException

            gillespy2_out = gillespy2.Reaction(name=r_out, 
                                               rate=list(self.listOfParameters.values())[i_out],
                                               reactants={out_reactant:1}, 
                                               products={})
            gillespy2_reactions.append(gillespy2_out)

        self.add_reaction(gillespy2_reactions)

    def set_timespan(self):
        '''
            Method to set the simulation timing of the gillespy2 Random DNA Chemistry
        '''
        self.timespan(numpy.linspace(start=self.period_start, 
                                     stop=self.period_end, 
                                     num=self.numel))

# if __name__ == '__main__':

#    import matplotlib.pyplot as plt
#    from random_dna_chem import RandomDNAStrandDisplacementCircuit
#    from params import input_params, time_params
#    from collections import OrderedDict

#     randomDNAChem = RandomDNAStrandDisplacementCircuit(input_params=input_params, 
#                                                        time_params=time_params)

#     color_array = ['#000000', '#0000FF', '#00FF00', '#00FFFF', '#000080',
#                    '#008000', '#008080', '#800000', '#800080', '#808000',
#                    '#808080', '#C0C0C0', '#FF0000', '#FF00FF', '#FFFF00',
#                    '#8B0000', '#006400', '#BDB76B', '#008B8B', '#191970']
    
#     plt.figure(figsize = (18,10))
#     plt.title('Plot of all molecules over a simulation time')
#     plt.xlabel('Time (s)')
#     plt.ylabel('Number of molecules')

#     gillespy2_results = []
#     num_trajectories = 1

#     # Non-perturb period
#     gillespy2_model = RandomDNAChemPerturbationGillespy2(non_gillespy2_chem=randomDNAChem,
#                                                          rate_in_timeIndex=0, # index 0 for time t=0
#                                                          period_start=randomDNAChem.time_params['time_array'][0],
#                                                          period_end=randomDNAChem.time_params['time_array'][1],
#                                                          numel=1001,
#                                                          previous_gillespy2_result=None)
#     gillespy2_result = gillespy2_model.run(number_of_trajectories=num_trajectories)
#     gillespy2_results.append(gillespy2_result)
#     for index in range(num_trajectories):
#         trajectory = gillespy2_result[index]
#         for species_index, species in enumerate(randomDNAChem.species_lookup['S']):
#             species_plot = plt.plot(trajectory['time'],
#                                     trajectory['{}'.format(species)],
#                                     color=color_array[species_index],
#                                     label=species)

#     # Perturb period
#     for time_index in range(1, len(randomDNAChem.time_params['time_array']) - 1):
#         time_offset = randomDNAChem.time_params['t_perturb'] + randomDNAChem.time_params['t_hold'] * (time_index - 1)
#         gillespy2_model = RandomDNAChemPerturbationGillespy2(non_gillespy2_chem=randomDNAChem,
#                                                              rate_in_timeIndex=time_index, # index 1 for time t=1 and so on
#                                                              period_start=randomDNAChem.time_params['time_array'][time_index] - time_offset,
#                                                              period_end=randomDNAChem.time_params['time_array'][time_index + 1] - time_offset,
#                                                              numel=1001,
#                                                              previous_gillespy2_result=gillespy2_result)
#         gillespy2_result = gillespy2_model.run(number_of_trajectories=num_trajectories)
#         gillespy2_results.append(gillespy2_result)
#         for index in range(num_trajectories):
#             trajectory = gillespy2_result[index]
#             for species_index, species in enumerate(randomDNAChem.species_lookup['S']):
#                 species_plot = plt.plot(trajectory['time'] + time_offset,
#                                         trajectory['{}'.format(species)],
#                                         color=color_array[species_index],
#                                         label=species)

#     handles, labels = plt.gca().get_legend_handles_labels()
#     by_label = OrderedDict(zip(labels, handles))
#     plt.legend(by_label.values(), by_label.keys(), loc='best')
#     # plot_name = 'random_dna_chem_with_perturb'
#     try:
#         plot_name
#     except NameError:
#         plt.show()
#     else:
#         plt.savefig('plots/' + plot_name + '.eps')
