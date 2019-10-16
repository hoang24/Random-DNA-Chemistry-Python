import gillespy2
import numpy
import matplotlib.pyplot as plt
from random_dna_chem import RandomDNAStrandDisplacementCircuit
from params import input_params, time_params


class RandomDNAChemConstructionGillespy2(gillespy2.Model):
    '''
        Class to construct the random DNA strand displacement circuit using gillespy2
    '''

    def __init__(self, non_gillespy2_chem):
        super().__init__(self)
        self.randomDNAChem = non_gillespy2_chem
        self.convert_rates()
        self.convert_species()
        self.convert_reactions()
        self.set_timespan()

    def convert_rates(self):
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
            gillespy2_kin = gillespy2.Parameter(name=name_in, expression=rate_in[0])
            gillespy2_rates.append(gillespy2_kin)

        for index_out, (r_out, rate_out) in enumerate(self.randomDNAChem.rateConst_lookup['rate_OUT'].items()):
            name_out = 'efflux{}'.format(index_out)
            gillespy2_kout = gillespy2.Parameter(name=name_out, expression=rate_out[0])
            gillespy2_rates.append(gillespy2_kout)

        self.add_parameter(gillespy2_rates)

    def convert_species(self):
        gillespy2_species = []

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

        self.add_species(gillespy2_species)

    def convert_reactions(self):
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
        self.timespan(numpy.linspace(start=self.randomDNAChem.time_params['t_start'], 
                                     stop=self.randomDNAChem.time_params['t_end'], 
                                     num=1001))

if __name__ == '__main__':
    randomDNAChem = RandomDNAStrandDisplacementCircuit(input_params=input_params, 
                                                       time_params=time_params)

    gillespy2_model = RandomDNAChemConstructionGillespy2(non_gillespy2_chem=randomDNAChem)

    num_trajectories = 10
    gillespy2_results = gillespy2_model.run(number_of_trajectories=num_trajectories)

    for index in range(num_trajectories):
        trajectory = gillespy2_results[index]
        for species in randomDNAChem.species_lookup['S']:
            plt.plot(trajectory['time'], trajectory['{}'.format(species)])

    plt.xlabel('Time (s)')
    plt.ylabel('Number of molecules')
    plt.show()
