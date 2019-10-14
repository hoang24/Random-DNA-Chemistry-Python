'''
    Script to simulate (create and run) the random DNA strand displacement circuit
'''

from random_dna_chem import RandomDNAStrandDisplacementCircuit
from params import input_params, time_params


class RandomDNAChemDisplay(RandomDNAStrandDisplacementCircuit):

    def __init__(self):
        super().__init__(input_params=input_params, time_params=time_params)
        print('\n')

    def display_species(self):
        print('Number of species: {}'.format(self.species_lookup['nS']))
        print('Species set: {}'.format(self.species_lookup['S']))
        print('\n')
        
    def display_reactions(self):
        print('Number of reactions: {}'.format(self.reaction_lookup['nR']))
        print('Reaction set: {}'.format(self.reaction_lookup['R']))
        print('\n')

    def display_species_concentration(self):
        print('Upper single strands and concentrations: ')
        for u, conu in self.concentration_lookup['conU'].items():
            print('{}    {}'.format(u, conu))
        print('\n')

        print('Lower single strands and concentrations: ')
        for l, conl in self.concentration_lookup['conL'].items():
            print('{}    {}'.format(l, conl))
        print('\n')

        print('Full double strands and concentrations: ')
        for f, conf in self.concentration_lookup['conF'].items():
            print('{}    {}'.format(f, conf))
        print('\n')

        print('Partial double strands and concentrations: ')
        for p, conp in self.concentration_lookup['conP'].items():
            print('{}    {}'.format(p, conp))
        print('\n')

    def display_order(self):
        print('Order lookup dictionary (higher number higher order): ')
        print(self.order_lookup)
        print('\n')

    def display_reactions_rates(self):
        print('Displacement reactions and rates:')
        for r_displace, rate_displace in self.rateConst_lookup['rate_DISPLACE'].items():
            print('{}    {}'.format(r_displace, rate_displace))
        print('\n')

        print('Binding reactions and rates:')
        for r_bind, rate_bind in self.rateConst_lookup['rate_BIND'].items():
            print('{}    {}'.format(r_bind, rate_bind))
        print('\n')

        print('Influx reactions and rates:')
        for r_in, rate_in in self.rateConst_lookup['rate_IN'].items():
            print('{}    {}'.format(r_in, rate_in))
        print('\n')

        print('Efflux reactions and rates:')
        for r_out, rate_out in self.rateConst_lookup['rate_OUT'].items():
            print('{}    {}'.format(r_out, rate_out))
        print('\n')

    def display_inputs(self):
        print('Input parameters: ')
        for iparam, iparam_value in self.input_params.items():
            if isinstance(iparam_value, dict):
                print('{}: '.format(iparam))
                for key, value in iparam_value.items():
                    print('  * {} = {}'.format(key, value))
            else:
                print('{} = {}'.format(iparam, iparam_value))
        print('\n')

    def display_times(self):
        print('Time parameters: ')
        for tparam, tparam_value in self.time_params.items():
            print('{} = {}'.format(tparam, tparam_value))
        print('\n')

if __name__ == '__main__':
    randomDNAChem = RandomDNAChemDisplay()
    randomDNAChem.display_species()
    randomDNAChem.display_reactions()
    randomDNAChem.display_species_concentration()
    randomDNAChem.display_order()
    randomDNAChem.display_reactions_rates()
    randomDNAChem.display_inputs()
    randomDNAChem.display_times()