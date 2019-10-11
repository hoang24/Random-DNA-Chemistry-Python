'''
    Script to simulate (create and run) the random DNA strand displacement circuit
'''

from Random_DNA_Chem import RandomDNAStrandDisplacementCircuit
from params import input_params, time_params


random_DNA_network = RandomDNAStrandDisplacementCircuit(input_params=input_params, time_params=time_params)

print('Number of species: {}'.format(random_DNA_network.species_lookup['nS']))
print('Species set: {}'.format(random_DNA_network.species_lookup['S']))

print('Order lookup dictionary: ')
print(random_DNA_network.order_lookup)
print(random_DNA_network.species_lookup['DS'])
print(random_DNA_network.species_lookup['SS']) 

for key, value in random_DNA_network.reaction_lookup.items():
    print('{} ~~ {}'.format(key, value))
    print('\n')

for key, value in random_DNA_network.rateConst_lookup.items():
    print('{} ~~ {}'.format(key, value))
    print('\n')

print('intial concentration: {}'.format(random_DNA_network.concentration_lookup))

print(random_DNA_network.input_params)
print(random_DNA_network.time_params)