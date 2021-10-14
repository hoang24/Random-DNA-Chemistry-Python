import numpy as np
import time
from string import Template


class RandomDSDChemistry():
    '''
        Random DNA Strand Displacement (DSD) Chemistry 2nd implementation
        Generate Visual DSD Script
        Attributes:
            input_params (dict): lookup dictionary for input parameters
            time_params (dict): lookup dictionary for timing parameters
            species_lookup (dict): lookup dictionary for species names and number of species
            reaction_lookup (dict): lookup dictionary for orders of partial double strands
            concentration_lookup (dict): lookup dictionary for all concentrations of all species
            rateConst_lookup (dict): lookup dictionary for reaction rate constants of all reactions
    '''
    def __init__(self, input_params, time_params):
        '''
            Python initialization method.
            Args:
                input_params (dict): input parameters
                time_params (dict): timing parameters
        '''
        # User-defined parameters
        self.input_params = input_params
        self.time_params = time_params

        # Method to construct the random DSD chemistry
        self.create_species_single()
        self.create_domain_lookup()
        self.create_species_double()
        self.create_species_sets()
        self.choose_influx()
        self.choose_efflux()
        self.create_species_lookup()
        self.create_initial_condition()
        self.create_species_count_lookup()
        self.create_initial_rateConst_binding()
        self.create_reaction_binding()
        self.create_reaction_displacement()
        self.create_initial_rateConst_displacement()
        self.create_initial_rateConst_influx()
        self.create_reaction_influx()
        self.create_initial_rateConst_efflux()
        self.create_reaction_efflux()
        self.create_reaction_sets()
        self.create_reaction_lookup()
        self.create_rateConst_lookup()
        self.update_time_params()
        self.update_influx_rate()
        self.create_perturbation_lookup()

    def create_species_single(self):
        '''
            Method to initialize the set of single DNA strands
            nU (int): number of upper single strands
            nL (int): number of lower single strands
            U (list): set of upper single strands
            L (list): set of lower single strands
        '''
        self.nL = int(round(self.input_params['n'] / (1 + self.input_params['p'])))
        self.L = []
        for nl in range(self.nL):
            self.L.append('L{}'.format(nl))

        self.nU = int(self.input_params['n'] - self.nL)
        self.U = []
        for nu in range(self.nU):
            self.U.append('U{}'.format(nu))

    def create_domain_lookup(self):
        '''
            Method to assign the switching domain types
            switching_domains (list of str): list of switching domain types (e.g. A, B, C, etc.)
            identity_domains (list of str): list of identity domain types (e.g. 1, 2, 3, 4, etc.)
        '''
        SWITCHING_DOMAINS = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'
        self.switching_domains = list(SWITCHING_DOMAINS[:self.input_params['s']])

        if self.nU >= self.nL:
            self.identity_domains = [i[1:len(i)] for i in self.U]
        else:
            self.identity_domains = [i[1:len(i)] for i in self.L]

        self.domain_lookup = {'U': {}, 'L': {}}
        for u in self.U:
            switching_domain = str(np.random.choice(self.switching_domains))
            self.domain_lookup['U'].update({u: ['M', 'M', switching_domain, u[-1]]})

        for l in self.L:
            switching_domain = str(np.random.choice(self.switching_domains))
            self.domain_lookup['L'].update({l: ['M', 'M', switching_domain, l[-1]]})

    def create_species_double(self):
        '''
            Method to initialize the set of double DNA strands
            nF (int): number of full double strands
            nPS (int): number of strong (2/3) partial double strands
            nPW (int): number of weak (1/3) partial double strands
            F (list): set of full double strands
            PS (list): set of strong (2/3) partial double strands
            PW (list): set of weak (1/3) partial double strands
        '''
        self.F = []
        self.PS = []
        self.PW = []
        for u in self.U:
            for l in self.L:
                if self.domain_lookup['U'][u][2] == self.domain_lookup['L'][l][2]: # check if switching domains are the same type for upper and lower strands
                    if self.domain_lookup['U'][u][3] == self.domain_lookup['L'][l][3]: # check if identity domains are the same type for upper and lower strands
                        self.F.append(u + l)
                    else:
                        self.PS.append(u + l)
                else:
                    self.PW.append(u + l)
        self.nF = len(self.F)
        self.nPS = len(self.PS)
        self.nPW = len(self.PW)

    def create_species_sets(self):
        '''
            Method to create species sets
            S (list): set of all species in the chemistry
            SS (list): set of all single strands (including upper and lower strands)
            P (list): set of all partial double strands (including the strong 2/3 double strands and the weak 1/3 double strands)
            DS (list): set of all double strands (inluding full double and strong+weak partial double strands)
            nS (int): total number of species in the chemistry
            nSS (int): number of single stranded species in the chemistry
            nP (int): number of partial double species in the chemistry (strong and weak types)
            nDS (int): number of double stranded species in the chemistry
        '''
        self.S = self.U + self.L + self.F + self.PS + self.PW
        self.SS = self.U + self.L
        self.P = self.PS + self.PW
        self.DS = self.F + self.P

        self.nS = len(self.S)
        self.nSS = len(self.SS)
        self.nP = len(self.P)
        self.nDS = len(self.DS)

    def choose_influx(self):
        '''
            Method to randomly choose a single stranded species to be influx species
        '''
        self.nI = self.input_params['a_in']
        self.I = list(np.random.choice(a=self.SS, size=self.nI, replace=False))

    def choose_efflux(self):
        '''
            Method to create a set of efflux species that include all species in the chemistry
        '''
        self.nO = self.nS
        self.O = self.S

    def create_species_lookup(self):
        '''
            Methods to create a lookup dictionary for sets of species and number of species
            Returns:
                species_lookup (dict): lookup dictionary for sets of species and number of species
        '''
        self.species_lookup = {
            'U': self.U,
            'L': self.L,
            'F': self.F,
            'PS': self.PS,
            'PW': self.PW,
            'SS': self.SS,
            'P': self.P,
            'DS': self.DS,
            'S': self.S,
            'I': self.I,
            'O': self.O,
            'nU': self.nU,
            'nL': self.nL,
            'nF': self.nF,
            'nPS': self.nPS,
            'nPW': self.nPW,
            'nP': self.nP,
            'nSS': self.nSS,
            'nDS': self.nDS,
            'nS': self.nS,
            'nI': self.nI,
            'nO': self.nO
        }

    def create_initial_condition(self):
        '''
            Method to create sets of initial condition for each DNA strand types.
            ic_U (liat): initial species count for all upper single strands
            ic_L (liat): initial species count for all lower single strands
            ic_F (liat): initial species count for all full double strands
            ic_PS (liat): initial species count for all strong (2/3) partial double strands
            ic_PW (liat): initial species count for all weak (1/3) partial double strands
        '''
        ic_U = np.random.randint(self.input_params['d']['lb'], self.input_params['d']['ub'], self.species_lookup['nU'])
        ic_L = np.random.randint(self.input_params['d']['lb'], self.input_params['d']['ub'], self.species_lookup['nL'])
        ic_F = np.random.randint(self.input_params['d']['lb'], self.input_params['d']['ub'], self.species_lookup['nF'])
        ic_PS = np.random.randint(self.input_params['d']['lb'], self.input_params['d']['ub'], self.species_lookup['nPS'])
        ic_PW = np.random.randint(self.input_params['d']['lb'], self.input_params['d']['ub'], self.species_lookup['nPW'])

        self.ic_U = [int(icu) for icu in ic_U]
        self.ic_L = [int(icl) for icl in ic_L]
        self.ic_F = [int(icf) for icf in ic_F]
        self.ic_PS = [int(icps) for icps in ic_PS]
        self.ic_PW = [int(icpw) for icpw in ic_PW]

    def create_species_count_lookup(self):
        '''
            Method to create a lookup dictionary for all species count for each species
                initial_species_count_lookup (dict): lookup dictionary for initial species count of each species
                species_count_lookup (dict): lookup dictionary for all species count of each species
        '''
        conU = {}
        for u_idx, u in enumerate(self.U):
            conU.update({u: [self.ic_U[u_idx]]})

        conL = {}
        for l_idx, l in enumerate(self.L):
            conL.update({l: [self.ic_L[l_idx]]})

        conF = {}
        for f_idx, f in enumerate(self.F):
            conF.update({f: [self.ic_F[f_idx]]})

        conPS = {}
        for ps_idx, ps in enumerate(self.PS):
            conPS.update({ps: [self.ic_PS[ps_idx]]})

        conPW = {}
        for pw_idx, pw in enumerate(self.PW):
            conPW.update({pw: [self.ic_PW[pw_idx]]})

        self.species_count_lookup = {
            'conU': conU,
            'conL': conL,
            'conF': conF,
            'conPS': conPS,
            'conPW': conPW
        }

    def decode_double(self, ds):
        '''
            Method to decode the double strands to the upper and lower strands 
            Args:
                ds (str): double strand
            Returns:
                u (str): upper strand that made up the double strand
                l (str): lower strand that made up the double strand
        '''
        for i, d in enumerate(ds):
            if d == 'U':
                u_index = i
            elif d == 'L':
                l_index = i
        u = ds[u_index:l_index]
        l = ds[l_index:len(ds)]
        return u, l

    def create_initial_rateConst_binding(self):
        '''
            Method to create sets reaction rate constants for all binding reactions.
            k_BIND (list): set of reaction rate constants for all binding reactions
        '''

        norm_dist_bind = np.random.normal(loc=self.input_params['theta']['mean'], 
                                          scale=np.sqrt(self.input_params['theta']['variance']),
                                          size=self.nDS)
        norm_dist_bind = np.abs(norm_dist_bind)
        self.input_params['theta'].update({'norm_dist_bind': norm_dist_bind})
        self.k_BIND = []
        for norm_dist_bind_value in norm_dist_bind:
            self.k_BIND.append(np.abs(norm_dist_bind_value))

    def create_reaction_binding(self):
        '''
            Method to create a set of binding reactions
            R_BIND (list): set of binding reactions
            nR_BIND (int): number of binding reactions
        '''
        self.R_BIND = []
        for i_ds, ds in enumerate(self.DS):
            u, l = self.decode_double(ds)
            R_bind = f"{u}() + {l}() -> {ds}()"
            self.R_BIND.append(R_bind)
        self.nR_BIND = len(self.R_BIND)

    def create_reaction_displacement(self):
        '''
            Method to create a set of displacement reactions
            R_DISPLACE (list): set of displacement reactions
            nR_DISPLACE (int): number of displacement reactions
        '''

        self.R_DISPLACE = []
        for u in self.U: # for each upper strand
            for ps in self.PS: # for each strong partial double strand
                du, dl = self.decode_double(ps)
                if u not in ps: # if single strand is not part of double strand
                    if (self.domain_lookup['U'][u][2] == self.domain_lookup['L'][dl][2]) and (self.domain_lookup['U'][u][3] == self.domain_lookup['L'][dl][3]): # if the switching and identity domain of the upper strand matches that of the lower part of the double strand
                        prod = u + dl
                        R_displace = f"{u}() + {ps}() -> {du}() + {prod}()"
                        self.R_DISPLACE.append(R_displace)

            for pw in self.PW: # for each weak partial double strand
                du, dl = self.decode_double(pw)
                if u not in pw: # if single strand is not part of double strand
                    if self.domain_lookup['U'][u][2] == self.domain_lookup['L'][dl][2]: # if the switching domain of the upper strand matches that of the lower part of the double strand
                        prod = u + dl
                        R_displace = f"{u}() + {pw}() -> {du}() + {prod}()"
                        self.R_DISPLACE.append(R_displace)

        for l in self.L: # for each lower strand
            for ps in self.PS: # for each strong partial double strand
                du, dl = self.decode_double(ps)
                if l not in ps: # if single strand is not part of double strand
                    if (self.domain_lookup['U'][du][2] == self.domain_lookup['L'][l][2]) and (self.domain_lookup['U'][du][3] == self.domain_lookup['L'][l][3]): # if the switching and identity domain of the upper strand matches that of the lower part of the double strand
                        prod = du + l
                        R_displace = f"{l}() + {ps}() -> {dl}() + {prod}()"
                        self.R_DISPLACE.append(R_displace)

            for pw in self.PW: # for each weak partial double strand
                du, dl = self.decode_double(pw)
                if l not in pw: # if single strand is not part of double strand
                    if self.domain_lookup['U'][du][2] == self.domain_lookup['L'][l][2]: # if the switching domain of the lower strand matches that of the upper part of the double strand
                        prod = du + l
                        R_displace = f"{l}() + {pw}() -> {dl}() + {prod}()"
                        self.R_DISPLACE.append(R_displace)

        self.nR_DISPLACE = len(self.R_DISPLACE)

    def create_initial_rateConst_displacement(self):
        '''
            Method to create sets reaction rate constants for all displacement reactions.
            k_DISPLACE (list): list of reaction rate constants for all displacement reactions
        '''

        norm_dist_displace = np.random.normal(loc=self.input_params['theta']['mean'], 
                                              scale=np.sqrt(self.input_params['theta']['variance']),
                                              size=self.nR_DISPLACE)
        norm_dist_displace = np.abs(norm_dist_displace)
        self.input_params['theta'].update({'norm_dist_displace': norm_dist_displace})
        self.k_DISPLACE = []
        for norm_dist_displace_value in norm_dist_displace:
            self.k_DISPLACE.append(np.abs(norm_dist_displace_value))

    def create_initial_rateConst_influx(self):
        '''
            Method to create sets reaction rate constants for all influx reactions.
            k_IN (list): list of reaction rate constants for all influx reactions
        '''
        norm_dist_in = np.random.normal(loc=self.input_params['theta_in']['mean'], 
                                        scale=np.sqrt(self.input_params['theta_in']['variance']),
                                        size=self.nI)
        norm_dist_in = np.abs(norm_dist_in)
        self.input_params['theta_in'].update({'norm_dist_in': norm_dist_in})
        self.k_IN = []
        for norm_dist_in_value in norm_dist_in:
            self.k_IN.append(np.abs(norm_dist_in_value))

    def create_reaction_influx(self):
        '''
            Method to create a set of influx reactions
            R_IN (tuple): set of influx reactions
            nR_IN (int): number of influx reactions
        '''
        self.R_IN = []
        for i_i, i in enumerate(self.I):
            R_in = f'-> {i}()'
            self.R_IN.append(R_in)
        self.nR_IN = len(self.R_IN)

    def create_initial_rateConst_efflux(self):
        '''
            Method to create sets reaction rate constants for all efflux reactions.
            k_OUT (list): list of reaction rate constants for all efflux reactions
        '''
        norm_dist_out = np.random.normal(loc=self.input_params['theta_out']['mean'], 
                                         scale=np.sqrt(self.input_params['theta_out']['variance']),
                                         size=1)
        norm_dist_out = np.abs(norm_dist_out)
        self.input_params['theta_out'].update({'norm_dist_out': norm_dist_out})
        self.k_OUT = []
        for norm_dist_out_value in norm_dist_out:
            self.k_OUT.append(np.abs(norm_dist_out_value))

    def create_reaction_efflux(self):
        '''
            Method to create a set of efflux reactions
            R_OUT (tuple): set of efflux reactions
            nR_OUT (int): number of efflux reactions
        '''
        self.R_OUT = []
        for i_o, o in enumerate(self.O):
            R_out = f'{o}() ->'
            self.R_OUT.append(R_out)
        self.nR_OUT = len(self.R_OUT)

    def create_reaction_sets(self):
        '''
            R (tuple): set of all reactions
            nR (int): total number of reactions in network
        '''
        self.R = self.R_BIND + self.R_DISPLACE + self.R_IN + self.R_OUT
        self.nR = len(self.R)

    def create_reaction_lookup(self):
        '''
            Methods to create a lookup dictionary for sets of reactions and number of reactions
            reaction_lookup (dict): lookup dictionary for sets of reactions and number of reactions
        '''
        self.reaction_lookup = {
            'R_BIND': self.R_BIND,
            'R_DISPLACE': self.R_DISPLACE,
            'R_IN': self.R_IN,
            'R_OUT': self.R_OUT,
            'R': self.R,
            'nR_BIND': self.nR_BIND,
            'nR_DISPLACE': self.nR_DISPLACE,
            'nR_IN': self.nR_IN,
            'nR_OUT': self.nR_OUT,
            'nR': self.nR
        }

    def create_rateConst_lookup(self):
        '''
            Method to create a lookup dictionary for rate constants for each types of reactions.
            Returns:
                rateConst_lookup (dict): lookup dictionary for all rate constants of each reactions.
        '''
        self.rate_BIND = {}
        for i_b, r_bind in enumerate(self.R_BIND):
            self.rate_BIND.update({r_bind: [self.k_BIND[i_b]]})

        self.rate_DISPLACE = {}
        for i_d, r_displace in enumerate(self.R_DISPLACE):
            self.rate_DISPLACE.update({r_displace: [self.k_DISPLACE[i_d]]})

        self.rate_IN = {}
        for i_i, r_in in enumerate(self.R_IN):
            self.rate_IN.update({r_in: [self.k_IN[i_i]]})

        self.rate_OUT = {}
        for i_o, r_out in enumerate(self.R_OUT):
            self.rate_OUT.update({r_out: [self.k_OUT[0]]})

        self.rateConst_lookup = {
            'rate_BIND': self.rate_BIND,
            'rate_DISPLACE': self.rate_DISPLACE,
            'rate_IN': self.rate_IN,
            'rate_OUT': self.rate_OUT
        }

    def update_time_params(self):
        '''
            Method to calculate the time-related variables and update the time_params dictionary
        '''

        t_start = self.time_params['t_start']
        t_end = self.time_params['t_end']
        t_perturb = self.time_params['t_perturb']
        t_hold = self.time_params['t_hold']
        num_perturb = int(np.ceil((t_end - t_perturb) / t_hold))
        self.time_params.update({'num_perturb': num_perturb})
        time_array = list(np.arange(start = self.time_params['t_start'] + self.time_params['t_perturb'],
                                    stop  = self.time_params['t_end'],
                                    step  = self.time_params['t_hold']))
        time_array = tuple([t_start] + time_array + [t_end])
        if len(time_array) != self.time_params['num_perturb'] + 2: # t_start and t_end elements
            raise Exception('Length of time array does not agree with the number of perturbations')
        self.time_params.update({'time_array': time_array})

    def update_influx_rate(self):
        '''
            Method to generate the influx rate at each perturbation
        '''

        for influx_reaction, influx_rates in self.rateConst_lookup['rate_IN'].items():
            base_IN = influx_rates[0]
            for perturb_index in range(self.time_params['num_perturb']):
                base_IN *= (1 + np.random.rand() - 0.5)
                self.rateConst_lookup['rate_IN'][influx_reaction].append(base_IN)

        for influx_reaction in self.rateConst_lookup['rate_IN'].keys():
            if len(self.rateConst_lookup['rate_IN'][influx_reaction]) != self.time_params['num_perturb'] + 1:
                raise Exception('Length of influx vector does not agree with the number of perturbations')

    def create_perturbation_lookup(self):
        '''
            Method to create a dictionary showing time and influx rate at each perturbation event
        '''

        self.perturbation_lookup = {}
        for time_index, time in enumerate(self.time_params['time_array'][0:-1]):
            self.perturbation_lookup.update({time: {}})
            for r_in, rate_ins in self.rateConst_lookup['rate_IN'].items():
                self.perturbation_lookup[time].update({r_in: rate_ins[time_index]})


if __name__ == '__main__':
    from params_dsd import input_params, time_params
    import code

    dsd = RandomDSDChemistry(input_params=input_params, time_params=time_params)
    code.interact(local=locals())
