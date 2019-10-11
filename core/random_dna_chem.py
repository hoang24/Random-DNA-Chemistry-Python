import numpy as np
import time


class RandomDNAStrandDisplacementCircuit(object):
    '''
        Random DNA Strand Displacement Circuit class to generate a random DNA strand diaplacement circuit
        Take parameters (input, etc.) from params.py
            Optimal random DNA chemistry using results from optimization methods (turning point measure and genetic search)
            (nL, nU, k) = (3, 4, 3)
            Global ordered type
        Attributes:
            input_params (dict): lookup dictionary for input parameters
            time_params (dict): lookup dictionary for timing parameters
            species_lookup (dict): lookup dictionary for species names and number of species
            order_lookup (dict): lookup dictionary for orders of partial double strands
            reaction_lookup (dict): lookup dictionary for orders of partial double strands
            concentration_lookup (dict): lookup dictionary for all concentrations of all species
            rateConst_lookup (dict): lookup dictionary for reaction rate constants of all reactions
        Methods:
            create_species_lookup()
            impose_order_partial()
            create_reaction_lookup()
            create_concentration_lookup()
            create_rateConst_lookup()
            update_time_params()
    '''

    def __init__(self, input_params={}, time_params={}):
        '''
            Python initialization method.
            Args:
                input_params (dict): input parameters
                time_params (dict): timing parameters
        '''

        self.input_params = input_params.copy()
        self.time_params = time_params.copy()
        self.species_lookup = self.create_species_lookup()
        self.order_lookup = self.impose_order_partial()
        self.reaction_lookup = self.create_reaction_lookup()
        self.concentration_lookup = self.create_concentration_lookup()
        self.rateConst_lookup = self.create_rateConst_lookup()
        self.update_time_params()

    def _create_species_single(self):
        '''
            Method to initialize the set of single DNA strands
            Returns:
                nU (int): number of upper single strands
                nL (int): number of lower single strands
                U (tuple): set of upper single strands
                L (tuple): set of lower single strands
        '''

        nL = int(round(self.input_params['n'] / (1 + self.input_params['p'])))
        L = []
        for nl in range(nL):
            L.append('L{}'.format(nl))
        L = tuple(L)

        nU = self.input_params['n'] - nL
        U = []
        for nu in range(nU):
            U.append('U{}'.format(nu))
        U = tuple(U)

        return nU, nL, U, L

    def _create_species_double(self, nU, nL, U, L):
        '''
            Method to initialize the set of single DNA strands
            Args:
                nU (int): number of upper single strands
                nL (int): number of lower single strands
                U (tuple): set of upper single strands
                L (tuple): set of lower single strands
            Returns:
                nF (int): number of full double strands
                nP (int): number of partial double strands
                F (tuple): set of full double strands
                P (tuple): set of partial double strands
        '''

        min_single = np.min([nL, nU])
        nF = int(self.input_params['y'] * min_single)
        index_full_list = np.random.choice(a=min_single, size=nF, replace=False)
        F = []
        for index_full in index_full_list:
            F.append(U[index_full] + L[index_full])

        # positive normal distribution of partial double strands per upper strand
        norm_dist = np.random.normal(loc=self.input_params['phi']['mean'], 
                                     scale=self.input_params['phi']['variance'],
                                     size=nU)
        self.input_params['phi'].update({'norm_dist': norm_dist})

        num_double_per_upper = []
        for norm_dist_value in norm_dist:
            num_double_per_upper.append(int(round(norm_dist_value)))
        if len(set(num_double_per_upper)) == 1:
            num_double_per_upper = num_double_per_upper[0] # k (int): number of double strands per upper strands
        else:
            raise BaseException

        P = F
        while len(set(F + P)) != len(F + P):
            P = []
            for upper in U: # for each upper strand
                # Choose randomly number of lower strand counterparts without repetitions
                lower_per_upper = np.random.choice(a=L, size=num_double_per_upper, replace=False)
                for lower in lower_per_upper:
                    partial = upper + lower
                    P.append(partial)
        nP = len(P)
        F = tuple(F)
        P = tuple(P)

        return nF, nP, F, P

    def mirror_selection(self):
        '''
            Method to perform mirror selection
        '''
        pass

    def _create_all_species_sets(self, U, nU, L, nL, F, nF, P, nP):
        '''
            Method to initialize the all-species sets
            Args:
                nU (int): number of upper single strands
                nL (int): number of lower single strands
                U (tuple): set of upper single strands
                L (tuple): set of lower single strands
                nF (int): number of full double strands
                nP (int): number of partial double strands
                F (tuple): set of full double strands
                P (tuple): set of partial double strands
            Returns:
                nS (int): number of all DNA strands in network
                nSS (int): number of all single DNA strands in network
                nDS (int): number of all double DNA strands in network
                S (tuple): set of all DNA strands in network
                SS (tuple): set of all single DNA strands in network
                DS (tuple): set of all double DNA strands in network
        '''

        S = U + L + F + P
        nS = len(S)
        SS = U + L
        nSS = len(SS)
        DS = F + P
        nDS = len(DS)

        if nS != (nU + nL + nF + nP):
            raise BaseException
        if nSS != (nU + nL):
            raise BaseException
        if nDS != (nF + nP):
            raise BaseException

        return nS, nSS, nDS, S, SS, DS

    def _create_influx_efflux(self, S, nS):
        '''
            Methods to create set of influx and efflux
            Args:
                S (tuple): set of all species in network
                nS (int): number of all species in network
            Returns:
                nI (int): number of influx species
                nO (int): number of efflux species
                I (tuple): set of influx species
                O (tuple): set of efflux species
        '''

        nI = int(round(self.input_params['a_in'] * nS))
        nO = int(round(self.input_params['a_out'] * nS))
        I = np.random.choice(a=S, size=nI, replace=False)
        O = np.random.choice(a=S, size=nO, replace=False)

        return nI, nO, I, O

    def create_species_lookup(self):
        '''
            Methods to create a lookup dictionary for sets of species and number of species
            Returns:
                species_lookup (dict): lookup dictionary for sets of species and number of species
        '''
        nU, nL, U, L = self._create_species_single()
        nF, nP, F, P = self._create_species_double(nU=nU, nL=nL, U=U, L=L)
        nS, nSS, nDS, S, SS, DS = self._create_all_species_sets(U=U, nU=nU, L=L, nL=nL, F=F, nF=nF, P=P, nP=nP)
        nI, nO, I, O = self._create_influx_efflux(S=S, nS=nS)
        species_lookup = {
            'U': U,
            'L': L,
            'F': F,
            'P': P,
            'SS': SS,
            'DS': DS,
            'S': S,
            'I': I,
            'O': O,
            'nU': nU,
            'nL': nL,
            'nF': nF,
            'nP': nP,
            'nSS': nSS,
            'nDS': nDS,
            'nS': nS,
            'nI': nI,
            'nO': nO
        }

        return species_lookup

    def impose_order_partial(self):
        '''
            Method to create a lookup dictionary to show ordering of the double strands
            Returns:
                order_lookup (dict): lookup dictionary for the orders of the partial double strands
        '''

        order_lookup = {}
        for p_index in range(self.species_lookup['nP']):
            partial_with_order = np.random.choice(a=self.species_lookup['P'])
            while partial_with_order in order_lookup.values():
                partial_with_order = np.random.choice(a=self.species_lookup['P'])
            order_lookup.update({p_index: partial_with_order})

        for f_index, f in enumerate(self.species_lookup['F']):
            order_lookup.update({self.species_lookup['nP'] + f_index: f})

        return order_lookup

    def _get_key_by_value_from_dict(self, dictionary, value):
        '''
            Method to retrieve a key given its value in a dictionary. Apply to dictionary with non-repeating values.
            Args:
                dictionary (dict): a dictionary
                value (int, float, string, etc): value in dictionary
            Returns:
                key (string): key of the value in the dictionary
        '''
        item_list = dictionary.items()
        key = -1
        for item in item_list:
            if value == item[1]:
                key = item[0]
        if key == -1:
            raise Exception('Value given does not exist in given dictionary. Value = {}'.format(value))

        return key

    def _create_reaction_binding(self, U, L, DS):
        '''
            Method to create a set of binding reactions
            Args:
                U (tuple): set of upper strand species
                L (tuple): set of lower strand species
                DS (tuple): set of all double strand species
            Returns:
                R_BIND (tuple): set of binding reactions
                nR_BIND (int): number of binding reactions
        '''

        R_BIND = []
        for u in U:
            for l in L:
                ds_result = u + l
                if ds_result in DS:
                    R_bind = u + ' + ' + l + ' --> ' + ds_result
                    R_BIND.append(R_bind)
        R_BIND = tuple(R_BIND)
        nR_BIND = len(R_BIND)

        return R_BIND, nR_BIND

    def _create_reaction_displacement(self, SS):
        '''
            Method to create a set of displacement reactions
            Args:
                SS (tuple): set of all single strand species
            Returns:
                R_DISPLACE (tuple): set of displacement reactions
                nR_DISPLACE (int): number of displacement reactions
        '''

        R_DISPLACE = []
        for reactant_ss in SS: # for each single strand
            for reactant_ds_order, reactant_ds in self.order_lookup.items(): # for each double strand in order
                if reactant_ss not in reactant_ds: # if single strand is not part of double strand
                    if 'U' in reactant_ss: # if single strand is upper strand
                        result_ds = reactant_ss + reactant_ds[-2:]
                        result_ss = reactant_ds[:2]
                    elif 'L' in reactant_ss: # if single strand is lower strand
                        result_ds = reactant_ds[:2] + reactant_ss
                        result_ss = reactant_ds[-2:]
                    else: # invalid single strand notation
                        raise Exception('Unknown DNA species type. Defined type: U for Upper strand, L for Lower strand.')
                    if result_ds in self.order_lookup.values(): # if the resulted double strand made out of the single strand and the reactant double strand exists
                        result_ds_order = self._get_key_by_value_from_dict(self.order_lookup, result_ds)
                        if result_ds_order > reactant_ds_order: # if the resulted double strand has higher order than the reactant double strand
                            R_displace = reactant_ss + ' + ' + reactant_ds + ' --> ' + result_ds + ' + ' + result_ss
                            R_DISPLACE.append(R_displace)
        R_DISPLACE = tuple(R_DISPLACE)
        nR_DISPLACE = len(R_DISPLACE)

        return R_DISPLACE, nR_DISPLACE

    def _create_reaction_influx(self, I):
        '''
            Method to create a set of influx reactions
            Args:
                I (tuple): set of influx species
            Returns:
                R_IN (tuple): set of influx reactions
                nR_IN (int): number of influx reactions
        '''
        R_IN = []
        for i in I:
            R_in = '0 --> ' + i
            R_IN.append(R_in)
        R_IN = tuple(R_IN)
        nR_IN = len(R_IN)

        return R_IN, nR_IN

    def _create_reaction_efflux(self, O):
        '''
            Method to create a set of efflux reactions
            Args:
                O (tuple): set of efflux species
            Returns:
                R_OUT (tuple): set of efflux reactions
                nR_OUT (int): number of efflux reactions
        '''
        R_OUT = []
        for o in O:
            R_out = o + ' --> 0'
            R_OUT.append(R_out)
        R_OUT = tuple(R_OUT)
        nR_OUT = len(R_OUT)

        return R_OUT, nR_OUT

    def _create_all_reaction_set(self, R_BIND, R_DISPLACE, R_IN, R_OUT, nR_BIND, nR_DISPLACE, nR_IN, nR_OUT):
        '''
            Method to create set of all reactions
            Args:
                R_BIND (tuple): set of binding reactions
                nR_BIND (int): number of binding reactions
                R_DISPLACE (tuple): set of displacement reactions
                nR_DISPLACE (int): number of displacement reactions
                R_IN (tuple): set of influx reactions
                nR_IN (int): number of influx reactions
                R_OUT (tuple): set of efflux reactions
                nR_OUT (int): number of efflux reactions
            Returns:
                R (tuple): set of all reactions
                nR (int): total number of reactions in network
        '''

        R = tuple(R_BIND + R_DISPLACE + R_IN + R_OUT)
        nR = len(R)
        if nR != (nR_BIND + nR_DISPLACE + nR_IN + nR_OUT):
            raise BaseException

        return R, nR

    def create_reaction_lookup(self):
        '''
            Methods to create a lookup dictionary for sets of reactions and number of reactions
            Returns:
                reaction_lookup (dict): lookup dictionary for sets of reactions and number of reactions
        '''

        R_DISPLACE, nR_DISPLACE = self._create_reaction_displacement(SS=self.species_lookup['SS'])
        R_BIND, nR_BIND = self._create_reaction_binding(U=self.species_lookup['U'], L=self.species_lookup['L'], DS=self.species_lookup['DS'])
        R_IN, nR_IN = self._create_reaction_influx(I=self.species_lookup['I'])
        R_OUT, nR_OUT = self._create_reaction_efflux(O=self.species_lookup['O'])
        R, nR = self._create_all_reaction_set(R_BIND=R_BIND, R_DISPLACE=R_DISPLACE, R_IN=R_IN, R_OUT=R_OUT, 
                                              nR_BIND=nR_BIND, nR_DISPLACE=nR_DISPLACE, nR_IN=nR_IN, nR_OUT=nR_OUT)
        reaction_lookup = {
            'R_DISPLACE': R_DISPLACE,
            'R_BIND': R_BIND,
            'R_IN': R_IN,
            'R_OUT': R_OUT,
            'R': R,
            'nR_DISPLACE': nR_DISPLACE,
            'nR_BIND': nR_BIND,
            'nR_IN': nR_IN,
            'nR_OUT': nR_OUT,
            'nR': nR
        }

        return reaction_lookup

    def _create_initial_concentration(self):
        '''
            Method to create sets of initial concentration for each DNA strand types.
            Returns:
                ic_U (tuple): initial concentrations for all upper single strands
                ic_L (tuple): initial concentrations for all lower single strands
                ic_F (tuple): initial concentrations for all full double strands
                ic_P (tuple): initial concentration for all partial double strands
        '''

        ic_U = tuple(np.random.randint(0, 1000, self.species_lookup['nU']))
        ic_L = tuple(np.random.randint(0, 1000, self.species_lookup['nL']))
        ic_F = tuple(np.random.randint(0, 1000, self.species_lookup['nF']))
        ic_P = tuple(np.random.randint(0, 1000, self.species_lookup['nP']))

        return ic_U, ic_L, ic_F, ic_P

    def create_concentration_lookup(self):
        '''
            Method to create a lookup dictionary for all concentrations for each species
            Returns:
                initial_concentration_lookup (dict): lookup dictionary for initial concentration of each species
                concentration_lookup (dict): lookup dictionary for all concentrations of each species
        '''

        ic_U, ic_L, ic_F, ic_P = self._create_initial_concentration()

        conU = {}
        for u_index, u in enumerate(self.species_lookup.get('U')):
            conU.update({'{}'.format(u): [ic_U[u_index]]})

        conL = {}
        for l_index, l in enumerate(self.species_lookup.get('L')):
            conL.update({'{}'.format(l): [ic_L[l_index]]})

        conF = {}
        for f_index, f in enumerate(self.species_lookup.get('F')):
            conF.update({'{}'.format(f): [ic_F[f_index]]})

        conP = {}
        for p_index, p in enumerate(self.species_lookup.get('P')):
            conP.update({'{}'.format(p): [ic_P[p_index]]})

        concentration_lookup = {
            'conU': conU,
            'conL': conL,
            'conF': conF,
            'conP': conP
        }

        return concentration_lookup

    def _create_initial_rateConst_binding(self, nR_BIND):
        '''
            Method to create sets reaction rate constants for all binding reactions.
            Args:
                nR_BIND (int): number of binding reactions
            Returns:
                k_BIND (tuple): list of reaction rate constants for all binding reactions
        '''

        norm_dist_bind = np.random.normal(loc=self.input_params['theta']['mean'], 
                                          scale=np.sqrt(self.input_params['theta']['variance']),
                                          size=nR_BIND)
        self.input_params['theta'].update({'norm_dist_bind': norm_dist_bind})
        k_BIND = []
        for norm_dist_bind_value in norm_dist_bind:
            k_BIND.append(np.abs(norm_dist_bind_value))
        k_BIND = tuple(k_BIND)

        return k_BIND

    def _create_initial_rateConst_displacement(self, nR_DISPLACE):
        '''
            Method to create sets reaction rate constants for all displacement reactions.
            Args:
                nR_DISPLACE (int): number of displacement reactions
            Returns:
                k_DISPLACE (tuple): list of reaction rate constants for all displacement reactions
        '''

        norm_dist_displace = np.random.normal(loc=self.input_params['theta']['mean'], 
                                              scale=np.sqrt(self.input_params['theta']['variance']),
                                              size=nR_DISPLACE)
        self.input_params['theta'].update({'norm_dist_displace': norm_dist_displace})
        k_DISPLACE = []
        for norm_dist_displace_value in norm_dist_displace:
            k_DISPLACE.append(np.abs(norm_dist_displace_value))
        k_DISPLACE = tuple(k_DISPLACE)

        return k_DISPLACE

    def _create_initial_rateConst_influx(self, nR_IN):
        '''
            Method to create sets reaction rate constants for all influx reactions.
            Args:
                nR_IN (int): number of influx reactions
            Returns:
                k_IN (list): list of reaction rate constants for all influx reactions
        '''

        norm_dist_in = np.random.normal(loc=self.input_params['theta_in']['mean'], 
                                        scale=np.sqrt(self.input_params['theta_in']['variance']),
                                        size=nR_IN)
        self.input_params['theta_in'].update({'norm_dist_in': norm_dist_in})
        k_IN = []
        for norm_dist_in_value in norm_dist_in:
            k_IN.append(np.abs(norm_dist_in_value))
        k_IN = tuple(k_IN)

        return k_IN

    def _create_initial_rateConst_efflux(self, nR_OUT):
        '''
            Method to create sets reaction rate constants for all efflux reactions.
            Args:
                nR_OUT (int): number of efflux reactions
            Returns:
                k_OUT (list): list of reaction rate constants for all efflux reactions
        '''

        norm_dist_out = np.random.normal(loc=self.input_params['theta_out']['mean'], 
                                         scale=np.sqrt(self.input_params['theta_out']['variance']),
                                         size=nR_OUT)
        self.input_params['theta_out'].update({'norm_dist_out': norm_dist_out})
        k_OUT = []
        for norm_dist_out_value in norm_dist_out:
            k_OUT.append(np.abs(norm_dist_out_value))
        k_OUT = tuple(k_OUT)

        return k_OUT

    def create_rateConst_lookup(self):
        '''
            Method to create a lookup dictionary for rate constants for each types of reactions.
            Returns:
                rateConst_lookup (dict): lookup dictionary for all rate constants of each reactions.
        '''

        k_BIND = self._create_initial_rateConst_binding(nR_BIND=self.reaction_lookup['nR_BIND'])
        k_DISPLACE = self._create_initial_rateConst_displacement(nR_DISPLACE=self.reaction_lookup['nR_DISPLACE'])
        k_IN = self._create_initial_rateConst_influx(nR_IN=self.reaction_lookup['nR_IN'])
        k_OUT = self._create_initial_rateConst_efflux(nR_OUT=self.reaction_lookup['nR_OUT'])

        rate_BIND = {}
        for bind_index, r_bind in enumerate(self.reaction_lookup.get('R_BIND')):
            rate_BIND.update({'{}'.format(r_bind): [k_BIND[bind_index]]})

        rate_DISPLACE = {}
        for displace_index, r_displace in enumerate(self.reaction_lookup.get('R_DISPLACE')):
            rate_DISPLACE.update({'{}'.format(r_displace): [k_DISPLACE[displace_index]]})

        rate_IN = {}
        for in_index, r_in in enumerate(self.reaction_lookup.get('R_IN')):
            rate_IN.update({'{}'.format(r_in): [k_IN[in_index]]})

        rate_OUT = {}
        for out_index, r_out in enumerate(self.reaction_lookup.get('R_OUT')):
            rate_OUT.update({'{}'.format(r_out): [k_OUT[out_index]]})

        rateConst_lookup = {
            'rate_BIND': rate_BIND,
            'rate_DISPLACE': rate_DISPLACE,
            'rate_IN': rate_IN,
            'rate_OUT': rate_OUT
        }

        return rateConst_lookup

    def update_time_params(self):
        '''
            Method to calculate the time-related variables and update the time_params dictionary

        '''

        t_end = self.time_params['t_end']
        t_perturb = self.time_params['t_perturb']
        t_hold = self.time_params['t_hold']
        num_perturb = int(np.ceil((t_end - t_perturb) / t_hold))
        self.time_params.update({'num_perturb': num_perturb})

    def run_chem(self):
        '''
            Method to run the random DNA strand displacement circuit chemistry
            Args:
                
            Returns:
                
        '''

        t = self.time_params['t_start']
        iteration = 0
        while t < self.time_params['t_end']:
            pass

    def Gillespie_initialization(self):
        '''
            Method for Gillespie algorithm's initialization step
            Args:
                
            Returns:
                
        '''

        pass

    def Gilespie_monte_carlo(self):
        '''
            Method for Gillespie algorithm's Monte Carlo step
            Args:
                
            Returns:
                
        '''

        pass

    def Gillespie_update(self):
        '''
            Method for Gillespie algorithm's update step
            Args:
                
            Returns:
                
        '''

        pass

    def Gillespie_iterate(self):
        '''
            Method for Gillespie algorithm's iterate step
            Args:
                
            Returns:
                
        '''

        pass
