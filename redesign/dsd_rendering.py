from string import Template


class DSDRenderProgram():
    def __init__(self, chemistry, initial_conditions, influx_index, dsd_filename, initial, final):
        '''
            chemistry (class): random DNA strand displacement circuit chemistry class
            initial_conditions (dict of int): initial number of species for a simulation period
            influx_index (int): index of perturbation period
            initial (float): start time of the simulation period
            final (float): end time of the simulation period
            filename (str): location and name of the rendered Visual DSD file
        '''
        # Constant parameters
        self.bind_rate = 0.003 # base value of binding rate
        self.unbind_rate = 0.1 # base value of unbinding rate
        self.points = 1001 # number of datapoints in the simulation period

        # User-defined parameters
        self.chemistry = chemistry
        self.initial_conditions = initial_conditions
        self.influx_index = influx_index
        self.dsd_filename = dsd_filename
        self.initial = initial
        self.final = final

        # Rendering methods
        self.render_species()
        self.render_domains()
        self.render_species_def_single()
        self.render_species_def_double()
        self.render_initial_conditions()
        self.render_reactions()
        self.render_DSD_program()

    def render_species(self):
        '''
            Render the species from Random DNA Strand Circuit species
            Args:
                species (tuple of str): all species in the Random DNA Strand Circuit
            Returns:
                rendered_species (list of str): rendered visualDSD species in the S() format
        '''
        self.rendered_species = []
        for s in self.chemistry.S:
            self.rendered_species.append(f'{s}()')

    def render_domains(self):
        '''
            Render the Domains definition 'dom' in the Visual DSD code
            rendered_domains (list of str): list of domain definitions for each species (dom t0 ..., dom t1..., etc.)
        '''
        self.rendered_domains = []

        self.rendered_domains.append(f'dom M = {{bind={self.bind_rate}; unbind={self.unbind_rate}; colour="black"}} // main domain')

        for switching_domain in self.chemistry.switching_domains:
            self.rendered_domains.append(f'dom {switching_domain} = {{bind={self.bind_rate}; unbind={self.unbind_rate}}} // switching domain')

        for identity_domain in self.chemistry.identity_domains:
            self.rendered_domains.append(f'dom t{identity_domain} = {{bind={self.bind_rate}; unbind={self.unbind_rate}}} // identity domain')

    def render_species_def_single(self):
        '''
            Render the Modules definition 'def' for the single strands (upper and lower) in the Visual DSD code
            rendered_single (list of str): list of upper and lower strand modules
        '''
        self.rendered_single = []
        for u in self.chemistry.U:
            doms = self.chemistry.domain_lookup['U'][u]
            self.rendered_single.append(f'def {u}() = <{doms[0]}^ {doms[1]}^ {doms[2]}^ t{doms[3]}^> // upper strand')

        for l in self.chemistry.L:        
            doms = self.chemistry.domain_lookup['L'][l]
            self.rendered_single.append(f'def {l}() = {{{doms[0]}^* {doms[1]}^* {doms[2]}^* t{doms[3]}^*}} // lower strand')

    def render_species_def_double(self):
        '''
            Render the Modules definition 'def' for the double strands (full and partial) in the Visual DSD code
            rendered_doubles (list of str): list of partial and full double strand modules
            rendered_singles (list of str): list of upper and lower strand modules
        '''
        self.rendered_double = []
        for f in self.chemistry.F:
            u, l = self.chemistry.decode_double(f)
            doms_u = self.chemistry.domain_lookup['U'][u]
            doms_l = self.chemistry.domain_lookup['L'][l]
            if doms_u != doms_l:
                raise Exception(f'{f} is not a full double strand')
            self.rendered_double.append(f'def {f}() = [{doms_u[0]}^ {doms_u[1]}^ {doms_u[2]}^ t{doms_u[3]}^] // full double strand')

        for ps in self.chemistry.PS:
            u, l = self.chemistry.decode_double(ps)
            doms_u = self.chemistry.domain_lookup['U'][u]
            doms_l = self.chemistry.domain_lookup['L'][l]
            if doms_u[2] != doms_l[2]:
                raise Exception(f'{ps} is not a strong partial double strand')
            self.rendered_double.append(f'def {ps}() = [{doms_u[0]}^ {doms_u[1]}^ {doms_u[2]}^]<t{doms_u[3]}^>{{t{doms_l[3]}^*}} // strong (2/3) partial double strand')

        for pw in self.chemistry.PW:
            u, l = self.chemistry.decode_double(pw)
            doms_u = self.chemistry.domain_lookup['U'][u]
            doms_l = self.chemistry.domain_lookup['L'][l]
            self.rendered_double.append(f'def {pw}() = [{doms_u[0]}^ {doms_u[1]}^]<{doms_u[2]}^ t{doms_u[3]}^>{{{doms_l[2]}^* t{doms_l[3]}^*}} // weak (1/3) partial double strand')

    def render_initial_conditions(self):
        '''
            Render the Initial Conditions of the Visual DSD Code
            rendered_initial_conditions (list of str): list of the rendered DSD initial conditions
        '''
        conU = self.initial_conditions['conU']
        conL = self.initial_conditions['conL']
        conF = self.initial_conditions['conF']
        conPS = self.initial_conditions['conPS']
        conPW = self.initial_conditions['conPW']

        self.rendered_initial_conditions = []
        for u, conu in conU.items():
            self.rendered_initial_conditions.append(f'{conu[0]} {u}()')
        for l, conl in conL.items():
            self.rendered_initial_conditions.append(f'{conl[0]} {l}()')
        for f, conf in conF.items():
            self.rendered_initial_conditions.append(f'{conf[0]} {f}()')
        for ps, conps in conPS.items():
            self.rendered_initial_conditions.append(f'{conps[0]} {ps}()')
        for pw, conpw in conPW.items():
            self.rendered_initial_conditions.append(f'{conpw[0]} {pw}()')

    def render_reactions(self):
        '''
            Render the Visual DSD reactions and rates.
            rendered_reactions (list of str): list of rendered DSD reactions and rates
        '''

        bind_rates = self.chemistry.rateConst_lookup['rate_BIND']
        displace_rates = self.chemistry.rateConst_lookup['rate_DISPLACE']
        influx_rates = self.chemistry.rateConst_lookup['rate_IN']
        efflux_rates = self.chemistry.rateConst_lookup['rate_OUT']

        self.rendered_reactions = []
        for r_bind, rate_bind in bind_rates.items():
            self.rendered_reactions.append(r_bind.replace('->', f'->{{{rate_bind[0]}}}') + ' // binding reaction')
        for r_displace, rate_displace in displace_rates.items():
            self.rendered_reactions.append(r_displace.replace('->', f'->{{{rate_displace[0]}}}') + ' // displacement reaction')
        for r_in, rate_in in influx_rates.items():
            self.rendered_reactions.append(r_in.replace('->', f'->{{{rate_in[self.influx_index]}}}') + ' // influx reaction')
        for r_out, rate_out in efflux_rates.items():
            self.rendered_reactions.append(r_out.replace('->', f'->{{{rate_out[0]}}}') + ' // efflux reaction')

    def render_DSD_program(self):
        '''
            Render the visual DSD program
        '''
        self.rendered_modules = self.rendered_single + self.rendered_double
        data_DSD = {
            'initial': self.initial,
            'final': self.final,
            'points': self.points,
            'species_list': '; '.join(self.rendered_species),
            'dom_list': '\n'.join(self.rendered_domains),
            'def_list': '\n'.join(self.rendered_modules),
            'initial_cond_list': '\n| '.join(self.rendered_initial_conditions),
            'reaction_list': '\n| '.join(self.rendered_reactions)
        }

        with open('dsd_template.txt', 'r') as template_file:
            template_DSD = Template(template_file.read())
        self.rendered_DSD = template_DSD.substitute(data_DSD)

        with open(self.dsd_filename, 'w') as rendered_file:
            rendered_file.write(self.rendered_DSD)

if __name__ == '__main__':
    from params import input_params, time_params
    from random_dsd_chem import RandomDSDChemistry
    import code

    dsd = RandomDSDChemistry(input_params=input_params, time_params=time_params)
    DSDRenderProgram(chemistry=dsd, initial_conditions=dsd.species_count_lookup, influx_index=0, dsd_filename='rendered.txt', initial=0, final=0.01)

    code.interact(local=locals())
