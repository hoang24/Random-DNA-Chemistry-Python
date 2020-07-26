from random_dna_chem import RandomDNAStrandDisplacementCircuit
from params_dsd import input_params, time_params
from string import Template


def get_toeholds(uppers, lowers):
    '''
        Extract toehold domains from Random DNA Strand Circuit species
        Args:
            uppers (tuple of str): list of all upper strands
            lowers (tuple of str): list of all lower strands
        Returns:
            toeholds (list of str): list of toehold domains
    '''
    toeholds = []
    for upper in uppers:
        toehold = 't' + upper[-1]
        toeholds.append(toehold)
    for lower in lowers:
        toehold = 't' + lower[-1]
        if toehold not in toeholds:
            toeholds.append(toehold)
    return toeholds

def render_species(species):
    '''
        Render the species from Random DNA Strand Circuit species
        Args:
            species (list of str): all species in the Random DNA Strand Circuit
        Returns:
            rendered_species (list of str): rendered visualDSD species in the S() format
    '''
    rendered_species = []
    for s in species:
        rendered_species.append(f'{s}()')
    return rendered_species

def render_rates(toeholds, base_bind, base_unbind):
    '''
        Render the binding and unbinding rate names and values
        Args:
            toeholds (list of str): list of toehold domains
            base_bind (float): base value of binding rate
            base_unbind (float): base value of unbinding rate
        Returns:
            bind_list (list of str): list of all binding rate names (k0, k1, k2, etc.)
            unbind_list (list of str): list of all unbinding rate names (u0, u1, u2, etc.)
            rendered_rates (list of str): list of all bind/unbind rate names and values (k0=..., u0=..., k1=..., u1=..., etc.)
    '''
    bind_list = []
    unbind_list = []
    rendered_rates = []
    for toehold in toeholds:
        bind_list.append(f'k{toehold[-1]}')
        unbind_list.append(f'u{toehold[-1]}')
    for bind, unbind in zip(bind_list, unbind_list):
        rendered_rates.append(f'{bind}={base_bind}; {unbind}={base_unbind}')
    return bind_list, unbind_list, rendered_rates

def render_doms(toeholds, bind_list, unbind_list):
    '''
        Render the Domains definition 'dom' in the Visual DSD code
        Args:
            toeholds (list of str): list of toehold domains
            bind_list (list of str): list of all binding rate names (k0, k1, k2, etc.)
            unbind_list (list of str): list of all unbinding rate names (u0, u1, u2, etc.)
        Returns:
            rendered_domains (list of str): list of domain definitions for each species (dom t0 ..., dom t1..., etc.)
    '''
    rendered_domains = []
    for toehold, bind, unbind in zip(toeholds, bind_list, unbind_list):
        rendered_domains.append(f'dom {toehold} = {{bind={bind}; unbind={unbind}}}')
    return rendered_domains

def render_singles(uppers, lowers, fulls):
    '''
        Render the Modules definition 'def' for the single strands (upper and lower) in the Visual DSD code
        Args:
            uppers (tuple of str): list of all upper strands
            lowers (tuple of str): list of all lower strands
            fulls (tuple of str): list of all full double strands
        Returns:
            rendered_singles (list of str): list of upper and lower strand modules
    '''
    uppers = list(uppers)
    lowers = list(lowers)
    rendered_singles = []
    transitions = ['tA', 'tB']
    # Render complementary strands
    val = 0
    for full in fulls:
        upper = full[:2]
        lower = full[2:]
        rendered_singles.append(f'def {upper}() = < main^ trans^ {transitions[val]}^ t{upper[-1]}^ >')
        rendered_singles.append(f'def {lower}() = {{ main^* trans^* {transitions[val]}^* t{lower[-1]}^* }}')
        val = (val + 1) % 2
    # Render non-complementary strands
    for full in fulls:
        uppers.remove(full[:2])
        lowers.remove(full[2:])
    val = 0
    for upper in uppers:
        rendered_singles.append(f'def {upper}() = < main^ trans^ {transitions[val]}^ t{upper[-1]}^ >')
        val = (val + 1) % 2
    val = 1
    for lower in lowers:
        rendered_singles.append(f'def {lower}() = {{ main^* trans^* {transitions[val]}^* t{lower[-1]}^* }}')
        val = (val + 1) % 2
    return rendered_singles

def render_doubles(rendered_singles):
    '''
        Render the Modules definition 'def' for the double strands (full and partial) in the Visual DSD code
        Args:
            uppers (tuple of str): list of all upper strands
            lowers (tuple of str): list of all lower strands
            fulls (tuple of str): list of all full double strands
            partials (list of str): list of all partial double strands
            rendered_singles (list of str): list of upper and lower strand modules
        Returns:
            rendered_doubles (list of str): list of partial and full double strand modules
            rendered_singles (list of str): list of upper and lower strand modules
    '''
    rendered_uppers = []
    rendered_lowers = []
    for rendered_single in rendered_singles:
        if 'def U' in rendered_single:
            rendered_uppers.append(rendered_single)
        elif 'def L' in rendered_single:
            rendered_lowers.append(rendered_single)
        else:
            raise Exception('Undefined definition for Visual DSD modules')

    rendered_doubles = []
    for rendered_upper in rendered_uppers:
        for rendered_lower in rendered_lowers:
            upper_name = rendered_upper.split(' ')[1][:2]
            lower_name = rendered_lower.split(' ')[1][:2]
            complex1 = ''
            complex2 = ''
            upper1 = ''
            upper2 = ''
            lower1 = ''
            lower2 = ''

            # render tA and tB toehold domains
            if rendered_upper.split(' ')[6] in rendered_lower.split(' ')[6]:
                complex1 = rendered_upper.split(' ')[6]
            else:
                upper1 = rendered_upper.split(' ')[6]
                lower1 = rendered_lower.split(' ')[6]

            # render t0, t1, t2, etc. toehold domains
            if rendered_upper.split(' ')[7] in rendered_lower.split(' ')[7]:
                if rendered_upper.split(' ')[6] in rendered_lower.split(' ')[6]:
                    complex2 = rendered_upper.split(' ')[7]
                else:
                    upper2 = rendered_upper.split(' ')[7]
                    lower2 = rendered_lower.split(' ')[7]
            else:
                upper2 = rendered_upper.split(' ')[7]
                lower2 = rendered_lower.split(' ')[7]

            rendered_doubles.append(f'def {upper_name}{lower_name}() = [main^ trans^ {complex1} {complex2}]<{upper1} {upper2}>{{{lower1} {lower2}}}')

    # clean up the rendered double strand strings
    for double_index in range(len(rendered_doubles)):
        rendered_doubles[double_index] = rendered_doubles[double_index].replace('< >{ }', '')
        rendered_doubles[double_index] = rendered_doubles[double_index].replace('  ', ' ')
    for double_index, rendered_double in enumerate(rendered_doubles):
        rendered_doubles[double_index] = rendered_doubles[double_index].replace(' ]', ']')
        rendered_doubles[double_index] = rendered_doubles[double_index].replace('< ', '<')
        rendered_doubles[double_index] = rendered_doubles[double_index].replace(' >', '>')
        rendered_doubles[double_index] = rendered_doubles[double_index].replace('{ ', '{')
        rendered_doubles[double_index] = rendered_doubles[double_index].replace(' }', '}')

    # clean up the rendered single strand strings
    for single_index in range(len(rendered_singles)):
        rendered_singles[single_index] = rendered_singles[single_index].replace('< ', '<')
        rendered_singles[single_index] = rendered_singles[single_index].replace(' >', '>')
        rendered_singles[single_index] = rendered_singles[single_index].replace('{ ', '{')
        rendered_singles[single_index] = rendered_singles[single_index].replace(' }', '}')

    return rendered_doubles, rendered_singles

def render_initial_conditions(initial_upper, initial_lower, initial_full, initial_partial):
    '''
        Render the Initial Conditions of the Visual DSD Code
        Args:
            initial_upper (dict of int): initial number of species for the upper strands
            initial_lower (dict of int): initial number of species for the lower strands
            initial_full (dict of int): initial number of species for the full double strands
            initial_partial (dict of int): initial number of species for the partial double strands
        Returns:
            rendered_initial_conditions (list of str): list of the rendered DSD initial conditions
    '''
    rendered_initial_conditions = []
    for upper, conu in initial_upper.items():
        rendered_initial_conditions.append(f'{conu[0]} {upper}()')
    for lower, conl in initial_lower.items():
        rendered_initial_conditions.append(f'{conl[0]} {lower}()')
    for full, conf in initial_full.items():
        rendered_initial_conditions.append(f'{conf[0]} {full}()')
    for partial, conp in initial_partial.items():
        rendered_initial_conditions.append(f'{conp[0]} {partial}()')
    return rendered_initial_conditions

def render_reactions(bind_reactions, displace_reactions):
    '''
        Render the Visual DSD reactions and rates.
        Args:
            bind_reactions (dict of str): for binding reactions, keys are reactions, values are rates
            displace_reactions (dict of str): for displacement reactions, keys are reactions, values are rates
        Returns:
            rendered_reactions (list of str): list of rendered DSD reactions and rates
    '''
    rendered_reactions = []
    for r_bind, rate_bind in bind_reactions.items():
        rendered_reactions.append(r_bind.replace('->', f'->{{{rate_bind[0]}}}'))
    for r_displace, rate_displace in displace_reactions.items():
        rendered_reactions.append(r_displace.replace('->', f'->{{{rate_displace[0]}}}'))
    return rendered_reactions

def render_DSD_program(final, rendered_species, rendered_rates, rendered_domains, rendered_singles, 
                       rendered_doubles, rendered_initial_conditions, rendered_reactions):
    '''
        Render the visual DSD program
        Args:
            final (float): the end time of the simulation
            rendered_species (list of str): rendered visualDSD species in the S() format
            rendered_rates (list of str): list of all bind/unbind rate names and values (k0=..., u0=..., k1=..., u1=..., etc.)
            rendered_domains (list of str): list of domain definitions for each species (dom t0 ..., dom t1..., etc.)
            rendered_singles (list of str): list of upper and lower strand modules
            rendered_doubles (list of str): list of partial and full double strand modules            
            rendered_initial_conditions (list of str): list of the rendered DSD initial conditions
    '''
    rendered_modules = rendered_singles + rendered_doubles # concatenate
    data_DSD = {
        'final': final,
        'species_list': '; '.join(rendered_species),
        'rate_list': '; '.join(rendered_rates),
        'dom_list': '\n'.join(rendered_domains),
        'def_list': '\n'.join(rendered_modules),
        'initial_cond_list': '\n| '.join(rendered_initial_conditions),
        'reaction_list': '\n| '.join(rendered_reactions)
    }
    with open('template.txt', 'r') as template_file:
        template_DSD = Template(template_file.read())
    rendered_DSD = template_DSD.substitute(data_DSD)
    with open('rendered.txt', 'w') as rendered_file:
        rendered_file.write(rendered_DSD)


if __name__ == '__main__':
    randomDNAChem = RandomDNAStrandDisplacementCircuit(input_params=input_params, time_params=time_params)
    final = randomDNAChem.time_params['t_end']
    S = randomDNAChem.species_lookup['S']
    U = randomDNAChem.species_lookup['U']
    L = randomDNAChem.species_lookup['L']
    F = randomDNAChem.species_lookup['F']
    P = randomDNAChem.species_lookup['P']
    conU = randomDNAChem.concentration_lookup['conU']
    conL = randomDNAChem.concentration_lookup['conL']
    conF = randomDNAChem.concentration_lookup['conF']
    conP = randomDNAChem.concentration_lookup['conP']
    bind_reactions = randomDNAChem.rateConst_lookup['rate_BIND']
    displace_reactions = randomDNAChem.rateConst_lookup['rate_DISPLACE']

    toeholds = get_toeholds(uppers=U, lowers=L)
    rendered_species = render_species(species=S)
    bind_list, unbind_list, rendered_rates = render_rates(toeholds=toeholds, base_bind=0.003, base_unbind=0.1)
    rendered_domains = render_doms(toeholds=toeholds, bind_list=bind_list, unbind_list=unbind_list)
    rendered_singles = render_singles(uppers=U, lowers=L, fulls=F)
    rendered_doubles, rendered_singles = render_doubles(rendered_singles=rendered_singles)
    rendered_initial_conditions = render_initial_conditions(initial_upper=conU, initial_lower=conL, 
                                                            initial_full=conF, initial_partial=conP)
    rendered_reactions = render_reactions(bind_reactions=bind_reactions, displace_reactions=displace_reactions)
    render_DSD_program(final, rendered_species, rendered_rates, rendered_domains, 
                       rendered_singles, rendered_doubles, rendered_initial_conditions, rendered_reactions)
