from random_dna_chem import RandomDNAStrandDisplacementCircuit
from params import input_params5, time_params1
from string import Template


def getDomain(Upper, Lower, Full, Partial):
    '''
        Convert the DNA Strand Displacement Circuit Strands to Domains
        Args:
            Upper (list): list of all upper strands
            Lower (list): list of all lower strands
            Full (list): list of all full double strands
            Partial (list): list of all partial double strands
        Returns:
            allDom (list): list of all recognition and toehold domains
            recognitionDom (list): list of all recognition domains (long enough to bind irreversibly to their complementary sequence)
            toeholdDom (list): list of all toehold domains (short enough to bind resersibly to their complementary sequence)
    '''

    allDom = []
    for u in Upper:
        allDom.append(u.lower())
        allDom.append('t' + u.lower())
    for l in Lower:
        allDom.append(l.lower())
        allDom.append('t' + l.lower())
    for f in Full:    
        allDom.remove('l' + f[-1])
        allDom.remove('tl' + f[-1])
    for p in Partial:
        try:
            allDom.remove('l' + p[-1])
        except ValueError:
            pass
    recognitionDom = []
    toeholdDom = []
    for dom in allDom:
        toeholdDom.append(dom) if 't' in dom else recognitionDom.append(dom)

    return allDom, recognitionDom, toeholdDom

def renderSpecies(species):
    '''
        Render the Species of the Visual DSD code
    '''
    rendered_species = []
    for s in species:
        rendered_species.append(f'{s}()')
    return rendered_species

def renderRate(toehold, base_bind, base_unbind):
    '''
        Render the Binding Rate of the Visual DSD code
    '''
    bindList = []
    unbindList = []
    for i in range(len(toehold)):
        bindList.append(f'k{i}')
        unbindList.append(f'u{i}')

    rendered_rates = []
    for k, u in zip(bindList, unbindList):
        rendered_rates.append(f'{k}={base_bind}; {u}={base_unbind}')

    return bindList, unbindList, rendered_rates

def renderDom(toehold, bindList, unbindList):
    '''
        Render the Domains definition 'dom' of the Visual DSD code
    '''
    rendered_domains = []
    for th, bind, unbind in zip(toehold, bindList, unbindList):
        rendered_domains.append(f'dom {th} = {{bind={bind}; unbind={unbind}}}')
    return rendered_domains

def renderDef(Upper, Lower, Full, Partial, recognition):
    '''
        Render the Modules definition 'def' of the Visual DSD code
    '''
    rendered_modules = []
    for u in Upper:
        dom = u.lower()
        rendered_modules.append(f'def {u}() = <{dom} t{dom}^>')
    for l in Lower:
        dom = l.lower()
        rendered_modules.append(f'def {l}() = {{{dom} t{dom}^}}')
    for f in Full:
        dom = f[:2].lower()
        rendered_modules.append(f'def {f}() = [{dom} t{dom}^]')
    for p in Partial:
        u_dom = p[:2].lower()
        l_dom = p[2:].lower()
        rendered_modules.append(f'def {p}() = [{u_dom}]<t{u_dom}^>{{t{l_dom}^}}')
    return rendered_modules

def renderInitialConditions(conU, conL, conF, conP):
    '''
        Render the Initial Conditions of the Visual DSD Code
    '''
    rendered_initial_conditions = []
    for u, conu in conU.items():
        rendered_initial_conditions.append(f'{conu} {u}()')
    for l, conl in conL.items():
        rendered_initial_conditions.append(f'{conl} {l}()')
    for f, conf in conF.items():
        rendered_initial_conditions.append(f'{conf} {f}()')
    for p, conp in conP.items():
        rendered_initial_conditions.append(f'{conp} {p}()')
    return rendered_initial_conditions

def renderDSD(rendered_species, rendered_rates, rendered_domains, rendered_modules, rendered_initial_conditions):
    data = {
        'final': randomDNAChem.time_params['t_end'],
        'species_list': '; '.join(rendered_species),
        'rate_list': '; '.join(rendered_rates),
        'dom_list': '\n'.join(rendered_domains),
        'def_list': '\n'.join(rendered_modules),
        'initial_cond_list': '\n| '.join(rendered_initial_conditions)
    }
    with open('template.txt', 'r') as f:
        src = Template(f.read())
    result = src.substitute(data)

    with open('rendered.txt', 'w') as f:
        f.write(result)


randomDNAChem = RandomDNAStrandDisplacementCircuit(input_params=input_params5, time_params=time_params1)
S = randomDNAChem.species_lookup['S']
U = randomDNAChem.species_lookup['U']
L = randomDNAChem.species_lookup['L']
F = randomDNAChem.species_lookup['F']
P = randomDNAChem.species_lookup['P']
conU = randomDNAChem.concentration_lookup['conU']
conL = randomDNAChem.concentration_lookup['conL']
conF = randomDNAChem.concentration_lookup['conF']
conP = randomDNAChem.concentration_lookup['conP']

allDom, recognitionDom, toeholdDom = getDomain(Upper=U, Lower=L, Full=F, Partial=P)
rendered_species = renderSpecies(species=S)
bindList, unbindList, rendered_rates = renderRate(toehold=toeholdDom, base_bind=0.003, base_unbind=0.1)
rendered_domains = renderDom(toehold=toeholdDom, bindList=bindList, unbindList=unbindList)
rendered_modules = renderDef(Upper=U, Lower=L, Full=F, Partial=P, recognition=recognitionDom)
rendered_initial_conditions = renderInitialConditions(conU=conU, conL=conL, conF=conF, conP=conP)
renderDSD(rendered_species, rendered_rates, rendered_domains, rendered_modules, rendered_initial_conditions)
