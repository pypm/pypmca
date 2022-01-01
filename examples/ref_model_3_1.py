# -*- coding: utf-8 -*-
"""
pyPM.ca reference model #3.1

To make the critical parameters less dependent of each other, break the measured
populations into separate paths.

Reduced the lag time between contagious and reported, (as compared to reference model #1)
to better match the situation in BC and elsewhere - transitions now occur where they should

Add a contact tracing population, to allow another mechanism to remove people from
the contagious population early. Initial value has p=0, so it does not act until it
is switched on.

Version 2.1: - turn on negative binomial for multiplier. Add parameter to control dispersion (set to 0.999: tiny
effect, at p = 0.5, variance is double poisson) - turn on reporting errors for reported. Add parameter to control
level (set to 1: no effect) - add capability to include reporting anomalies: new pop = tested positive receives
contact and symp, reduces contagious. added to reported new pop = reporting anomalies: propagates to reported. Add
new transition: injector to reporting anomalies

Version 2.2:
 - add 3 additional rate transitions
 - add an additional outbreak
 - add an addition reporting anomaly
 - include parameter to specify low edge of fraction of report backlog that is cleared
 - include parameter to specify days of week which have no reporting (cases are included the next day)

Version 2.3:
 - separate the time scales for contagious_removal from the reporting process: have a separate propagator from
   contagious to removed (removed from circulation), which subtracts from contagious. This allows an independent
   measurement of this time scale from the turn-around process
 - the removal delay default is mean, sigma = 7, 3
 - change default latent period (delay from infected to contagious) was (mean,sigma = 2,1) now: (5,3)
 - change default icu_fraction to 0. (since many jurisdictions have no ICU data)
 - change default delay for non_icu: was (mean, sigma = 12,3) now: (6,3)
 - change default delay for icu: was (14,1) now: (4,2)
 - change default alpha_0 and alpha_1-n: was (0.385, 0.062) to: (0.4, 0.1)
 - change default transition time for alpha_0 -> alpha_1 from 16 to 20
 - change default cont_0 from 55. to 10.
 - change initial boot population from 1. to 0.1
 - change default death delay sigma from 5. to 10.

Version 2.5:
 - add additional transitions
 - implement large scale vaccination, with limited supply:
    - immunization occurs for a fraction of those vaccinated (vaccine effectiveness)
    - only vaccinate the population not having a positive test result and not already vaccinated:
      "vaccination candidates" (or "vaccan" for short)
        - reduce vaccan for each vaccination
        - reduce vaccan for those reported who were vaccination candidates N(rep-vac-can):
          those reported who were NOT vaccination candidates are those who were vaccinated but did not gain immunity,
          N(rep-not-vac-can) = N(rep) * vac/tot * (1-vaccine effectiveness)
          N(rep-vac-can) = N(rep) * (1 - (1-ve)*vac/tot)
    - the size of the susceptible population is reduced as a result of both vaccinations and infections
        - to calculate the reduction from vaccination: must track the sub-population of vaccine candidates
          who are susceptible (or "susvaccan" for short). The expected reduction in susceptible population
          depends on the ratio of the sizes of the "susvaccan" and "vaccan" populations and the vaccine effectiveness
        - each infection reduces the susceptible population by 1. The expected reduction of "susvaccan" would be
          the ratio of the sizes of the susvaccan and susceptible populations

Version 2.6:
 - move most delay distributions from 'norm' to 'gamma'
 - add linear transition for reporting fraction
 - add reporting noise for deaths and hospitalization
 - set death reporting noise to be weekly

Version 2.7:
 - more transitions
 - change linear modification to be start_value end_value, rather than start_value slope

Version 2.8:
 - more testing fraction transitions
 - add second infection cycle for B117 variant

Version 2.9:
 - more vaccination transitions
 - more variant alpha transitions
 - add third infection cycle for B.1.617.2 variant (use the letter w in place of v for naming)
 - add second outbreaks for the two variants
 - include waning immunity (separately for vaccination and natural):
   - add natural immunized population
   - propagate (with long time delays) from immunized and naturally immunized to "lost immunity" that
     then increases susceptible

Version 3.0:
 - reported does not necessarily lead to removal from vaccine candidates
   (new parameter, frac_report_novacc, has default=1 to be backward compatible)
 - more wariant transitions (delta)
 - more vaccine transitions
 - 4th strain, for omicron: w-> x
   - since natural and vaccine escape is suspected, special susceptible populations are setup for those
 - allow breakthrough populations to be tracked (reported, hospitalized, deaths)
   - breakthrough is any person vaccinated who is not immunized (waning parameters determine susceptible size)
   - share the alphas of the naive populations (adjust waning parameters to incorporate altered immunity)
   - the original susceptible includes all. bt_susceptible produces bt_infected, which propagate to both bt_contagious
     and contagious
 - include boosters

 Version 3.1:
  - To allow for different hospitalization fractions etc, use independent populations
  - Use the new collector connector, to use with collector populations (such as reported, susceptible, etc)
  - Use the operator to combine parameters for calculating fractions etc
  - This requires pypmca version 0.2.24 or later
  - remove all contact tracing and all references to ventilation (reduce number of populations to track)
  - by changing the version number to 4, the breakthrough cycle is removed (speed up processing)

@author: karlen
"""

from pypmca import Model, Population, Delay, Parameter, Multiplier, Propagator, \
    Splitter, Adder, Subtractor, Modifier, Injector, Collector, Operator

# This is not listed in package requirements:
# Users generally do not need to run this code.
import matplotlib
import numpy as np


def complementary_color(c):
    my_hex = matplotlib.colors.to_hex(c)
    rgb = (my_hex[1:3], my_hex[3:5], my_hex[5:7])
    comp = ['%02X' % (255 - int(a, 16)) for a in rgb]
    return '#' + ''.join(comp)


# define rotated color: i_rot is in range(n_rot)
# Currently 4 replications of populations:
# 0:collector, 1: bt (breakthrough), 2: ve (vaccine escape), 3: ne (natural escape), 4: os (original susceptible)
n_rot = 5


def rotated_color(i_rot, c):
    cos = np.cos(2. * np.pi * i_rot / n_rot)
    sin = np.sin(2. * np.pi * i_rot / n_rot)
    term = (1. - cos) / 3.
    sqrt = np.sqrt(1. / 3.)
    rm = []

    row = [cos + term, term - sqrt * sin, term + sqrt * sin]
    rm.append(row)
    row = [term + sqrt * sin, cos + term, term - sqrt * sin]
    rm.append(row)
    row = [term - sqrt * sin, term + sqrt * sin, cos + term]
    rm.append(row)

    my_hex = matplotlib.colors.to_hex(c)
    rgb = (my_hex[1:3], my_hex[3:5], my_hex[5:7])
    rgb_vec = [int(a, 16) for a in rgb]

    rot_vec = []
    for irow in range(3):
        tot = 0.
        for icol in range(3):
            tot += rgb_vec[icol] * rm[irow][icol]
        if tot < 0:
            tot = 0
        if tot > 255:
            tot = 255
        tot = int(tot + 0.5)

        rot_vec.append('%02X' % tot)

    return '#' + ''.join(rot_vec)


# Reference model for BC

version = 4
subversion = 1

# no_bt turns off the breakthrough. Turned off in version 4.
no_bt = version == 4

bc_model = Model('ref_model_' + str(version) + '_' + str(subversion))
bc_model.set_t0(2020, 3, 1)

fast_delay = Delay('fast', 'fast', model=bc_model)

initial_pop_par = Parameter('N_0', 5000000, 5000, 50000000,
                            'Population of the region at t0', 'int')

total_pop = Population('total', initial_pop_par,
                       'total population of the region', color='black')

initial_contagious_par = Parameter('cont_0', 10., 0., 5000.,
                                   'Number of contagious people at t0',
                                   hidden=False)

# Initialization of independent populations and their collectors
# **** SETUP THE INITIAL VALUES OF THE COLLECTORS TO MATCH THE
# THE SUM OF THE INITIAL VALUES OF THE POPULATIONS ****
# ==============================================================

if not no_bt:
    cycles = {'os': 'original susceptible', 'bt': 'breakthrough', 've': 'vaccine escape', 'ne': 'natural escape'}
else:
    # removes the breakthrough cycle
    cycles = {'os': 'original susceptible', 've': 'vaccine escape', 'ne': 'natural escape'}

variants = {'o': 'original', 'v': 'variant', 'w': 'wariant', 'x': 'xariant'}
# color rotation
irc = {'os': 1, 'bt': 2, 've': 3, 'ne': 4}

variants_in_cycle = {'os': ['o', 'v', 'w', 'x'],
                     'bt': ['o', 'v', 'w', 'x'],
                     've': ['x'],
                     'ne': ['x']}

collector_populations = {}

# susceptibles
color = 'cornflowerblue'

susceptible_pops = {}
for cycle in cycles:
    susceptible_pops[cycle] = Population(cycle + '_susceptible', 0, 'Susceptible:' + cycles[cycle],
                                         color=rotated_color(irc[cycle], color))

susceptible_pops['os'].set_initial_value(initial_pop_par)

susceptible_pop = Population('susceptible', initial_pop_par, 'total number of people who could become infected',
                             color=color)
collector_populations[susceptible_pop] = susceptible_pops

# infected
colors = ['orange', 'darkgoldenrod', 'sienna', 'sandybrown']

infected_pops = {}
for variant, color in zip(variants, colors):
    infected_pops_variant = {}
    for cycle in cycles:
        if variant in variants_in_cycle[cycle]:
            infected_pops_variant[cycle] = Population(cycle + '_infected_' + variant, 0,
                                                      cycles[cycle] + ' infected by type: ' + variants[variant],
                                                      color=rotated_color(irc[cycle], color))
    infected_pops[variant] = infected_pops_variant

    infected_pop = Population('infected_' + variant, 0, 'infected by type: ' + variants[variant],
                              hidden=True, color=color)
    collector_populations[infected_pop] = infected_pops_variant

# contagious
colors = ['red', 'rosybrown', 'maroon', 'lightcoral']

contagious_pops = {}
contagious_by_variant = {}
for variant, color in zip(variants, colors):
    contagious_pops_variant = {}
    for cycle in cycles:
        if variant in variants_in_cycle[cycle]:
            contagious_pops_variant[cycle] = Population(cycle + '_contagious_' + variant, 0,
                                                        cycles[cycle] + ' contagious with type: ' + variants[variant],
                                                        color=rotated_color(irc[cycle], color))
    contagious_pops[variant] = contagious_pops_variant

    # if variant == 'o':
    #    contagious_pops_variant['os'].set_initial_value(initial_contagious_par)

    contagious_pop = Population('contagious_' + variant, 0,
                                'contagious with type: ' + variants[variant], hidden=True, color=color)
    # done below: collector_populations[contagious_pop] = contagious_pops_variant
    contagious_by_variant[variant] = contagious_pop

contagious_by_variant['o'].set_initial_value(initial_contagious_par)
original_contagious_pop = contagious_by_variant['o']

# oooooooooooooooooooooooooooo
# Define the infection cycles
# oooooooooooooooooooooooooooo

neg_binom_par = Parameter('neg_binom_p', 0.5, 0.001, 0.999,
                          'Dispersion parameter p for neg binom')

contagious_frac = Parameter('cont_frac', 0.9, 0., 1.,
                            'fraction of infected people that become contagious',
                            hidden=False)
contagious_delay_pars = {
    'mean': Parameter('cont_delay_mean', 5., 0., 50.,
                      'mean time from being infected to becoming contagious'),
    'sigma': Parameter('cont_delay_sigma', 3., 0.01, 20.,
                       'standard deviation of times from being infected to becoming contagious')
}

contagious_delay = Delay('cont_delay', 'norm', contagious_delay_pars, bc_model)

trans_rates = {}
for variant in variants:
    par_name = 'alpha'
    if variant != 'o':
        par_name = 'alpha_' + variant

    trans_rate = Parameter(par_name, 0.4, 0., 2., 'mean number of people that a ' + variants[variant] +
                           ' contagious person infects per day', hidden=True)
    trans_rates[variant] = trans_rate

    for cycle in cycles:
        if variant in variants_in_cycle[cycle]:
            bc_model.add_connector(
                Multiplier(cycle + '_infection cycle_' + variant,
                           [susceptible_pops[cycle], contagious_by_variant[variant], total_pop],
                           infected_pops[variant][cycle], trans_rate, fast_delay, bc_model,
                           distribution='nbinom', nbinom_par=neg_binom_par))

            bc_model.add_connector(
                Propagator(cycle + ' infected to contagious' + variant, infected_pops[variant][cycle],
                           contagious_pops[variant][cycle], contagious_frac, contagious_delay))

# --------------------------------------------
# keep track of naturally immunized population - a fraction of those infected gain immunity
# --------------------------------------------

nat_immunity_fraction = Parameter('nat_immunity_frac', 1., 0., 1.,
                                  'fraction of infections that lead to natural immunity')

color = 'gold'
nat_immunized_pops = {}
for cycle in cycles:
    nat_immunized_pops[cycle] = Population(cycle + '_nat_immunized', 0,
                                           'people who gained immunity through ' + cycles[cycle] + ' infection',
                                           hidden=True, color=rotated_color(irc[cycle], color))

    for variant in variants:
        if variant in variants_in_cycle[cycle]:
            # all strains of infection lead to natural immunity (that eventually wanes)
            bc_model.add_connector(
                Adder(cycle + ' natural immunity ' + variant, infected_pops[variant][cycle],
                      nat_immunized_pops[cycle], scale_factor=nat_immunity_fraction))

            # remove the same fraction from the corresponding susceptible population
            bc_model.add_connector(
                Subtractor('subtract infected from susceptible: ' + cycle + '_' + variant, susceptible_pops[cycle],
                           infected_pops[variant][cycle]))

nat_immunized_pop = Population('nat_immunized', 0,
                               'people who gained immunity through infection',
                               hidden=True, color=color)
collector_populations[nat_immunized_pop] = nat_immunized_pops

# a fraction of all strain infections will be susceptible to natural escape (except xariant)
# ------------------------------------------------------------------------------------------

nat_escape_fraction = Parameter('ne_frac', 0.5, 0., 1.,
                                'fraction of natural immunizations that can escape')
if not no_bt:
    cycle_list = ['os', 'bt']
else:
    cycle_list = ['os']

for cycle in cycle_list:
    for variant in variants:
        if variant != 'x':
            bc_model.add_connector(
                Propagator(cycle + ' natural escape ' + variant, infected_pops[variant][cycle],
                           susceptible_pops['ne'], nat_escape_fraction, fast_delay))

# waning immunity
# ---------------

# Specified by a delay distribution F(t) and a fraction of the immunizations that wane, f.
# Waning is not considered to be equivalent to returning to a naive state, but instead
# the immunity is not as effective after waning.
# The product f*F(t) is the cumulative probability of waning immunity, which accounts
# for the combined effect of:
# - increasing probability for an individual to wane
# - the reducing effectiveness of the immunity of waned individuals

nat_waned_delay_pars = {
    'mean': Parameter('nat_waned_delay_mean', 365., 1., 10000., 'mean time from natural immunity to waned immunity'),
    'sigma': Parameter('nat_waned_delay_sigma', 50., 1., 500.,
                       'standard deviation of times from natural immunity to waned immunity')
}
nat_waned_delay = Delay('nat_waned_delay', 'gamma', nat_waned_delay_pars, bc_model)

nat_waned_fraction = Parameter('nat_waned_frac', 1., 0., 1.,
                               'fraction of natural immunizations that wane')

color = 'mediumvioletred'

for cycle in cycle_list:
    waned_nat_immunity_pop = Population(cycle + ' waned nat immunity', 0,
                                        'people who lost nat immunity some time after gaining nat immunity: ' +
                                        cycles[cycle], hidden=True, color=rotated_color(irc[cycle], color))

    bc_model.add_connector(
        Propagator(cycle + ' nat waned immunity', nat_immunized_pops[cycle], waned_nat_immunity_pop,
                   nat_waned_fraction, nat_waned_delay))

    # Once waned, the susceptible population increases

    bc_model.add_connector(
        Adder(cycle + ' waned nat to susceptible', waned_nat_immunity_pop, susceptible_pops[cycle]))

# Include vaccination
# -------------------

# A multiplier is used. Once the vaccination candidate population goes below zero, it stops.

daily_vaccinated_pop = Population('daily vaccinated', 0, 'number of people vaccinated each day')

vaccinated_pop = Population('vaccinated', 0,
                            'people vaccinated', color='skyblue')
vaccan_pop = Population('vacc cand', initial_pop_par,
                        'vaccination candidates', color='brown')
susvaccan_pop = Population('sus vacc cand', initial_pop_par,
                           'vaccination candidates who are susceptible', color='seagreen')

vac_fraction = Parameter('vac_frac', 1., 0., 1.,
                         'fraction of vaccinations included == 1')

bc_model.add_connector(
    Multiplier('vaccination', [vaccan_pop, daily_vaccinated_pop, vaccan_pop],
               vaccinated_pop, vac_fraction, fast_delay, bc_model))

# Only unvaccinated people are vaccinated, a fraction of those were susceptible (designate those as usefully vaccinated)

usefully_vaccinated_pop = Population('usefully vaccinated', 0,
                                     'people who were vaccinated when susceptible', color='navy')

bc_model.add_connector(
    Adder('vaccinating susceptibles', vaccinated_pop, usefully_vaccinated_pop,
          ratio_populations=[susvaccan_pop, vaccan_pop]))

# All the susceptibles who were vaccinated become breakthrough candidates

if not no_bt:
    bt_candidate_pop = Population('bt_candidate', 0,
                                  'number of susceptible people who became eligible as a bt candidate')

    bc_model.add_connector(
        Propagator('breakthrough queue', usefully_vaccinated_pop, bt_candidate_pop, vac_fraction, fast_delay))

    # This means moving susceptible to bt susceptible

    bc_model.add_connector(
        Adder('bt_candidate adds to bt susceptible', bt_candidate_pop, susceptible_pops['bt']))

    bc_model.add_connector(
        Subtractor('bt_candidate removes susceptible', susceptible_pops['os'], bt_candidate_pop))

# Now immunize some of those susceptible people who were vaccinated

immunized_pop = Population('immunized', 0, 'number of susceptible people who were immunized by vaccine')

vaccine_effectiveness = Parameter('vaccine_eff', 0.8, 0., 1.,
                                  'probability that a susceptible person gains immunity when vaccinated')

immunized_delay_pars = {
    'mean': Parameter('immunized_delay_mean', 32., 5., 100., 'mean time from vaccination to immunity'),
    'sigma': Parameter('immunized_delay_sigma', 10., 1., 50.,
                       'standard deviation of times from vaccination to immunity')
}
immunized_delay = Delay('immunized_delay', 'gamma', immunized_delay_pars, bc_model)

bc_model.add_connector(
    Propagator('immunization', usefully_vaccinated_pop, immunized_pop, vaccine_effectiveness, immunized_delay))

# immunizations reduce the breakthrough susceptible (or original susceptible)
if not no_bt:
    bc_model.add_connector(
        Subtractor('subtract immunized from bt susceptible', susceptible_pops['bt'], immunized_pop))
else:
    bc_model.add_connector(
        Subtractor('subtract immunized from os susceptible', susceptible_pops['os'], immunized_pop))

# a fraction of immunizations will be susceptible to vaccine escape

vac_escape_fraction = Parameter('ve_frac', 0.5, 0., 1.,
                                'fraction of vaccine immunizations that can escape')

bc_model.add_connector(
    Propagator('vaccine escape', immunized_pop, susceptible_pops['ve'], vac_escape_fraction, fast_delay))

# Waning vaccination immunity

waned_vac_immunity_pop = Population('waned vac immunity', 0,
                                    'people who lost vac immunity some time after gaining vac immunity',
                                    hidden=True, color='mediumpurple')

vac_waned_delay_pars = {
    'mean': Parameter('vac_waned_delay_mean', 365., 1., 10000.,
                      'mean time from vaccination immunity to waned immunity'),
    'sigma': Parameter('vac_waned_delay_sigma', 50., 1., 500.,
                       'standard deviation of times from vaccination immunity to waned immunity')
}
vac_waned_delay = Delay('vac_waned_delay', 'gamma', vac_waned_delay_pars, bc_model)

vac_waned_fraction = Parameter('vac_waned_frac', 1., 0., 1.,
                               'fraction of vaccination immunizations that wane')

bc_model.add_connector(
    Propagator('vac waned immunity', immunized_pop, waned_vac_immunity_pop, vac_waned_fraction, vac_waned_delay))

# Once waned, the bt susceptible population increases and the ve susceptible population decreases

if not no_bt:
    bc_model.add_connector(
        Adder('waned to bt susceptible', waned_vac_immunity_pop, susceptible_pops['bt']))
else:
    bc_model.add_connector(
        Adder('waned to os susceptible', waned_vac_immunity_pop, susceptible_pops['os']))

bc_model.add_connector(
    Subtractor('waned removes ve susceptible', susceptible_pops['ve'], waned_vac_immunity_pop,
               scale_factor=vac_escape_fraction))

#####################
# Include boosters:
#####################
# Grow the population of booster candidates (delay following first dose)
# This requires the breakthrough process to be activated

if not no_bt:
    boostcan_pop = Population('boost cand', 0,
                              'booster candidates', color=rotated_color(1, 'brown'))

    boosting_delay_pars = {
        'mean': Parameter('boosting_delay_mean', 200., 120., 300., 'mean time from vaccination to booster eligibility'),
        'sigma': Parameter('boosting_delay_sigma', 20., 1., 200.,
                           'standard deviation of times from vaccination to booster eligibility')
    }
    boosting_delay = Delay('boosting_delay', 'gamma', boosting_delay_pars, bc_model)

    boost_fraction = Parameter('boost_frac', 1., 0., 1.,
                               'fraction of boosters included == 1')

    bc_model.add_connector(
        Propagator('boosting eligibility', vaccinated_pop, boostcan_pop, boost_fraction, boosting_delay))

    # A multiplier is used. Once the booster candidate population goes below zero, it stops.

    daily_boosted_pop = Population('daily boosted', 0, 'number of people boosted each day')

    boosted_pop = Population('boosted', 0, 'people boosted', color=rotated_color(1, 'skyblue'))

    bc_model.add_connector(
        Multiplier('boosting', [boostcan_pop, daily_boosted_pop, boostcan_pop],
                   boosted_pop, boost_fraction, fast_delay, bc_model))

    # Only a fraction of those were susceptible (designate those as usefully boosted)

    usefully_boosted_pop = Population('usefully boosted', 0,
                                      'people who were boosted when susceptible', color=rotated_color(1, 'navy'))

    bc_model.add_connector(
        Adder('boosting bt susceptibles', boosted_pop, usefully_boosted_pop,
              ratio_populations=[susceptible_pops['bt'], boostcan_pop]))

    # Now re-immunize some of those susceptible people who were boosted

    reimmunized_pop = Population('reimmunized', 0, 'number of bt susceptible people who were immunized by booster')
    booster_effectiveness = Parameter('booster_eff', 0.8, 0., 1.,
                                      'probability that a bt susceptible person gains immunity when boosted')

    bc_model.add_connector(
        Propagator('reimmunization', usefully_boosted_pop, reimmunized_pop, booster_effectiveness, immunized_delay))

    # a fraction of re-immunizations will be susceptible to vaccine escape

    bc_model.add_connector(
        Propagator('vaccine escape reimmunized', reimmunized_pop, susceptible_pops['ve'], vac_escape_fraction,
                   fast_delay))

    # Waning booster immunity

    waned_boost_immunity_pop = Population('waned booster immunity', 0,
                                          'people who lost immunity some time after gaining immunity through a booster',
                                          hidden=True, color=rotated_color(1, 'mediumpurple'))

    bc_model.add_connector(
        Propagator('boost waned immunity', reimmunized_pop, waned_boost_immunity_pop, vac_waned_fraction,
                   vac_waned_delay))

    # Once waned, the bt susceptible population increases and ve susceptible decreases

    bc_model.add_connector(
        Adder('booster waned to bt susceptible', waned_boost_immunity_pop, susceptible_pops['bt']))

    bc_model.add_connector(
        Subtractor('booster waned reduces ve susceptible', susceptible_pops['ve'], waned_boost_immunity_pop,
                   scale_factor=vac_escape_fraction))

# The contagious either recover or die
# This split is only used to track the deaths.
# The removal of recoveries/positive tests etc done in a separate path.
# Interest in health outcomes for naive and immunized treated separately
# Track deaths (and hosp) likewise.
# oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo

# The contagious_by_variant populations need to be collected before the propagation to deaths

for variant in variants:
    contagious_pop_variant = contagious_by_variant[variant]
    populations = []
    for key in contagious_pops[variant]:
        populations.append(contagious_pops[variant][key])

    bc_model.add_connector(
        Collector(contagious_pop_variant.name, populations, contagious_pop_variant))

# add death reporting noise that is reported out weekly
# to increase weekly death reporting variance, set death_backlog near 0, and set death_noise below 1.
death_noise_par = Parameter('death_noise', 1., 0., 1.,
                            'Death report noise parameter')
death_backlog_par = Parameter('death_backlog', 0., 0., 1.,
                              'Death report backlog parameter')
death_report_days = Parameter('death_report_days', 127, 0, 127, 'days of week with reporting (bit encoded)',
                              parameter_type='int')

recovered_pop = Population('recovered', 0,
                           'People who have recovered from the illness ', color='limegreen')

recover_fraction = Parameter('recover_frac', 0.99, 0., 1.,
                             'fraction of (o,v,w) infected who recover', hidden=False)

recover_delay_pars = {
    'mean': Parameter('recover_delay_mean', 14., 0., 50., 'mean time from infection to recovery'),
    'sigma': Parameter('recover_delay_sigma', 5., 0.01, 20.,
                       'standard deviation of times from infection to recovery')
}
recover_delay = Delay('recover_delay', 'gamma', recover_delay_pars, bc_model)

death_delay_pars = {
    'mean': Parameter('death_delay_mean', 21., 0., 50.,
                      'mean time from infection to death', hidden=False),
    'sigma': Parameter('death_delay_sigma', 10., 0.01, 20.,
                       'standard deviation of times from infection to death')
}
death_delay = Delay('death_delay', 'gamma', death_delay_pars, bc_model)

colors = ['mediumpurple', 'darkslateblue', 'mediumslateblue', 'darkviolet']
deaths_pops = {}
for variant, color in zip(variants, colors):
    if variant == 'x':
        deaths_pops_x = {}
        for cycle in cycles:
            deaths_pops_x[cycle] = Population(cycle + ' deaths ' + variant, 0,
                                              cycles[cycle] + ' who have died from ' + variant, hidden=True,
                                              color=rotated_color(irc[cycle], color), show_sim=True,
                                              report_noise=True, report_noise_par=death_noise_par,
                                              report_backlog_par=death_backlog_par, report_days=death_report_days,
                                              report_noise_weekly=True)
            # include these in the overall deaths collector
            deaths_pops[variant + cycle] = deaths_pops_x[cycle]

            recover_fraction_x = Parameter(cycle + '_recover_frac_' + variant, 0.99, 0., 1.,
                                           cycles[cycle] + ' fraction of x infected who recover', hidden=True)

            bc_model.add_connector(
                Splitter(cycle + ' recovery ' + variant, contagious_pops[variant][cycle],
                         [recovered_pop, deaths_pops_x[cycle]],
                         [recover_fraction_x], [recover_delay, death_delay]))

        deaths_pop_x = Population('deaths ' + variant, 0,
                                  'all those who have died from ' + variant, hidden=True,
                                  color=color, show_sim=True,
                                  report_noise=True, report_noise_par=death_noise_par,
                                  report_backlog_par=death_backlog_par, report_days=death_report_days,
                                  report_noise_weekly=True
                                  )
        collector_populations[deaths_pop_x] = deaths_pops_x
    else:
        deaths_pops[variant] = Population('deaths ' + variant, 0,
                                          'all those who have died from ' + variant, hidden=True,
                                          color=rotated_color(0, color), show_sim=True,
                                          report_noise=True, report_noise_par=death_noise_par,
                                          report_backlog_par=death_backlog_par, report_days=death_report_days,
                                          report_noise_weekly=True)

        bc_model.add_connector(
            Splitter('recovery_' + variant, contagious_by_variant[variant], [recovered_pop, deaths_pops[variant]],
                     [recover_fraction], [recover_delay, death_delay]))

deaths_pop = Population('deaths', 0,
                        'all those who have died', hidden=False,
                        color='indigo', show_sim=True,
                        report_noise=True, report_noise_par=death_noise_par,
                        report_backlog_par=death_backlog_par, report_days=death_report_days,
                        report_noise_weekly=True)

collector_populations[deaths_pop] = deaths_pops

###############################
# removal from contagious pools
###############################

# The removal of people from the contagious population is treated separately from
# reporting, to allow the reported number to be independent of delay in removal

removed_frac = Parameter('removed_frac', 1., 0., 1.,
                         'fraction of contagious people eventually removed == 1')
removed_delay_pars = {
    'mean': Parameter('removed_delay_mean', 7., 0., 20.,
                      'mean time from becoming contagious to being removed from contagious', hidden=False),
    'sigma': Parameter('removed_delay_sigma', 3., 0.1, 20.,
                       'std dev of times from becoming contagious to being removed from contagious',
                       hidden=False)
}
removed_delay = Delay('removed_delay', 'norm', removed_delay_pars, bc_model)

removed_pops = {}
color = ['blue', 'darkslateblue', 'dodgerblue', 'lightskyblue']
for variant, color in zip(variants, colors):
    removed_pop_variant = Population('removed_' + variant, 0,
                                     'People removed from the contagious population', color=color)
    for cycle in cycles:
        if variant in variants_in_cycle[cycle]:
            bc_model.add_connector(
                Propagator(cycle + ' contagious to removed: ' + variant, contagious_pops[variant][cycle],
                           removed_pop_variant, removed_frac, removed_delay))

    # Note that there is no removal from os, bt, ve and ne contagious... those populations are not part of the
    # infection cycle. As a result, the bt (ve and ne) contagious will be cumulative, not current.
    bc_model.add_connector(
        Subtractor('removal from contagious ' + variant, contagious_by_variant[variant], removed_pop_variant))

######################################
# Symptoms, testing, hospitalizations
######################################

# The newly contagious are split into two groups:
# symptomatic and non-symptomatic.

asymptomatic_recovered_pop = Population('asymptomatic recovered', 0,
                                        'People who have not shown symptoms', color='silver')

symptomatic_delay_pars = {
    'mean': Parameter('symptomatic_delay_mean', 2., 0., 50.,
                      'mean time from becoming contagious to having symptoms'),
    'sigma': Parameter('symptomatic_delay_sigma', 1., 0.01, 20.,
                       'std dev of times from becoming contagious to having symptoms')
}
symptomatic_delay = Delay('symptomatic_delay', 'gamma', symptomatic_delay_pars,
                          bc_model)

asymptomatic_delay_pars = {
    'mean': Parameter('asymp_rec_delay_mean', 12., 0., 50.,
                      'mean time from becoming contagious to recovery (without symptoms)'),
    'sigma': Parameter('asymp_rec_delay_sigma', 4., 0.01, 20.,
                       'std dev of times from becoming contagious to recovery with no symptoms')
}
asymptomatic_delay = Delay('asymp_rec_delay', 'gamma', asymptomatic_delay_pars,
                           bc_model)

# The nominal fraction who become symptomatic: 'os' and 'o'
symptomatic_fraction = Parameter('symptomatic_frac', 0.9, 0., 1.,
                                 'nominal fraction of contagious people who become symptomatic', hidden=False)

# Allow the fraction who become symptomatic to depend on cycles and variants

symptomatic_scale_cycle = {}
for cycle in cycles:
    if cycle != 'os':
        symptomatic_scale_cycle[cycle] = Parameter(cycle + '_symptomatic_scale', 1., 0., 1., cycles[cycle] +
                                                   ' scale factor for fraction of contagious people who become ' +
                                                   'symptomatic', hidden=True)

symptomatic_scale_variant = {}
for variant in variants:
    if variant != 'o':
        symptomatic_scale_variant[variant] = Parameter('symptomatic_scale_' + variant, 1., 0., 1., variants[variant] +
                                                       ' scale factor for fraction of contagious people who become ' +
                                                       'symptomatic', hidden=True)

colors = ['chocolate', 'peru', 'burlywood', 'blanchedalmond']

symptomatic_pops = {}
for variant, color in zip(variants, colors):
    symptomatic_pops_variant = {}
    for cycle in cycles:
        if variant in variants_in_cycle[cycle]:
            symptomatic_pops_variant[cycle] = Population(cycle + ' symptomatic ' + variant, 0, cycles[cycle] + ' ' +
                                                         variants[variant] + ' with symptoms',
                                                         color=rotated_color(irc[cycle], color))

            # scale the symptomatic fraction as needed
            if variant == 'o':
                if cycle == 'os':
                    symptomatic_frac = symptomatic_fraction
                else:
                    symptomatic_frac = Operator([symptomatic_fraction, symptomatic_scale_cycle[cycle]], '*')
            else:
                if cycle == 'os':
                    symptomatic_frac = Operator([symptomatic_fraction, symptomatic_scale_variant[variant]], '*')
                else:
                    symptomatic_frac = Operator([symptomatic_fraction, symptomatic_scale_cycle[cycle],
                                                 symptomatic_scale_variant[variant]], '**')

            bc_model.add_connector(
                Splitter(cycle + 'symptoms_' + variant, contagious_pops[variant][cycle],
                         [symptomatic_pops_variant[cycle], asymptomatic_recovered_pop],
                         [symptomatic_frac],
                         [symptomatic_delay, asymptomatic_delay]))

    symptomatic_pops[variant] = symptomatic_pops_variant

###########
# reporting
###########

# propagate a fraction of symptomatics to population with positive test results
# the reported population includes the reporting anomalies

report_noise_par = Parameter('report_noise', 1., 0., 1., 'Report noise parameter')

report_backlog_par = Parameter('report_backlog', 1., 0., 1., 'Report backlog parameter')

report_days = Parameter('report_days', 127, 0, 127, 'days of week with reporting (bit encoded)',
                        parameter_type='int')

reported_fraction = Parameter('reported_frac', 0.8, 0., 1.,
                              'fraction of symptomatic people who will receive a positive report')
reported_delay_pars = {
    'mean': Parameter('reported_delay_mean', 3., 0., 50.,
                      'mean time from becoming symptomatic to getting positive report', hidden=False),
    'sigma': Parameter('reported_delay_sigma', 1., 0.01, 20.,
                       'standard deviation of times from having symptoms to getting positive report')
}
reported_delay = Delay('reported_delay', 'gamma', reported_delay_pars, bc_model)

colors = ['springgreen', 'olive', 'tomato', 'slateblue']
reported_pops = {}
for variant, color in zip(variants, colors):
    reported_pops_variant = {}
    for cycle in cycles:
        if variant in variants_in_cycle[cycle]:
            reported_pops_variant[cycle] = Population(cycle + '_reported_' + variant, 0,
                                                      cycles[cycle] + ' cases_' + variant,
                                                      hidden=True, color=rotated_color(irc[cycle], color),
                                                      show_sim=True,
                                                      report_noise=True, report_noise_par=report_noise_par,
                                                      report_backlog_par=report_backlog_par, report_days=report_days)
            bc_model.add_connector(
                Propagator(cycle + '_testing_' + variant, symptomatic_pops[variant][cycle],
                           reported_pops_variant[cycle],
                           reported_fraction, reported_delay))
    reported_pops[variant] = reported_pops_variant

    reported_pop_variant = Population('reported_' + variant, 0,
                                      variants[variant] + ' cases',
                                      hidden=True, color=color, show_sim=True,
                                      report_noise=True, report_noise_par=report_noise_par,
                                      report_backlog_par=report_backlog_par, report_days=report_days)

    # do this collection now so it is in sync with all_reported
    populations = []
    for key in reported_pops_variant:
        populations.append(reported_pops_variant[key])

    bc_model.add_connector(
        Collector(reported_pop_variant.name, populations, reported_pop_variant))

# setup collectors for reported by cycle and all reports

color = 'forestgreen'
all_reported_pops = {}
for cycle in cycles:
    reported_pops_cycle = {}
    for variant in variants:
        if variant in variants_in_cycle[cycle]:
            reported_pops_cycle[variant] = reported_pops[variant][cycle]
            all_reported_pops[cycle + variant] = reported_pops[variant][cycle]

    reported_pop_cycle = Population(cycle + '_reported', 0,
                                    cycles[cycle] + ' cases',
                                    hidden=True, color=rotated_color(irc[cycle], color), show_sim=True,
                                    report_noise=True, report_noise_par=report_noise_par,
                                    report_backlog_par=report_backlog_par, report_days=report_days)

    # do this collection now so it is in sync with all_reported
    populations = []
    for key in reported_pops_cycle:
        populations.append(reported_pops_cycle[key])

    bc_model.add_connector(
        Collector(reported_pop_cycle.name, populations, reported_pop_cycle))

reported_pop = Population('reported', 0,
                          'Positives and reporting anomalies',
                          hidden=False, color='forestgreen', show_sim=True,
                          report_noise=True, report_noise_par=report_noise_par,
                          report_backlog_par=report_backlog_par, report_days=report_days)

# collector_populations[reported_pop] = all_reported_pops
# The reported population is subsequently used to adjust the vaccine candidates
# So, the collection has to be done now, rather than at the end

populations = []
for key in all_reported_pops:
    populations.append(all_reported_pops[key])

bc_model.add_connector(
    Collector(reported_pop.name, populations, reported_pop))

# include reporting anomalies (in overall reports only)
# ---------------------------

report_anomalies_pop = Population('report anomalies', 0,
                                  'Represents anomalous batches of reports',
                                  color='lightcyan')

anomaly_fraction = Parameter('anomaly_frac', 1., 0., 1.,
                             'fraction of anomalous reports included == 1')

anomaly_delay_pars = {
    'mean': Parameter('anomaly_delay_mean', 7., 0., 50.,
                      'mean time from injection of anomalies to being reported'),
    'sigma': Parameter('anomaly_delay_sigma', 3., 0.01, 20.,
                       'standard deviation of times from injection of anomalies to being reported')
}
anomaly_delay = Delay('anomaly_delay', 'gamma', anomaly_delay_pars, bc_model)

bc_model.add_connector(
    Propagator('report anomalies to reported', report_anomalies_pop,
               reported_pop, anomaly_fraction, anomaly_delay))

#############################
# SYMPTOMS -> HOSPITALIZATION
#############################
#
# To ensure independence, keep track of the different types of hospitalization
# as separate streams. Since the total probability for entering hospital is low,
# the independent treatment is approximately correct
#
# - Severity may depend on vaccination and previous infection status
# hospitalization fractions are specified separately for each cycle
# - Severity may depend on strain
# - To make this backward compatible, include new parameters that
# scale the nominal fractions (nominal: 'os' and 'o')
#
#
# 2 independent hospitalized streams:
# -> non-ICU hospitalized
# -> ICU

# population to track total number in hospital (non_icu + icu admissions)

hosp_noise_par = Parameter('hosp_noise', 0.1, 0., 1.,
                           'Hospital report noise parameter')
hosp_backlog_par = Parameter('hosp_backlog', 0.5, 0., 1.,
                             'Hospital report backlog parameter')
hosp_report_days = Parameter('hosp_report_days', 127, 0, 127, 'days of week with reporting (bit encoded)',
                             parameter_type='int')

hospitalized_pop = Population('hospitalized', 0,
                              'Total hospitalization cases', color='slategrey', show_sim=True,
                              report_noise=True, report_noise_par=hosp_noise_par,
                              report_backlog_par=hosp_backlog_par, report_days=hosp_report_days)

icu_pop = Population('icu admissions', 0, 'Total icu admissions', color='deeppink', show_sim=False)

# Non-ICU
# -------

non_icu_hospitalized_delay_pars = {
    'mean': Parameter('non_icu_hosp_delay_mean', 6., 0., 50.,
                      'mean time from symptoms to non_icu hospitalization', hidden=False),
    'sigma': Parameter('non_icu_hosp_delay_sigma', 3., 0.01, 20.,
                       'standard deviation of times from symptoms to non_icu hospitalization')
}
non_icu_hospitalized_delay = Delay('non_icu_hosp_delay', 'gamma', non_icu_hospitalized_delay_pars, bc_model)

non_icu_released_pop = Population('non_icu_rel', 0,
                                  'Hospitalized not needing ICU and released',
                                  hidden=True, color='mistyrose', show_sim=False)

non_icu_delay_pars = {
    'mean': Parameter('non_icu_rel_delay_mean', 10., 0., 50.,
                      'mean time from non_icu hospital to release', hidden=False),
    'sigma': Parameter('non_icu_rel_delay_sigma', 3., 0.01, 20.,
                       'standard deviation of times from non_icu hospital to release')
}
non_icu_delay = Delay('non_icu_rel_delay', 'gamma', non_icu_delay_pars, bc_model)

non_icu_delay_pars_x = {
    'mean': Parameter('non_icu_rel_delay_mean_x', 10., 0., 50.,
                      'mean time from non_icu hospital to release: xariant', hidden=False),
    'sigma': Parameter('non_icu_rel_delay_sigma_x', 3., 0.01, 20.,
                       'standard deviation of times from non_icu hospital to release: xariant')
}
non_icu_delay_x = Delay('non_icu_rel_delay_x', 'gamma', non_icu_delay_pars_x, bc_model)

release_fraction = Parameter('release_frac', 1., 1., 1.,
                             'fraction for all released == 1', hidden=True)

# nominal hospitalized fraction ('os' and 'o)
nih_fraction = Parameter('non_icu_hosp_frac', 0.2, 0., 1.,
                         'nominal fraction of those with symptoms who will be admitted to non_icu hospital',
                         hidden=False)

nih_scale_cycle = {}
for cycle in cycles:
    if cycle != 'os':
        nih_scale_cycle[cycle] = Parameter(cycle + '_nih_scale', 1., 0., 1.,
                                           cycles[cycle] + ' scale factor for fraction of symptomatic who ' +
                                           'are admitted to non_icu hospital', hidden=True)

nih_scale_variant = {}
for variant in variants:
    if variant != 'o':
        nih_scale_variant[variant] = Parameter('nih_scale_' + variant, 1., 0., 1.,
                                               variants[variant] + ' scale factor for fraction of symptomatic who ' +
                                               'are admitted to non_icu hospital', hidden=True)

colors = ['indianred', 'hotpink', 'rosybrown', 'brown']
non_icu_hospitalized_pops = {}
for variant, color in zip(variants, colors):
    for cycle in cycles:
        if variant in variants_in_cycle[cycle]:
            non_icu_hospitalized_pops[cycle + variant] = Population(cycle + '_non_icu_hospitalized_' + variant, 0,
                                                                    cycles[cycle] + ' non_icu hospitalization ' +
                                                                    variants[variant] +
                                                                    ' cases', color=rotated_color(irc[cycle], color),
                                                                    show_sim=False)

            # scale the nih fraction as needed
            if variant == 'o':
                if cycle == 'os':
                    nih_frac = nih_fraction
                else:
                    nih_frac = Operator([nih_fraction, nih_scale_cycle[cycle]], '*')
            else:
                if cycle == 'os':
                    nih_frac = Operator([nih_fraction, nih_scale_variant[variant]], '*')
                else:
                    nih_frac = Operator([nih_fraction, nih_scale_cycle[cycle],
                                         nih_scale_variant[variant]], '**')

            bc_model.add_connector(
                Propagator(cycle + ' symptomatic to non_icu hospital: ' + variant, symptomatic_pops[variant][cycle],
                           non_icu_hospitalized_pops[cycle + variant], nih_frac, non_icu_hospitalized_delay))

            bc_model.add_connector(
                Adder('include non_icu in hospitalized:' + cycle + '_' + variant,
                      non_icu_hospitalized_pops[cycle + variant], hospitalized_pop))

            # Those with non_icu hospitalization get released eventually
            # Data from SA suggests this is shorter for omicron
            if variant !='x':
                bc_model.add_connector(
                    Propagator(cycle + ' non_icu hospital to released ' + variant,
                               non_icu_hospitalized_pops[cycle + variant],
                               non_icu_released_pop, release_fraction, non_icu_delay))
            else:
                bc_model.add_connector(
                    Propagator(cycle + ' non_icu hospital to released: ' + variant,
                               non_icu_hospitalized_pops[cycle + variant],
                               non_icu_released_pop, release_fraction, non_icu_delay_x))

non_icu_hospitalized_pop = Population('non_icu_hospitalized', 0,
                                      'Total non_icu hospitalization cases', color='dimgrey', show_sim=False)
collector_populations[non_icu_hospitalized_pop] = non_icu_hospitalized_pops

# ICU
# ---

to_icu_delay_pars = {
    'mean': Parameter('to_icu_delay_mean', 4., 0., 50.,
                      'mean time from symptoms to icu', hidden=False),
    'sigma': Parameter('to_icu_delay_sigma', 2., 0.01, 20.,
                       'standard deviation of times from symptoms to icu')
}
to_icu_delay = Delay('to_icu_delay', 'gamma', to_icu_delay_pars, bc_model)

# for backwards compatability, these parameters refer to icu as non_vent_icu

non_vent_icu_delay_pars = {
    'mean': Parameter('non_vent_icu_delay_mean', 14., 0., 50.,
                      'mean time from non-vent icu admission to release', hidden=False),
    'sigma': Parameter('non_vent_icu_delay_sigma', 5., 0.01, 20.,
                       'standard deviation of times from non-vent icu admission to release')
}
non_vent_icu_delay = Delay('in_icu_delay', 'gamma', non_vent_icu_delay_pars, bc_model)

non_vent_icu_delay_pars_x = {
    'mean': Parameter('non_vent_icu_delay_mean_x', 14., 0., 50.,
                      'mean time from non-vent icu admission to release: xariant', hidden=False),
    'sigma': Parameter('non_vent_icu_delay_sigma_x', 5., 0.01, 20.,
                       'standard deviation of times from non-vent icu admission to release: xariant')
}
non_vent_icu_delay_x = Delay('in_icu_delay_x', 'gamma', non_vent_icu_delay_pars_x, bc_model)

# nominal icu fraction ('os' and 'o)
icu_fraction = Parameter('icu_frac', 0., 0., 1.,
                         'fraction of symptomatic people who go to ' +
                         'icu', hidden=False)

icu_scale_cycle = {}
for cycle in cycles:
    if cycle != 'os':
        icu_scale_cycle[cycle] = Parameter(cycle + '_icu_scale', 1., 0., 1.,
                                           cycles[cycle] + ' scale factor for fraction of symptomatic who ' +
                                           'are admitted to icu', hidden=True)

icu_scale_variant = {}
for variant in variants:
    if variant != 'o':
        icu_scale_variant[variant] = Parameter('icu_scale_' + variant, 1., 0., 1.,
                                               variants[variant] + ' scale factor for fraction of symptomatic who ' +
                                               'are admitted to icu', hidden=True)

icu_released_pop = Population('icu_rel', 0,
                              'admitted to ICU and released:',
                              hidden=True, color='lightgrey', show_sim=False)

colors = ['deeppink', 'hotpink', 'rosybrown', 'brown']
icu_pops = {}
for variant, color in zip(variants, colors):
    for cycle in cycles:
        if variant in variants_in_cycle[cycle]:
            icu_pops[cycle + variant] = Population(cycle + ' icu admissions: ' + variant, 0, 'People admitted to ICU',
                                                   hidden=True, color=rotated_color(irc[cycle], color), show_sim=True)

            # scale the icu fraction as needed
            if variant == 'o':
                if cycle == 'os':
                    icu_frac = icu_fraction
                else:
                    icu_frac = Operator([icu_fraction, icu_scale_cycle[cycle]], '*')
            else:
                if cycle == 'os':
                    icu_frac = Operator([icu_fraction, icu_scale_variant[variant]], '*')
                else:
                    icu_frac = Operator([icu_fraction, icu_scale_cycle[cycle],
                                         icu_scale_variant[variant]], '**')

            bc_model.add_connector(
                Propagator(cycle + ' symptomatic to icu ' + variant, symptomatic_pops[variant][cycle],
                           icu_pops[cycle + variant], icu_frac, to_icu_delay))

            bc_model.add_connector(
                Adder('include icu in hospitalized ' + cycle + '_' + variant, icu_pops[cycle + variant],
                      hospitalized_pop))
            bc_model.add_connector(
                Adder('include this in total icu ' + cycle + '_' + variant, icu_pops[cycle + variant], icu_pop))

            # Those admitted to icu get released eventually
            # data from SA suggest reduced times
            if variant != 'x':
                bc_model.add_connector(
                    Propagator(cycle + ' icu to released ' + variant,
                               icu_pops[cycle + variant],
                               icu_released_pop, release_fraction, non_vent_icu_delay))
            else:
                bc_model.add_connector(
                    Propagator(cycle + ' icu to released ' + variant,
                               icu_pops[cycle + variant],
                               icu_released_pop, release_fraction, non_vent_icu_delay_x))

# make a copy of hospital admissions to keep track of how many remain in hospital
# oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo

in_hospital_pop = Population('in_hospital', 0,
                             'People currently in hospital',
                             hidden=False, color='darkcyan', show_sim=True)
bc_model.add_connector(
    Adder('copy hospitalizations', hospitalized_pop, in_hospital_pop))

# make a copy of icu admissions to keep track of how many remain in icu
# ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo

in_icu_pop = Population('in_icu', 0,
                        'People currently in ICU', hidden=False, color='hotpink', show_sim=True)

bc_model.add_connector(
    Adder('copy icu admissions', icu_pop, in_icu_pop))

# do subtractions to make the in_hospital etc correct
# oooooooooooooooooooooooooooooooooooooooooooooooooooo

bc_model.add_connector(
    Subtractor('remove non-ICU released', in_hospital_pop, non_icu_released_pop))

bc_model.add_connector(
    Subtractor('remove non-vent released from icu', in_icu_pop, icu_released_pop))
bc_model.add_connector(
    Subtractor('remove non-vent released from hospital', in_hospital_pop, icu_released_pop))

# adjust other populations as required
# ooooooooooooooooooooooooooooooooooooo

bc_model.add_connector(
    Subtractor('subtract deaths from total', total_pop, deaths_pop))

# some with a positive test will decide not to get vaccinated
frac_report_novacc_par = Parameter('frac_report_novacc', 1.0, 0., 1.,
                                   'Fraction of those reported who decide against vaccination')

# vaccination and booster accounting

bc_model.add_connector(
    Subtractor('subtract vaccinated from vacc cand', vaccan_pop, vaccinated_pop))

if not no_bt:
    bc_model.add_connector(
        Subtractor('subtract boosted from booster cand', boostcan_pop, boosted_pop))

bc_model.add_connector(
    Subtractor('subtract usefully vaccinated from the sus vacc cand', susvaccan_pop, usefully_vaccinated_pop))

for variant in variants:
    for cycle in cycles:
        if variant in variants_in_cycle[cycle]:
            bc_model.add_connector(
                Subtractor('subtract some of infected from sus vacc cand: ' + cycle + '_' + variant, susvaccan_pop,
                           infected_pops[variant][cycle], ratio_populations=[susvaccan_pop, susceptible_pop]))

rep_frac_pop = Population('frac', 0.)

bc_model.add_connector(
    Adder('calc rep frac', reported_pop, rep_frac_pop, ratio_populations=[vaccinated_pop, total_pop]))

bc_model.add_connector(
    Adder('add vaccinated reported', rep_frac_pop, vaccan_pop))

bc_model.add_connector(
    Subtractor('subtract immunized reported', vaccan_pop, rep_frac_pop, scale_factor=vaccine_effectiveness))

bc_model.add_connector(
    Subtractor('subtract reported', vaccan_pop, reported_pop, scale_factor=frac_report_novacc_par))

# Do all collections at the end

for collector_population in collector_populations:
    populations = []
    for key in collector_populations[collector_population]:
        populations.append(collector_populations[collector_population][key])

    bc_model.add_connector(
        Collector(collector_population.name, populations, collector_population)
    )

# transitional aspects of the model
# oooooooooooooooooooooooooooooooooo

n_rate_transitions = {'o': 14, 'v': 6, 'w': 12, 'x': 4}
n_rt_visible = {'o': 4, 'v': 0, 'w': 0, 'x': 0}

for variant in variants:
    suffix = ''
    if variant != 'o':
        suffix = '_' + variant

    times = []
    alphas = [Parameter('alpha_0' + suffix, 0.4, 0., 2.,
                        'initial transmission rate', hidden=False)]
    for i in range(n_rate_transitions[variant]):
        j = i + 1
        hidden = j > n_rt_visible[variant]
        enabled = i == 0
        times.append(Parameter('trans_rate_' + str(j) + '_time' + suffix, 20 + i, 0, 600,
                               'day of trans rate change ' + 'str(j)',
                               parameter_type='int', hidden=hidden))
        alphas.append(Parameter('alpha_' + str(j) + suffix, 0.1, 0., 2.,
                                'alpha' + suffix + ' after transition ' + str(j), hidden=hidden))
        bc_model.add_transition(
            Modifier('trans_rate_' + str(j) + suffix, 'rel_days', times[i], trans_rates[variant],
                     alphas[i], alphas[j], enabled=enabled, model=bc_model))

# vaccination start

n_vacc_periods = 14

for i in range(n_vacc_periods):
    j = i + 1
    vaccination_time = Parameter('vacc_time_' + str(j), 75 + i, 0, 600,
                                 'first day of vaccination period ' + str(j),
                                 parameter_type='int')
    vaccination_number = Parameter('vacc_number_' + str(j), 10., 0., 5000000.,
                                   'change in number vaccinated each day for period ' + str(j))
    bc_model.add_transition(
        Injector('vaccination_' + str(j), 'rel_days', vaccination_time, daily_vaccinated_pop,
                 vaccination_number, enabled=False, model=bc_model))

# booster scheduls

if not no_bt:
    n_boost_periods = 4

    for i in range(n_boost_periods):
        j = i + 1
        boost_time = Parameter('boost_time_' + str(j), 75 + i, 0, 600,
                               'first day of booster period ' + str(j),
                               parameter_type='int')
        boost_number = Parameter('boost_number_' + str(j), 10., 0., 5000000.,
                                 'change in number boosted each day for period ' + str(j))
        bc_model.add_transition(
            Injector('booster_' + str(j), 'rel_days', boost_time, daily_boosted_pop,
                     boost_number, enabled=False, model=bc_model))

# outbreaks (in the 'os' susceptible population)

outbreak_fraction = Parameter('outbreak_frac', 1., 0., 1.,
                              'fraction of infected in outbreak active ==1')

n_outbreaks = {'o': 7, 'v': 3, 'w': 4, 'x': 3}

for variant in variants:
    suffix = ''
    if variant != 'o':
        suffix = '_' + variant

    outbreak_delay_pars = {
        'mean': Parameter('outbreak' + suffix + '_delay_mean', 7., 0., 50.,
                          'mean delay time for outbreak' + suffix, hidden=False),
        'sigma': Parameter('outbreak' + suffix + '_delay_sigma', 1., 0.01, 20.,
                           'standard deviation of outbreak' + suffix + ' times',
                           hidden=False)
    }
    outbreak_delay = Delay('outbreak' + suffix + '_delay', 'gamma', outbreak_delay_pars, bc_model)

    outbreak_pop = Population('outbreaks' + suffix, 0, variants[variant] + ' infection outbreaks')

    for i in range(n_outbreaks[variant]):
        vid = ''
        if i > 0:
            vid = str(i)
        if variant != 'o':
            vid = variant + vid

        outbreak_time = Parameter('outbreak_' + vid + '_time', 14, 0, 800,
                                  'number of days since t0 when outbreak_' + vid + ' established',
                                  parameter_type='int', hidden=False)

        outbreak_number = Parameter('outbreak_' + vid + '_number', 10., 0., 50000.,
                                    'number of infections in outbreak_' + vid,
                                    hidden=False)

        bc_model.add_transition(
            Injector('outbreak_' + vid, 'rel_days', outbreak_time, outbreak_pop,
                     outbreak_number, enabled=False, model=bc_model))

    bc_model.add_connector(
        Propagator('outbreaks to os infected' + suffix, outbreak_pop, infected_pops[variant]['os'],
                   outbreak_fraction, outbreak_delay))

# reporting anomalies

anomaly_1_time = Parameter('rep_anomaly_1_time', 40, 0, 600,
                           'number of days since t0 when reporting anomaly_1 occurs',
                           parameter_type='int', hidden=False)

anomaly_1_number = Parameter('rep_anomaly_1_number', 50., 0., 50000.,
                             'number of anomalous reports in anomaly_1',
                             hidden=False)

bc_model.add_transition(
    Injector('rep_anomaly_1', 'rel_days', anomaly_1_time, report_anomalies_pop,
             anomaly_1_number, enabled=False, model=bc_model))

anomaly_2_time = Parameter('rep_anomaly_2_time', 80, 0, 600,
                           'number of days since t0 when reporting anomaly_2 occurs',
                           parameter_type='int')

anomaly_2_number = Parameter('rep_anomaly_2_number', 50., 0., 50000.,
                             'number of anomalous reports in anomaly_2')

bc_model.add_transition(
    Injector('rep_anomaly_2', 'rel_days', anomaly_2_time, report_anomalies_pop,
             anomaly_2_number, enabled=False, model=bc_model))

anomaly_3_time = Parameter('rep_anomaly_3_time', 80, 0, 600,
                           'number of days since t0 when reporting anomaly_3 occurs',
                           parameter_type='int')

anomaly_3_number = Parameter('rep_anomaly_3_number', 50., 0., 50000.,
                             'number of anomalous reports in anomaly_3')

bc_model.add_transition(
    Injector('rep_anomaly_3', 'rel_days', anomaly_3_time, report_anomalies_pop,
             anomaly_3_number, enabled=False, model=bc_model))

# Linear modification to testing fraction:

reported_fractions = [Parameter('report_frac_0', 0.95, 0., 1.,
                                'initial reported fraction', hidden=False)]

for i in range(4):
    j = i + 1
    hidden = j != 1
    reported_fractions.append(Parameter('report_frac_' + str(j), 0.96, 0., 1.,
                                        'reported fraction after mod ' + str(j), hidden=hidden))

    mod_reported_fraction_time = Parameter('report_frac_time_' + str(j), 50, 0, 800,
                                           'time of mod reported frac ' + str(j),
                                           parameter_type='int', hidden=hidden)

    mod_reported_fraction_nstep = Parameter('report_frac_nstep_' + str(j), -1, -1, 800,
                                            'number of days to apply report_frac ' + str(j) + ' linear modification',
                                            parameter_type='int', hidden=hidden)

    bc_model.add_transition(
        Modifier('mod_report_frac_' + str(j), 'rel_days', mod_reported_fraction_time, reported_fraction,
                 reported_fractions[i], reported_fractions[j], enabled=False, model=bc_model,
                 linear=True, n_step=mod_reported_fraction_nstep))

# define boot parameters
# ooooooooooooooooooooooo

bc_model.boot_setup(original_contagious_pop, 0.1,
                    exclusion_populations=[total_pop, susceptible_pop, susceptible_pops['os'],
                                           susvaccan_pop, vaccan_pop])

bc_model.save_file('ref_model_' + str(version) + '_' + str(subversion) + '.pypm')
