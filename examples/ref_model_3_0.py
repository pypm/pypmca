# -*- coding: utf-8 -*-
"""
pyPM.ca reference model #3.0

To make the critical parameters less dependent of each other, break the measured
populations into separate paths.

Reduced the lag time between contagious and reported, (as compared to reference model #1)
to better match the situation in BC and elsewhere - transitions now occur where they should

Add a contact tracing population, to allow another mechanism to remove people from
the contagious population early. Initial value has p=0, so it does not act until it
is switched on.

Version 2.1:
 - turn on negative binomial for multiplier. Add parameter to control dispersion (set to 0.999: tiny effect, at p = 0.5, variance is double poisson)
 - turn on reporting errors for reported. Add parameter to control level (set to 1: no effect)
 - add capability to include reporting anomalies: new pop = tested positive receives contact and symp, reduces contagious. added to reported
   new pop = reporting anomalies: propagates to reported. Add new transition: injector to reporting anomalies

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
 - add third infection cycle for B.1.617.2 variant (use the letter w in place of v for naminng)
 - add second outbreaks for the two variants
 - include waning immunity (separately for vaccination and natural):
   - add natural immunized population
   - propagate (with long time delays) from immunized and naturally immunized to "lost immunity" that
     then increases susceptible

Version 3.0:
 - reported does not necessarily lead to removal from vaccine candidates
   (new parameter, frac_report_novacc, has default=1 to be backward compatible)
 - more wariant transitions (delta)
 - 4th strain, in case omicron is real: w-> x
 - allow breakthrough populations to be tracked (reported, hospitalized, deaths)
   - breakthrough is any person vaccinated who is not immunized (waning parameters determine susceptible size)
   - share the alphas of the naive populations (adjust waning parameters to incorporate altered immunity)
   - the original susceptible includes all. bt_susceptible produces bt_infected, which propogate to both bt_contagious
     and contagious
 - include boosters

@author: karlen
"""

from pypmca import Model, Population, Delay, Parameter, Multiplier, Propagator, \
    Splitter, Adder, Subtractor, Chain, Modifier, Injector

# This is not listed in package requirements:
# Users generally do not need to run this code.
import matplotlib


def complementary_color(c):
    my_hex = matplotlib.colors.to_hex(c)
    rgb = (my_hex[1:3], my_hex[3:5], my_hex[5:7])
    comp = ['%02X' % (255 - int(a, 16)) for a in rgb]
    return '#' + ''.join(comp)


# Reference model for BC

bc_model = Model('ref_model_3_0')
bc_model.set_t0(2020, 3, 1)

# Initialization

initial_pop_par = Parameter('N_0', 5000000, 5000, 50000000,
                            'Population of the region at t0', 'int')

total_pop = Population('total', initial_pop_par,
                       'total population of the region', color='black')
susceptible_pop = Population('susceptible', initial_pop_par,
                             'number of people who could become infected', color='cornflowerblue')
infected_pop = Population('infected', 0,
                          'total number of people ever infected', color='orange')

# Define the infection cycle
# oooooooooooooooooooooooooooo

initial_contagious_par = Parameter('cont_0', 10., 0., 5000.,
                                   'Number of contagious people at t0',
                                   hidden=False)

contagious_pop = Population('contagious', initial_contagious_par,
                            'number of people that can cause someone to become infected',
                            hidden=False, color='red')

# this value is only used if the transition is removed
trans_rate = Parameter('alpha', 0.4, 0., 2.,
                       'mean number of people that a contagious person infects ' +
                       'per day', hidden=True)
fast_delay = Delay('fast', 'fast', model=bc_model)

neg_binom_par = Parameter('neg_binom_p', 0.5, 0.001, 0.999,
                          'Dispersion parameter p for neg binom')

bc_model.add_connector(
    Multiplier('infection cycle', [susceptible_pop, contagious_pop, total_pop],
               infected_pop, trans_rate, fast_delay, bc_model,
               distribution='nbinom', nbinom_par=neg_binom_par))

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

bc_model.add_connector(
    Propagator('infected to contagious', infected_pop,
               contagious_pop, contagious_frac, contagious_delay))

# Include a second infection cycle for another variant (eg. B.1.1.7)
# ------------------------------------------------------------------

infected_pop_v = Population('infected_v', 0,
                            'total number of people ever infected with variant',
                            hidden=True, color='darkgoldenrod')

contagious_pop_v = Population('contagious_v', 0,
                              'number of people that can cause variant infections',
                              hidden=True, color='rosybrown')

# this value is only used if the transition is removed:
trans_rate_v = Parameter('alpha_v', 0.4, 0., 2.,
                         'mean number of people that a contagious person infects ' +
                         'per day', hidden=True)

bc_model.add_connector(
    Multiplier('infection cycle_v', [susceptible_pop, contagious_pop_v, total_pop],
               infected_pop_v, trans_rate_v, fast_delay, bc_model,
               distribution='nbinom', nbinom_par=neg_binom_par))

bc_model.add_connector(
    Propagator('infected to contagious_v', infected_pop_v,
               contagious_pop_v, contagious_frac, contagious_delay))

# end of second infection cycle

# Include a third infection cycle for another variant (eg. B.1.617.2) - refer as wariant
# ------------------------------------------------------------------

infected_pop_w = Population('infected_w', 0,
                            'total number of people ever infected with wariant',
                            hidden=True, color='sienna')

contagious_pop_w = Population('contagious_w', 0,
                              'number of people that can cause wariant infections',
                              hidden=True, color='maroon')

# this value is only used if the transition is removed:
trans_rate_w = Parameter('alpha_w', 0.4, 0., 2.,
                         'mean number of people that a contagious person infects ' +
                         'per day', hidden=True)

bc_model.add_connector(
    Multiplier('infection cycle_w', [susceptible_pop, contagious_pop_w, total_pop],
               infected_pop_w, trans_rate_w, fast_delay, bc_model,
               distribution='nbinom', nbinom_par=neg_binom_par))

bc_model.add_connector(
    Propagator('infected to contagious_w', infected_pop_w,
               contagious_pop_w, contagious_frac, contagious_delay))

# end of third infection cycle

# Include a fourth infection cycle for another variant (eg. omicron) - refer as xariant
# ------------------------------------------------------------------

infected_pop_x = Population('infected_x', 0,
                            'total number of people ever infected with xariant',
                            hidden=True, color='sandybrown')

contagious_pop_x = Population('contagious_x', 0,
                              'number of people that can cause xariant infections',
                              hidden=True, color='lightcoral')

# this value is only used if the transition is removed:
trans_rate_x = Parameter('alpha_x', 0.4, 0., 2.,
                         'mean number of people that a contagious person infects ' +
                         'per day', hidden=True)

bc_model.add_connector(
    Multiplier('infection cycle_x', [susceptible_pop, contagious_pop_x, total_pop],
               infected_pop_x, trans_rate_x, fast_delay, bc_model,
               distribution='nbinom', nbinom_par=neg_binom_par))

bc_model.add_connector(
    Propagator('infected to contagious_x', infected_pop_x,
               contagious_pop_x, contagious_frac, contagious_delay))

# end of fourth infection cycle

# ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
# Define the infection cycles for those vaccinated but no longer immunized (breakthrough)
# ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo

bt_susceptible_pop = Population('bt_susceptible', 0,
                                'number of people who could become breakthrough infected',
                                color=complementary_color('cornflowerblue'))
bt_infected_pop = Population('bt_infected', 0,
                             'total number of people ever infected as a breakthough',
                             color=complementary_color('orange'))

# This is a sub population of the overall contagious population:
bt_contagious_pop = Population('bt_contagious', 0,
                               'number of people that became contagious as a breakthrough',
                               hidden=True, color=complementary_color('red'))

# The overall contagious_pop is involved in producing breakthrough infections:
bc_model.add_connector(
    Multiplier('bt_infection cycle', [bt_susceptible_pop, contagious_pop, total_pop],
               bt_infected_pop, trans_rate, fast_delay, bc_model,
               distribution='nbinom', nbinom_par=neg_binom_par))

bc_model.add_connector(
    Propagator('bt_infected to bt_contagious', bt_infected_pop,
               bt_contagious_pop, contagious_frac, contagious_delay))

bc_model.add_connector(
    Adder('include bt_contagious in contagious', bt_contagious_pop, contagious_pop))

# Include a second infection cycle for another variant (eg. B.1.1.7)
# ------------------------------------------------------------------

bt_infected_pop_v = Population('bt_infected_v', 0,
                               'total number of people ever infected with variant as a breakthrough',
                               hidden=True, color=complementary_color('darkgoldenrod'))

bt_contagious_pop_v = Population('bt_contagious_v', 0,
                                 'number of people that can cause variant infections as a breakthrough',
                                 hidden=True, color=complementary_color('rosybrown'))

bc_model.add_connector(
    Multiplier('bt_infection cycle_v', [bt_susceptible_pop, contagious_pop_v, total_pop],
               bt_infected_pop_v, trans_rate_v, fast_delay, bc_model,
               distribution='nbinom', nbinom_par=neg_binom_par))

bc_model.add_connector(
    Propagator('bt_infected_v to bt_contagious_v', bt_infected_pop_v,
               bt_contagious_pop_v, contagious_frac, contagious_delay))

bc_model.add_connector(
    Adder('include bt_contagious_v in contagious_v', bt_contagious_pop_v, contagious_pop_v))

# end of second infection cycle

# Include a third infection cycle for another variant (eg. B.1.617.2) - refer as wariant
# ------------------------------------------------------------------

bt_infected_pop_w = Population('bt_infected_w', 0,
                               'total number of people ever infected with wariant as a breakthrough',
                               hidden=True, color=complementary_color('sienna'))

bt_contagious_pop_w = Population('bt_contagious_w', 0,
                                 'number of people that can cause wariant infections as a breakthrough',
                                 hidden=True, color=complementary_color('maroon'))

bc_model.add_connector(
    Multiplier('bt_infection cycle_w', [bt_susceptible_pop, contagious_pop_w, total_pop],
               bt_infected_pop_w, trans_rate_w, fast_delay, bc_model,
               distribution='nbinom', nbinom_par=neg_binom_par))

bc_model.add_connector(
    Propagator('bt_infected_w to bt_contagious_w', bt_infected_pop_w,
               bt_contagious_pop_w, contagious_frac, contagious_delay))

bc_model.add_connector(
    Adder('include bt_contagious_w in contagious_w', bt_contagious_pop_w, contagious_pop_w))

# end of third infection cycle

# Include a fourth infection cycle for another variant (eg. omicron) - refer as xariant
# ------------------------------------------------------------------

bt_infected_pop_x = Population('bt_infected_x', 0,
                               'total number of people ever infected with xariant as a breakthrough',
                               hidden=True, color=complementary_color('sandybrown'))

bt_contagious_pop_x = Population('bt_contagious_x', 0,
                                 'number of people that can cause xariant infections as a breakthrough',
                                 hidden=True, color=complementary_color('lightcoral'))

bc_model.add_connector(
    Multiplier('bt_infection cycle_x', [bt_susceptible_pop, contagious_pop_x, total_pop],
               bt_infected_pop_x, trans_rate_x, fast_delay, bc_model,
               distribution='nbinom', nbinom_par=neg_binom_par))

bc_model.add_connector(
    Propagator('bt_infected_x to bt_contagious_x', bt_infected_pop_x,
               bt_contagious_pop_x, contagious_frac, contagious_delay))

bc_model.add_connector(
    Adder('include bt_contagious_x in contagious_x', bt_contagious_pop_x, contagious_pop_x))

# end of fourth infection cycle


# keep track of naturally immunized population

nat_immunized_pop = Population('nat immunized', 0,
                               'people who gained immunity through infection',
                               hidden=True, color='gold')

# all forms of infection lead to natural immunity (that eventually wanes)

bc_model.add_connector(
    Adder('natural immunity', infected_pop, nat_immunized_pop))

bc_model.add_connector(
    Adder('natural immunity v', infected_pop_v, nat_immunized_pop))

bc_model.add_connector(
    Adder('natural immunity w', infected_pop_w, nat_immunized_pop))

bc_model.add_connector(
    Adder('natural immunity x', infected_pop_x, nat_immunized_pop))

# breakthrough infection also leads to natural immunity (that eventually wanes again)

bt_nat_immunized_pop = Population('bt nat immunized', 0,
                                  'people who gained immunity through breakthrough infection',
                                  hidden=True, color=complementary_color('gold'))

bc_model.add_connector(
    Adder('bt natural immunity', bt_infected_pop, bt_nat_immunized_pop))

bc_model.add_connector(
    Adder('bt natural immunity v', bt_infected_pop_v, bt_nat_immunized_pop))

bc_model.add_connector(
    Adder('bt natural immunity w', bt_infected_pop_w, bt_nat_immunized_pop))

bc_model.add_connector(
    Adder('bt natural immunity x', bt_infected_pop_x, bt_nat_immunized_pop))

# waning immunity

# Specified by a delay distribution F(t) and a fraction of the immunizations that wane, f.
# Waning is not considered to be equivalent to returning to a naive state, but instead
# the immunity is not as effective after waning.
# The product f*F(t) is the cumulative probability of waning immunity, which accounts
# for the combined effect of:
# - increasing probability for an individual to wane
# - the reducing effectiveness of the immunity of waned individuals

waned_nat_immunity_pop = Population('waned nat immunity', 0,
                                    'people who lost nat immunity some time after gaining nat immunity',
                                    hidden=True, color='mediumvioletred')

nat_waned_delay_pars = {
    'mean': Parameter('nat_waned_delay_mean', 365., 1., 10000., 'mean time from natural immunity to waned immunity'),
    'sigma': Parameter('nat_waned_delay_sigma', 50., 1., 500.,
                       'standard deviation of times from natural immunity to waned immunity')
}
nat_waned_delay = Delay('nat_waned_delay', 'gamma', nat_waned_delay_pars, bc_model)

nat_waned_fraction = Parameter('nat_waned_frac', 1., 0., 1.,
                               'fraction of natural immunizations that wane')

bc_model.add_connector(
    Propagator('nat waned immunity', nat_immunized_pop, waned_nat_immunity_pop, nat_waned_fraction, nat_waned_delay))

# Once waned, the susceptible population increases

bc_model.add_connector(
    Adder('waned nat to susceptible', waned_nat_immunity_pop, susceptible_pop))

# And again for the breakthroughs

bt_waned_nat_immunity_pop = Population('bt waned nat immunity', 0,
                                       'people who lost nat immunity some time after gaining bt nat immunity',
                                       hidden=True, color=complementary_color('mediumvioletred'))

bc_model.add_connector(
    Propagator('bt nat waned immunity', bt_nat_immunized_pop, bt_waned_nat_immunity_pop,
               nat_waned_fraction, nat_waned_delay))

# Once waned, the bt_susceptible population increases

bc_model.add_connector(
    Adder('bt waned nat to bt susceptible', bt_waned_nat_immunity_pop, bt_susceptible_pop))

# Include vaccination:
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

# Once waned, the bt susceptible population increases

bc_model.add_connector(
    Adder('waned to susceptible', waned_vac_immunity_pop, bt_susceptible_pop))

#####################
# Include boosters:
#####################
# Grow the population of booster candidates (delay following first dose)

boostcan_pop = Population('boost cand', 0,
                          'booster candidates', color=complementary_color('brown'))

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

boosted_pop = Population('boosted', 0, 'people boosted', color=complementary_color('skyblue'))

bc_model.add_connector(
    Multiplier('boosting', [boostcan_pop, daily_boosted_pop, boostcan_pop],
               boosted_pop, boost_fraction, fast_delay, bc_model))

# Only a fraction of those were susceptible (designate those as usefully boosted)

usefully_boosted_pop = Population('usefully boosted', 0,
                                  'people who were boosted when susceptible', color=complementary_color('navy'))

bc_model.add_connector(
    Adder('boosting bt susceptibles', boosted_pop, usefully_boosted_pop,
          ratio_populations=[bt_susceptible_pop, boostcan_pop]))

# Now re-immunize some of those susceptible people who were boosted

reimmunized_pop = Population('reimmunized', 0, 'number of bt susceptible people who were immunized by booster')
booster_effectiveness = Parameter('booster_eff', 0.8, 0., 1.,
                                  'probability that a bt susceptible person gains immunity when boosted')

bc_model.add_connector(
    Propagator('reimmunization', usefully_boosted_pop, reimmunized_pop, booster_effectiveness, immunized_delay))

# Waning booster immunity

waned_boost_immunity_pop = Population('waned booster immunity', 0,
                                    'people who lost immunity some time after gaining immunity through a booster',
                                    hidden=True, color=complementary_color('mediumpurple'))

bc_model.add_connector(
    Propagator('boost waned immunity', reimmunized_pop, waned_boost_immunity_pop, vac_waned_fraction, vac_waned_delay))

# Once waned, the bt susceptible population increases

bc_model.add_connector(
    Adder('booster waned to susceptible', waned_boost_immunity_pop, bt_susceptible_pop))

# The contagious either recover or die
# This split is only used to track the deaths.
# The removal of recoveries/positive tests etc done in a separate path
# oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo

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

deaths_pop = Population('deaths', 0,
                        'people who have died from the illness', hidden=False,
                        color='indigo', show_sim=True,
                        report_noise=True, report_noise_par=death_noise_par,
                        report_backlog_par=death_backlog_par, report_days=death_report_days,
                        report_noise_weekly=True)

recover_fraction = Parameter('recover_frac', 0.99, 0., 1.,
                             'fraction of infected people who recover', hidden=False)
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

bc_model.add_connector(
    Splitter('recovery', contagious_pop, [recovered_pop, deaths_pop],
             [recover_fraction], [recover_delay, death_delay]))

bc_model.add_connector(
    Splitter('recovery_v', contagious_pop_v, [recovered_pop, deaths_pop],
             [recover_fraction], [recover_delay, death_delay]))

bc_model.add_connector(
    Splitter('recovery_w', contagious_pop_w, [recovered_pop, deaths_pop],
             [recover_fraction], [recover_delay, death_delay]))

bc_model.add_connector(
    Splitter('recovery_x', contagious_pop_x, [recovered_pop, deaths_pop],
             [recover_fraction], [recover_delay, death_delay]))

# Recovery and death following breakthrough infections

bt_recovered_pop = Population('bt recovered', 0, 'People who have recovered from a breakthrough illness ',
                              color=complementary_color('limegreen'))

bt_deaths_pop = Population('bt deaths', 0,
                           'people who have died from a breakthrough illness', hidden=True,
                           color=complementary_color('indigo'), show_sim=True,
                           report_noise=True, report_noise_par=death_noise_par,
                           report_backlog_par=death_backlog_par, report_days=death_report_days,
                           report_noise_weekly=True)

bc_model.add_connector(
    Splitter('bt recovery', bt_contagious_pop, [bt_recovered_pop, bt_deaths_pop],
             [recover_fraction], [recover_delay, death_delay]))

bc_model.add_connector(
    Splitter('bt recovery_v', bt_contagious_pop_v, [bt_recovered_pop, bt_deaths_pop],
             [recover_fraction], [recover_delay, death_delay]))

bc_model.add_connector(
    Splitter('bt recovery_w', bt_contagious_pop_w, [bt_recovered_pop, bt_deaths_pop],
             [recover_fraction], [recover_delay, death_delay]))

bc_model.add_connector(
    Splitter('bt recovery_x', bt_contagious_pop_x, [bt_recovered_pop, bt_deaths_pop],
             [recover_fraction], [recover_delay, death_delay]))

# The newly contagious are split into three groups:
# contact-traced, symptomatic and non-symptomatic.

contact_traced_pop = Population('contact traced', 0,
                                'People identified through contact tracing', color='coral')
contact_traced_fraction = Parameter('contact_traced_frac', 0., 0., 1.,
                                    'fraction of contagious people who are identified first ' +
                                    'through contact tracing', hidden=False)
contact_traced_delay_pars = {
    'mean': Parameter('contact_traced_delay_mean', 2., 0., 50.,
                      'mean time from becoming contagious to being contact traced', hidden=False),
    'sigma': Parameter('contact_traced_delay_sigma', 1., 0.01, 20.,
                       'std dev of times from becoming contagious to being contact traced')
}
contact_traced_delay = Delay('contact_traced_delay', 'gamma', contact_traced_delay_pars,
                             bc_model)

symptomatic_pop = Population('symptomatic', 0,
                             'People who have shown symptoms', color='chocolate')
symptomatic_fraction = Parameter('symptomatic_frac', 0.9, 0., 1.,
                                 'fraction of contagious people who become ' +
                                 'symptomatic', hidden=False)
symptomatic_delay_pars = {
    'mean': Parameter('symptomatic_delay_mean', 2., 0., 50.,
                      'mean time from becoming contagious to having symptoms'),
    'sigma': Parameter('symptomatic_delay_sigma', 1., 0.01, 20.,
                       'std dev of times from becoming contagious to having symptoms')
}
symptomatic_delay = Delay('symptomatic_delay', 'gamma', symptomatic_delay_pars,
                          bc_model)

asymptomatic_recovered_pop = Population('asymptomatic recovered', 0,
                                        'People who have not shown symptoms', color='silver')

asymptomatic_delay_pars = {
    'mean': Parameter('asymp_rec_delay_mean', 12., 0., 50.,
                      'mean time from becoming contagious to recovery (without symptoms)'),
    'sigma': Parameter('asymp_rec_delay_sigma', 4., 0.01, 20.,
                       'std dev of times from becoming contagious to recovery with no symptoms')
}
asymptomatic_delay = Delay('asymp_rec_delay', 'gamma', asymptomatic_delay_pars,
                           bc_model)

bc_model.add_connector(
    Splitter('symptoms', contagious_pop, [contact_traced_pop, symptomatic_pop, asymptomatic_recovered_pop],
             [contact_traced_fraction, symptomatic_fraction],
             [contact_traced_delay, symptomatic_delay, asymptomatic_delay]))

bc_model.add_connector(
    Splitter('symptoms_v', contagious_pop_v, [contact_traced_pop, symptomatic_pop, asymptomatic_recovered_pop],
             [contact_traced_fraction, symptomatic_fraction],
             [contact_traced_delay, symptomatic_delay, asymptomatic_delay]))

bc_model.add_connector(
    Splitter('symptoms_w', contagious_pop_w, [contact_traced_pop, symptomatic_pop, asymptomatic_recovered_pop],
             [contact_traced_fraction, symptomatic_fraction],
             [contact_traced_delay, symptomatic_delay, asymptomatic_delay]))

bc_model.add_connector(
    Splitter('symptoms_x', contagious_pop_x, [contact_traced_pop, symptomatic_pop, asymptomatic_recovered_pop],
             [contact_traced_fraction, symptomatic_fraction],
             [contact_traced_delay, symptomatic_delay, asymptomatic_delay]))

# The newly bt contagious are split into two groups: symptomatic and non-symptomatic.
# No contact tracing compartments for bt (until needed)

bt_symptomatic_pop = Population('bt symptomatic', 0,
                                'People with breakthrough infections who have shown symptoms',
                                color=complementary_color('chocolate'))

bt_symptomatic_fraction = Parameter('bt_symptomatic_frac', 0.5, 0., 1.,
                                    'fraction of breakthrough contagious people who become ' +
                                    'symptomatic', hidden=False)

bt_asymptomatic_recovered_pop = Population('bt asymptomatic recovered', 0,
                                           'People with breakthrough infections who have not shown symptoms',
                                           color=complementary_color('silver'))

bc_model.add_connector(
    Splitter('bt symptoms', bt_contagious_pop,
             [bt_symptomatic_pop, bt_asymptomatic_recovered_pop],
             [symptomatic_fraction],
             [symptomatic_delay, asymptomatic_delay]))

bc_model.add_connector(
    Splitter('bt symptoms_v', bt_contagious_pop_v,
             [bt_symptomatic_pop, bt_asymptomatic_recovered_pop],
             [symptomatic_fraction],
             [symptomatic_delay, asymptomatic_delay]))

bc_model.add_connector(
    Splitter('bt symptoms_w', bt_contagious_pop_w,
             [bt_symptomatic_pop, bt_asymptomatic_recovered_pop],
             [symptomatic_fraction],
             [symptomatic_delay, asymptomatic_delay]))

bc_model.add_connector(
    Splitter('bt symptoms_x', bt_contagious_pop_x,
             [bt_symptomatic_pop, bt_asymptomatic_recovered_pop],
             [symptomatic_fraction],
             [symptomatic_delay, asymptomatic_delay]))

# reporting parameters

report_noise_par = Parameter('report_noise', 1., 0., 1., 'Report noise parameter')

report_backlog_par = Parameter('report_backlog', 1., 0., 1., 'Report backlog parameter')

report_days = Parameter('report_days', 127, 0, 127, 'days of week with reporting (bit encoded)',
                        parameter_type='int')

# additional reporting line from contagious_v to reported_v
# This allows breakdown of different strains... since it is direct (not going through symptoms, need to
# match the fraction and delay parameters with the combination of those for overall reporting

reported_pop_v = Population('reported_v', 0,
                            'variant cases reported',
                            hidden=True, color='olive', show_sim=True,
                            report_noise=True, report_noise_par=report_noise_par,
                            report_backlog_par=report_backlog_par, report_days=report_days)

reported_fraction_v = Parameter('reported_frac_v', 0.8, 0., 1.,
                                'fraction of variant contagious who will ' + \
                                'receive a variant case report')
reported_delay_pars_v = {
    'mean': Parameter('reported_delay_mean_v', 5., 0., 50.,
                      'mean time from becoming contagious to getting variant case report', hidden=True),
    'sigma': Parameter('reported_delay_sigma_v', 2., 0.01, 20.,
                       'standard deviation of times from contagious to getting variant case report')
}
reported_delay_v = Delay('reported_delay_v', 'gamma', reported_delay_pars_v, bc_model)

bc_model.add_connector(
    Propagator('testing_v', contagious_pop_v, reported_pop_v, reported_fraction_v, reported_delay_v))

# additional reporting line from contagious_w to reported_w

reported_pop_w = Population('reported_w', 0,
                            'wariant cases reported',
                            hidden=True, color='tomato', show_sim=True,
                            report_noise=True, report_noise_par=report_noise_par,
                            report_backlog_par=report_backlog_par, report_days=report_days)

reported_fraction_w = Parameter('reported_frac_w', 0.8, 0., 1.,
                                'fraction of wariant contagious who will ' + \
                                'receive a wariant case report')
reported_delay_pars_w = {
    'mean': Parameter('reported_delay_mean_w', 5., 0., 50.,
                      'mean time from becoming contagious to getting wariant case report', hidden=True),
    'sigma': Parameter('reported_delay_sigma_w', 2., 0.01, 20.,
                       'standard deviation of times from contagious to getting wariant case report')
}
reported_delay_w = Delay('reported_delay_w', 'gamma', reported_delay_pars_w, bc_model)

bc_model.add_connector(
    Propagator('testing_w', contagious_pop_w, reported_pop_w, reported_fraction_w, reported_delay_w))

# additional reporting line from contagious_x to reported_x

reported_pop_x = Population('reported_x', 0,
                            'xariant cases reported',
                            hidden=True, color='tomato', show_sim=True,
                            report_noise=True, report_noise_par=report_noise_par,
                            report_backlog_par=report_backlog_par, report_days=report_days)

reported_fraction_x = Parameter('reported_frac_x', 0.8, 0., 1.,
                                'fraction of xariant contagious who will ' + \
                                'receive a xariant case report')
reported_delay_pars_x = {
    'mean': Parameter('reported_delay_mean_x', 5., 0., 50.,
                      'mean time from becoming contagious to getting xariant case report', hidden=True),
    'sigma': Parameter('reported_delay_sigma_x', 2., 0.01, 20.,
                       'standard deviation of times from contagious to getting xariant case report')
}
reported_delay_x = Delay('reported_delay_x', 'gamma', reported_delay_pars_x, bc_model)

bc_model.add_connector(
    Propagator('testing_x', contagious_pop_x, reported_pop_x, reported_fraction_x, reported_delay_x))

# The removal of people from the contagious population is treated separately from
# reporting, to allow the reported number to be independent of delay in removal

removed_pop = Population('removed', 0,
                         'People removed from the contagious population', color='blue')
removed_pop_v = Population('removed_v', 0,
                           'People removed from the variant contagious population', color='darkslateblue')
removed_pop_w = Population('removed_w', 0,
                           'People removed from the wariant contagious population', color='dodgerblue')
removed_pop_x = Population('removed_x', 0,
                           'People removed from the xariant contagious population', color='lightskyblue')

removed_frac = Parameter('removed_frac', 1., 0., 1.,
                         'fraction of contagious people eventually removed == 1')
removed_delay_pars = {
    'mean': Parameter('removed_delay_mean', 7., 0., 20.,
                      'mean time from becoming contagious to being removed from contagious', hidden=False),
    'sigma': Parameter('removed_delay_sigma', 3., 0.1, 20.,
                       'std dev of times from becoming contagious to being removed from contagious', hidden=False)
}
removed_delay = Delay('removed_delay', 'norm', removed_delay_pars, bc_model)

bc_model.add_connector(
    Propagator('contagious to removed', contagious_pop,
               removed_pop, removed_frac, removed_delay))

bc_model.add_connector(
    Subtractor('removal from contagious', contagious_pop, removed_pop))

# same for variant

removed_delay_pars_v = {
    'mean': Parameter('removed_delay_mean_v', 7., 0., 20.,
                      'variant mean time from becoming contagious to being removed from contagious', hidden=True),
    'sigma': Parameter('removed_delay_sigma_v', 3., 0.1, 20.,
                       'variant std dev of times from becoming contagious to being removed from contagious',
                       hidden=True)
}

removed_delay_v = Delay('removed_delay_v', 'gamma', removed_delay_pars_v, bc_model)

bc_model.add_connector(
    Propagator('contagious_v to removed_v', contagious_pop_v,
               removed_pop_v, removed_frac, removed_delay_v))

bc_model.add_connector(
    Subtractor('removal from contagious_v', contagious_pop_v, removed_pop_v))

# now for wariant

removed_delay_pars_w = {
    'mean': Parameter('removed_delay_mean_w', 7., 0., 20.,
                      'variant mean time from becoming contagious to being removed from contagious', hidden=True),
    'sigma': Parameter('removed_delay_sigma_w', 3., 0.1, 20.,
                       'variant std dev of times from becoming contagious to being removed from contagious',
                       hidden=True)
}

removed_delay_w = Delay('removed_delay_w', 'gamma', removed_delay_pars_w, bc_model)

bc_model.add_connector(
    Propagator('contagious_w to removed_w', contagious_pop_w,
               removed_pop_w, removed_frac, removed_delay_w))

bc_model.add_connector(
    Subtractor('removal from contagious_w', contagious_pop_w, removed_pop_w))

# now for xariant

removed_delay_pars_x = {
    'mean': Parameter('removed_delay_mean_x', 7., 0., 20.,
                      'variant mean time from becoming contagious to being removed from contagious', hidden=True),
    'sigma': Parameter('removed_delay_sigma_x', 3., 0.1, 20.,
                       'variant std dev of times from becoming contagious to being removed from contagious',
                       hidden=True)
}

removed_delay_x = Delay('removed_delay_x', 'gamma', removed_delay_pars_x, bc_model)

bc_model.add_connector(
    Propagator('contagious_x to removed_x', contagious_pop_x,
               removed_pop_x, removed_frac, removed_delay_x))

bc_model.add_connector(
    Subtractor('removal from contagious_x', contagious_pop_x, removed_pop_x))

# Note that there is no removal from bt contagious... that population is not part of the
# infection cycle. All incoming to bt contagious need to be passed to contagious.
# As a result, the bt contagious will be cumulative, not current.

# The symptomatic are sampled by two independent paths:  testing and hospitalization

# TESTING - REPORTING
# propagate a fraction of symptomatics to population with positive test results
# the reported population includes the reporting anomalies

positive_pop = Population('positives', 0,
                          'People who received a positive test report',
                          color='cyan')

reported_pop = Population('reported', 0,
                          'Positives and reporting anomalies',
                          hidden=False, color='forestgreen', show_sim=True,
                          report_noise=True, report_noise_par=report_noise_par,
                          report_backlog_par=report_backlog_par, report_days=report_days)

reported_fraction = Parameter('reported_frac', 0.8, 0., 1.,
                              'fraction of symptomatic people who will ' + \
                              'receive a positive report')
reported_delay_pars = {
    'mean': Parameter('reported_delay_mean', 3., 0., 50.,
                      'mean time from becoming symptomatic to getting positive report', hidden=False),
    'sigma': Parameter('reported_delay_sigma', 1., 0.01, 20.,
                       'standard deviation of times from having symptoms to getting positive report')
}
reported_delay = Delay('reported_delay', 'gamma', reported_delay_pars, bc_model)

bc_model.add_connector(
    Propagator('testing', symptomatic_pop, positive_pop, reported_fraction, reported_delay))

# include those being contact_traced as getting a positive report
bc_model.add_connector(
    Adder('report contact_traced', contact_traced_pop, positive_pop))

# two sources for reported_pop: normal positives and reporting anomalies

bc_model.add_connector(
    Adder('positives reported', positive_pop, reported_pop))

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

# TESTING - REPORTING for breakthroughs
# propagate a fraction of symptomatics to reported population
# no contact tracing or report anomalies -> therefore no "positives" intermediate population
# use same reporting fraction as for non-bt... if necessary, can reduce by reducing symptomatic fraction

bt_reported_pop = Population('bt reported', 0,
                             'Positives and reporting anomalies for bt',
                             hidden=True, color=complementary_color('forestgreen'), show_sim=True,
                             report_noise=True, report_noise_par=report_noise_par,
                             report_backlog_par=report_backlog_par, report_days=report_days)

bc_model.add_connector(
    Propagator('bt_testing', bt_symptomatic_pop, bt_reported_pop, reported_fraction, reported_delay))

# SYMPTOMS -> HOSPITALIZATION
#
# To ensure independence, keep track of the different types of hospitalization
# as separate streams. Since the total probability for entering hospital is low,
# the independent treatment is approximately correct
#
# 2 independent hospitalized streams:
# -> non-ICU hospitalized
# -> ICU
#
# Some in the ICU will get ventilated


non_icu_hospitalized_pop = Population('non_icu_hospitalized', 0,
                                      'Total non_icu hospitalization cases', color='dimgrey', show_sim=False)
non_icu_hospitalized_fraction = Parameter('non_icu_hosp_frac', 0.2, 0., 1.,
                                          'fraction of those with symptoms who will ' + \
                                          'be admitted to non_icu hospital', hidden=False)
non_icu_hospitalized_delay_pars = {
    'mean': Parameter('non_icu_hosp_delay_mean', 6., 0., 50.,
                      'mean time from symptoms to non_icu hospitalization', hidden=False),
    'sigma': Parameter('non_icu_hosp_delay_sigma', 3., 0.01, 20.,
                       'standard deviation of times from symptoms to non_icu hospitalization')
}
non_icu_hospitalized_delay = Delay('non_icu_hosp_delay', 'gamma', non_icu_hospitalized_delay_pars, bc_model)

bc_model.add_connector(
    Propagator('symptomatic to non_icu hospital', symptomatic_pop,
               non_icu_hospitalized_pop, non_icu_hospitalized_fraction, non_icu_hospitalized_delay))

# Those with non_icu hospitalization get released eventually

non_icu_released_pop = Population('non_icu_rel', 0,
                                  'Hospitalized not needing ICU and released',
                                  hidden=True, color='hotpink', show_sim=False)

non_icu_delay_pars = {
    'mean': Parameter('non_icu_rel_delay_mean', 10., 0., 50.,
                      'mean time from non_icu hospital to release', hidden=False),
    'sigma': Parameter('non_icu_rel_delay_sigma', 3., 0.01, 20.,
                       'standard deviation of times from non_icu hospital to release')
}
non_icu_delay = Delay('non_icu_rel_delay', 'gamma', non_icu_delay_pars, bc_model)

release_fraction = Parameter('release_frac', 1., 1., 1.,
                             'fraction for all released == 1', hidden=True)

bc_model.add_connector(
    Propagator('non_icu hospital to released', non_icu_hospitalized_pop,
               non_icu_released_pop, release_fraction, non_icu_delay))

# SYMPTOMS -> HOSPITALIZATION for breakthroughs

bt_non_icu_hospitalized_pop = Population('bt_non_icu_hospitalized', 0,
                                         'Breakthrough non_icu hospitalization cases',
                                         color=complementary_color('dimgrey'), show_sim=False)
bt_non_icu_hospitalized_fraction = Parameter('bt_non_icu_hosp_frac', 0.2, 0., 1.,
                                             'fraction of breakthroughs with symptoms who will ' + \
                                             'be admitted to non_icu hospital', hidden=False)

bc_model.add_connector(
    Propagator('bt symptomatic to non_icu hospital', bt_symptomatic_pop,
               bt_non_icu_hospitalized_pop, bt_non_icu_hospitalized_fraction, non_icu_hospitalized_delay))

# Those with non_icu hospitalization get released eventually

bt_non_icu_released_pop = Population('bt_non_icu_rel', 0,
                                     'Hospitalized breakthrough not needing ICU and released',
                                     hidden=True, color=complementary_color('hotpink'), show_sim=False)

bc_model.add_connector(
    Propagator('bt non_icu hospital to released', bt_non_icu_hospitalized_pop,
               bt_non_icu_released_pop, release_fraction, non_icu_delay))

# Symptoms -> ICU admission
# Keep track of how many currently in ICU
##########################################

icu_pop = Population('icu admissions', 0,
                     'People admitted to ICU', hidden=True, color='deeppink', show_sim=True)

to_icu_fraction = Parameter('icu_frac', 0., 0., 1.,
                            'fraction of symptomatic people who go to ' +
                            'icu', hidden=False)
to_icu_delay_pars = {
    'mean': Parameter('to_icu_delay_mean', 4., 0., 50.,
                      'mean time from symptoms to icu', hidden=False),
    'sigma': Parameter('to_icu_delay_sigma', 2., 0.01, 20.,
                       'standard deviation of times from symptoms to icu')
}
to_icu_delay = Delay('to_icu_delay', 'gamma', to_icu_delay_pars, bc_model)

bc_model.add_connector(
    Propagator('symptomatic to icu', symptomatic_pop,
               icu_pop, to_icu_fraction, to_icu_delay))

# breakthrough ICU admissions
# ----------------------------

bt_icu_pop = Population('bt icu admissions', 0,
                        'People with breakthrough infections admitted to ICU',
                        hidden=True, color=complementary_color('deeppink'), show_sim=True)

bt_to_icu_fraction = Parameter('bt_icu_frac', 0., 0., 1.,
                               'fraction of symptomatic bt people who go to ' +
                               'icu', hidden=False)

bc_model.add_connector(
    Propagator('bt symptomatic to icu', bt_symptomatic_pop,
               bt_icu_pop, bt_to_icu_fraction, to_icu_delay))

# ICU -> VENTILATOR
###################

ventilated_pop = Population('ventilated', 0,
                            'People who received ICU ventilator', color='mediumorchid')

icu_vent_fraction = Parameter('vent_frac', 0.3, 0., 1.,
                              'fraction of those in ICU who need ventillation')

to_vent_delay_pars = {
    'mean': Parameter('to_vent_delay_mean', 4., 0., 50.,
                      'mean time from icu admission to ventilator', hidden=False),
    'sigma': Parameter('to_vent_delay_sigma', 2., 0.01, 20.,
                       'standard deviation of times from icu admission to ventilator')
}
to_vent_delay = Delay('to_vent_delay', 'gamma', to_vent_delay_pars, bc_model)

non_ventilated_rel_pop = Population('non_ventilated_rel', 0,
                                    'ICU non-vent released', color='palevioletred')

non_vent_icu_delay_pars = {
    'mean': Parameter('non_vent_icu_delay_mean', 14., 0., 50.,
                      'mean time from non-vent icu admission to release', hidden=False),
    'sigma': Parameter('non_vent_icu_delay_sigma', 5., 0.01, 20.,
                       'standard deviation of times from non-vent icu admission to release')
}
non_vent_icu_delay = Delay('in_icu_delay', 'gamma', non_vent_icu_delay_pars, bc_model)

bc_model.add_connector(
    Splitter('ventilator', icu_pop, [ventilated_pop, non_ventilated_rel_pop],
             [icu_vent_fraction], [to_vent_delay, non_vent_icu_delay]))

# VENTILATOR -> RELEASED
########################

in_vent_delay_pars = {
    'mean': Parameter('in_vent_delay_mean', 10., 0., 50.,
                      'mean time from ventilator admission to departure', hidden=False),
    'sigma': Parameter('in_vent_delay_sigma', 5., 0.01, 20.,
                       'standard deviation of times from ventilator admission to departure')
}
in_vent_delay = Delay('in_vent_delay', 'gamma', in_vent_delay_pars, bc_model)

ventilated_rel_pop = Population('ventilated_rel', 0,
                                'ICU ventilated released', color='aqua')

vent_rel_fraction = Parameter('vent_rel_frac', 1., 1., 1.,
                              'fraction of those on ventialators eventually released == 1')

bc_model.add_connector(
    Propagator('ventilator to released', ventilated_pop,
               ventilated_rel_pop, vent_rel_fraction, in_vent_delay))

# No ventilation populations for breakthroughs
# All bt ICU releases are non-ventilated

bt_non_ventilated_rel_pop = Population('bt_non_ventilated_rel', 0,
                                       'ICU non-vent released', color=complementary_color('palevioletred'))

bc_model.add_connector(
    Propagator('bt icu released', bt_icu_pop, bt_non_ventilated_rel_pop,
               vent_rel_fraction, non_vent_icu_delay))

# Need new populations to track total number in hospital (non_icu + icu admissions)
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

bc_model.add_connector(
    Adder('include non_icu in hospitalized', non_icu_hospitalized_pop, hospitalized_pop))
bc_model.add_connector(
    Adder('include icu in hospitalized', icu_pop, hospitalized_pop))

# breakthroughs

bt_hospitalized_pop = Population('bt hospitalized', 0,
                                 'Breakthough hospitalization cases',
                                 color=complementary_color('slategrey'), show_sim=True,
                                 report_noise=True, report_noise_par=hosp_noise_par,
                                 report_backlog_par=hosp_backlog_par, report_days=hosp_report_days)

bc_model.add_connector(
    Adder('include bt non_icu in hospitalized', bt_non_icu_hospitalized_pop, bt_hospitalized_pop))
bc_model.add_connector(
    Adder('include bt icu in hospitalized', bt_icu_pop, bt_hospitalized_pop))

# make a copy of hospital admissions to keep track of how many remain in hospital
# ooooooooooooooooooooooooooooooooooooooooooooooooooooooooo

in_hospital_pop = Population('in_hospital', 0,
                             'People currently in hospital',
                             hidden=False, color='darkcyan', show_sim=True)
bc_model.add_connector(
    Adder('copy hospitalizations', hospitalized_pop, in_hospital_pop))

bt_in_hospital_pop = Population('bt_in_hospital', 0,
                                'Breakthroughs currently in hospital',
                                hidden=True, color=complementary_color('darkcyan'), show_sim=True)
bc_model.add_connector(
    Adder('copy bt hospitalizations', bt_hospitalized_pop, bt_in_hospital_pop))

# make a copy of icu admissions to keep track of how many remain in icu
# ooooooooooooooooooooooooooooooooooooooooooooooooooooooooo

in_icu_pop = Population('in_icu', 0,
                        'People currently in ICU', hidden=False, color='hotpink', show_sim=True)

bc_model.add_connector(
    Adder('copy icu admissions', icu_pop, in_icu_pop))

bt_in_icu_pop = Population('bt_in_icu', 0,
                           'Breakthroughs currently in ICU', hidden=True,
                           color=complementary_color('hotpink'), show_sim=True)

bc_model.add_connector(
    Adder('copy bt icu admissions', bt_icu_pop, bt_in_icu_pop))

# make a copy of ventilator admissions to keep track of how many remain on ventilator
# ooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
on_ventilator_pop = Population('on_ventilator', 0,
                               'People currently on ICU ventilator', color='mediumpurple', show_sim=True)

bc_model.add_connector(
    Adder('copy ventilator admissions', ventilated_pop, on_ventilator_pop))

# do subtractions to make the in_hospital etc correct
# oooooooooooooooooooooooooooooooooooooooooooooooooooo

bc_model.add_connector(
    Subtractor('remove non-ICU released', in_hospital_pop, non_icu_released_pop))

bc_model.add_connector(
    Subtractor('remove non-vent released from icu', in_icu_pop, non_ventilated_rel_pop))
bc_model.add_connector(
    Subtractor('remove non-vent released from hospital', in_hospital_pop, non_ventilated_rel_pop))

bc_model.add_connector(
    Subtractor('remove vent released from on_ventilator', on_ventilator_pop, ventilated_rel_pop))
bc_model.add_connector(
    Subtractor('remove vent released from icu', in_icu_pop, ventilated_rel_pop))
bc_model.add_connector(
    Subtractor('remove vent released from hospital', in_hospital_pop, ventilated_rel_pop))

# breakthroughs
# -------------

bc_model.add_connector(
    Subtractor('remove bt non-ICU released', bt_in_hospital_pop, bt_non_icu_released_pop))

bc_model.add_connector(
    Subtractor('remove bt non-vent released from icu', bt_in_icu_pop, bt_non_ventilated_rel_pop))
bc_model.add_connector(
    Subtractor('remove bt non-vent released from hospital', bt_in_hospital_pop, bt_non_ventilated_rel_pop))

# adjust other populations as required
# ooooooooooooooooooooooooooooooooooooo

bc_model.add_connector(
    Subtractor('subtract deaths from total', total_pop, deaths_pop))

bc_model.add_connector(
    Subtractor('subtract infected from susceptible', susceptible_pop, infected_pop))

bc_model.add_connector(
    Subtractor('subtract infected_v from susceptible', susceptible_pop, infected_pop_v))

bc_model.add_connector(
    Subtractor('subtract infected_w from susceptible', susceptible_pop, infected_pop_w))

bc_model.add_connector(
    Subtractor('subtract infected_x from susceptible', susceptible_pop, infected_pop_x))

# breakthroughs (infections reducing bt susceptible)

bc_model.add_connector(
    Subtractor('subtract bt infected from bt susceptible', bt_susceptible_pop, bt_infected_pop))

bc_model.add_connector(
    Subtractor('subtract bt infected_v from bt susceptible', bt_susceptible_pop, bt_infected_pop_v))

bc_model.add_connector(
    Subtractor('subtract bt infected_w from bt susceptible', bt_susceptible_pop, bt_infected_pop_w))

bc_model.add_connector(
    Subtractor('subtract bt infected_x from sbt usceptible', bt_susceptible_pop, bt_infected_pop_x))

# some with a positive test will decide not to get vaccinated
frac_report_novacc_par = Parameter('frac_report_novacc', 1.0, 0., 1.,
                                   'Fraction of those reported who decide against vaccination')

bc_model.add_connector(
    Subtractor('subtract vaccinated from vacc cand', vaccan_pop, vaccinated_pop,
               scale_factor=frac_report_novacc_par))

bc_model.add_connector(
    Subtractor('subtract usefully vaccinated from the sus vacc cand', susvaccan_pop, usefully_vaccinated_pop))

bc_model.add_connector(
    Subtractor('subtract immunized from susceptible', susceptible_pop, immunized_pop))

bc_model.add_connector(
    Subtractor('subtract some of infected from sus vacc cand', susvaccan_pop, infected_pop,
               ratio_populations=[susvaccan_pop, susceptible_pop]))

bc_model.add_connector(
    Subtractor('subtract some of infected_v from sus vacc cand', susvaccan_pop, infected_pop_v,
               ratio_populations=[susvaccan_pop, susceptible_pop]))

bc_model.add_connector(
    Subtractor('subtract some of infected_w from sus vacc cand', susvaccan_pop, infected_pop_w,
               ratio_populations=[susvaccan_pop, susceptible_pop]))

bc_model.add_connector(
    Subtractor('subtract some of infected_x from sus vacc cand', susvaccan_pop, infected_pop_x,
               ratio_populations=[susvaccan_pop, susceptible_pop]))

rep_frac_pop = Population('frac', 0.)

bc_model.add_connector(
    Adder('calc rep frac', reported_pop, rep_frac_pop, ratio_populations=[vaccinated_pop, total_pop]))

bc_model.add_connector(
    Adder('add vaccinated reported', rep_frac_pop, vaccan_pop))

bc_model.add_connector(
    Subtractor('subtract immunized reported', vaccan_pop, rep_frac_pop, scale_factor=vaccine_effectiveness))

bc_model.add_connector(
    Subtractor('subtract reported', vaccan_pop, reported_pop))

# transitional aspects of the model
# oooooooooooooooooooooooooooooooooo

n_rate_transitions = 14
n_rt_visible = 4

times = []
alphas = [Parameter('alpha_0', 0.4, 0., 2.,
                    'initial transmission rate', hidden=False)]
for i in range(n_rate_transitions):
    j = i + 1
    hidden = j > n_rt_visible
    enabled = i == 0
    times.append(Parameter('trans_rate_' + str(j) + '_time', 20 + i, 0, 600,
                           'day of trans rate change ' + 'str(j)',
                           parameter_type='int', hidden=hidden))
    alphas.append(Parameter('alpha_' + str(j), 0.1, 0., 2.,
                            'alpha after transition ' + str(j), hidden=hidden))
    bc_model.add_transition(
        Modifier('trans_rate_' + str(j), 'rel_days', times[i], trans_rate,
                 alphas[i], alphas[j], enabled=enabled, model=bc_model))

n_rate_transitions = 6
n_rt_visible = 0

times_v = []
alphas_v = [Parameter('alpha_0_v', 0.4, 0., 2.,
                      'initial variant transmission rate', hidden=False)]
for i in range(n_rate_transitions):
    j = i + 1
    hidden = j > n_rt_visible
    enabled = False
    times_v.append(Parameter('trans_rate_' + str(j) + '_time_v', 20 + i, 0, 600,
                             'day of variant trans rate change ' + 'str(j)',
                             parameter_type='int', hidden=hidden))
    alphas_v.append(Parameter('alpha_' + str(j) + '_v', 0.1, 0., 2.,
                              'alpha_v after transition ' + str(j), hidden=hidden))
    bc_model.add_transition(
        Modifier('trans_rate_' + str(j) + '_v', 'rel_days', times_v[i], trans_rate_v,
                 alphas_v[i], alphas_v[j], enabled=enabled, model=bc_model))

n_rate_transitions = 12

times_w = []
alphas_w = [Parameter('alpha_0_w', 0.4, 0., 2.,
                      'initial wariant transmission rate', hidden=False)]
for i in range(n_rate_transitions):
    j = i + 1
    hidden = j > n_rt_visible
    enabled = False
    times_w.append(Parameter('trans_rate_' + str(j) + '_time_w', 200 + i, 0, 600,
                             'day of variant trans rate change ' + 'str(j)',
                             parameter_type='int', hidden=hidden))
    alphas_w.append(Parameter('alpha_' + str(j) + '_w', 0.1, 0., 2.,
                              'alpha_w after transition ' + str(j), hidden=hidden))
    bc_model.add_transition(
        Modifier('trans_rate_' + str(j) + '_w', 'rel_days', times_w[i], trans_rate_w,
                 alphas_w[i], alphas_w[j], enabled=enabled, model=bc_model))

n_rate_transitions = 4

times_x = []
alphas_x = [Parameter('alpha_0_x', 0.4, 0., 2.,
                      'initial xariant transmission rate', hidden=False)]
for i in range(n_rate_transitions):
    j = i + 1
    hidden = j > n_rt_visible
    enabled = False
    times_x.append(Parameter('trans_rate_' + str(j) + '_time_x', 200 + i, 0, 600,
                             'day of variant trans rate change ' + 'str(j)',
                             parameter_type='int', hidden=hidden))
    alphas_x.append(Parameter('alpha_' + str(j) + '_x', 0.1, 0., 2.,
                              'alpha_x after transition ' + str(j), hidden=hidden))
    bc_model.add_transition(
        Modifier('trans_rate_' + str(j) + '_x', 'rel_days', times_x[i], trans_rate_x,
                 alphas_x[i], alphas_x[j], enabled=enabled, model=bc_model))

# contact tracing

trans_traced_1_time = Parameter('trans_trace_1_time', 100, 0, 600,
                                'number of days before contact traced fraction changes',
                                parameter_type='int', hidden=False)

trans_traced_0 = Parameter('trans_traced_0', 0., 0., 1.,
                           'initial contact traced fraction')

trans_traced_1 = Parameter('trans_traced_1', 0.1, 0., 1.,
                           'contact traced fraction after transition', hidden=False)

bc_model.add_transition(
    Modifier('trans_traced_1', 'rel_days', trans_traced_1_time, contact_traced_fraction,
             trans_traced_0, trans_traced_1, enabled=False, model=bc_model))

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

# outbreaks

outbreak_pop = Population('outbreaks', 0,
                          'Infection outbreaks')

outbreak_1_time = Parameter('outbreak_1_time', 14, 0, 600,
                            'number of days since t0 when outbreak_1 established',
                            parameter_type='int', hidden=False)

outbreak_1_number = Parameter('outbreak_1_number', 10., 0., 50000.,
                              'number of infections in outbreak_1',
                              hidden=False)

bc_model.add_transition(
    Injector('outbreak_1', 'rel_days', outbreak_1_time, outbreak_pop,
             outbreak_1_number, enabled=False, model=bc_model))

outbreak_fraction = Parameter('outbreak_frac', 1., 0., 1.,
                              'fraction of infected in outbreak active ==1')
outbreak_delay_pars = {
    'mean': Parameter('outbreak_delay_mean', 7., 0., 50.,
                      'mean delay time for outbreak', hidden=False),
    'sigma': Parameter('outbreak_delay_sigma', 1., 0.01, 20.,
                       'standard deviation of outbreak times',
                       hidden=False)
}
outbreak_delay = Delay('outbreak_delay', 'gamma', outbreak_delay_pars, bc_model)

bc_model.add_connector(
    Propagator('outbreaks to infected', outbreak_pop, infected_pop,
               outbreak_fraction, outbreak_delay))

outbreak_2_time = Parameter('outbreak_2_time', 21, 0, 600,
                            'number of days since t0 when outbreak_2 established',
                            parameter_type='int', hidden=False)

outbreak_2_number = Parameter('outbreak_2_number', 10., 0., 50000.,
                              'number of infections in outbreak_2',
                              hidden=False)

bc_model.add_transition(
    Injector('outbreak_2', 'rel_days', outbreak_2_time, outbreak_pop,
             outbreak_2_number, enabled=False, model=bc_model))

outbreak_3_time = Parameter('outbreak_3_time', 41, 0, 600,
                            'number of days since t0 when outbreak_3 established',
                            parameter_type='int', hidden=False)

outbreak_3_number = Parameter('outbreak_3_number', 10., 0., 50000.,
                              'number of infections in outbreak_3',
                              hidden=False)

bc_model.add_transition(
    Injector('outbreak_3', 'rel_days', outbreak_3_time, outbreak_pop,
             outbreak_3_number, enabled=False, model=bc_model))

outbreak_4_time = Parameter('outbreak_4_time', 61, 0, 600,
                            'number of days since t0 when outbreak_4 established',
                            parameter_type='int')

outbreak_4_number = Parameter('outbreak_4_number', 10., 0., 50000.,
                              'number of infections in outbreak_4')

bc_model.add_transition(
    Injector('outbreak_4', 'rel_days', outbreak_4_time, outbreak_pop,
             outbreak_4_number, enabled=False, model=bc_model))

outbreak_5_time = Parameter('outbreak_5_time', 61, 0, 600,
                            'number of days since t0 when outbreak_5 established',
                            parameter_type='int')

outbreak_5_number = Parameter('outbreak_5_number', 10., 0., 50000.,
                              'number of infections in outbreak_5')

bc_model.add_transition(
    Injector('outbreak_5', 'rel_days', outbreak_5_time, outbreak_pop,
             outbreak_5_number, enabled=False, model=bc_model))

outbreak_6_time = Parameter('outbreak_6_time', 61, 0, 600,
                            'number of days since t0 when outbreak_6 established',
                            parameter_type='int')

outbreak_6_number = Parameter('outbreak_6_number', 10., 0., 50000.,
                              'number of infections in outbreak_6')

bc_model.add_transition(
    Injector('outbreak_6', 'rel_days', outbreak_6_time, outbreak_pop,
             outbreak_6_number, enabled=False, model=bc_model))

# variant outbreak: necessary to introduce infections...

outbreak_pop_v = Population('outbreaks_v', 0,
                            'Variant infection outbreaks')

outbreak_v_time = Parameter('outbreak_v_time', 14, 0, 600,
                            'number of days since t0 when outbreak_v established',
                            parameter_type='int', hidden=True)

outbreak_v_number = Parameter('outbreak_v_number', 10., 0., 50000.,
                              'number of infections in outbreak_v',
                              hidden=True)

bc_model.add_transition(
    Injector('outbreak_v', 'rel_days', outbreak_v_time, outbreak_pop_v,
             outbreak_v_number, enabled=False, model=bc_model))

outbreak_v_delay_pars = {
    'mean': Parameter('outbreak_v_delay_mean', 7., 0., 50.,
                      'mean delay time for outbreak_v'),
    'sigma': Parameter('outbreak_v_delay_sigma', 1., 0.01, 20.,
                       'standard deviation of outbreak_v times')
}
outbreak_v_delay = Delay('outbreak_v_delay', 'gamma', outbreak_v_delay_pars, bc_model)

bc_model.add_connector(
    Propagator('outbreaks to infected_v', outbreak_pop_v, infected_pop_v,
               outbreak_fraction, outbreak_v_delay))

outbreak_v2_time = Parameter('outbreak_v2_time', 210, 0, 600,
                             'number of days since t0 when outbreak_v2 established',
                             parameter_type='int', hidden=False)

outbreak_v2_number = Parameter('outbreak_v2_number', 10., 0., 50000.,
                               'number of infections in outbreak_v2',
                               hidden=False)

bc_model.add_transition(
    Injector('outbreak_v2', 'rel_days', outbreak_v2_time, outbreak_pop_v,
             outbreak_v2_number, enabled=False, model=bc_model))

# wariant outbreak: necessary to introduce infections...

outbreak_pop_w = Population('outbreaks_w', 0,
                            'Wariant infection outbreaks')

outbreak_w_time = Parameter('outbreak_w_time', 140, 0, 600,
                            'number of days since t0 when outbreak_w established',
                            parameter_type='int', hidden=True)

outbreak_w_number = Parameter('outbreak_w_number', 10., 0., 50000.,
                              'number of infections in outbreak_w',
                              hidden=True)

bc_model.add_transition(
    Injector('outbreak_w', 'rel_days', outbreak_w_time, outbreak_pop_w,
             outbreak_w_number, enabled=False, model=bc_model))

outbreak_w_delay_pars = {
    'mean': Parameter('outbreak_w_delay_mean', 7., 0., 50.,
                      'mean delay time for outbreak_w'),
    'sigma': Parameter('outbreak_w_delay_sigma', 1., 0.01, 20.,
                       'standard deviation of outbreak_w times')
}
outbreak_w_delay = Delay('outbreak_w_delay', 'gamma', outbreak_w_delay_pars, bc_model)

bc_model.add_connector(
    Propagator('outbreaks to infected_w', outbreak_pop_w, infected_pop_w,
               outbreak_fraction, outbreak_w_delay))

outbreak_w2_time = Parameter('outbreak_w2_time', 210, 0, 600,
                             'number of days since t0 when outbreak_w2 established',
                             parameter_type='int', hidden=False)

outbreak_w2_number = Parameter('outbreak_w2_number', 10., 0., 50000.,
                               'number of infections in outbreak_w2',
                               hidden=False)

bc_model.add_transition(
    Injector('outbreak_w2', 'rel_days', outbreak_w2_time, outbreak_pop_w,
             outbreak_w2_number, enabled=False, model=bc_model))

outbreak_w3_time = Parameter('outbreak_w3_time', 510, 0, 900,
                             'number of days since t0 when outbreak_w3 established',
                             parameter_type='int', hidden=False)

outbreak_w3_number = Parameter('outbreak_w3_number', 10., 0., 50000.,
                               'number of infections in outbreak_w3',
                               hidden=False)

bc_model.add_transition(
    Injector('outbreak_w3', 'rel_days', outbreak_w3_time, outbreak_pop_w,
             outbreak_w3_number, enabled=False, model=bc_model))

# xariant outbreak: necessary to introduce infections...

outbreak_pop_x = Population('outbreaks_x', 0,
                            'xariant infection outbreaks')

outbreak_x_time = Parameter('outbreak_x_time', 140, 0, 600,
                            'number of days since t0 when outbreak_x established',
                            parameter_type='int', hidden=True)

outbreak_x_number = Parameter('outbreak_x_number', 10., 0., 50000.,
                              'number of infections in outbreak_x',
                              hidden=True)

bc_model.add_transition(
    Injector('outbreak_x', 'rel_days', outbreak_x_time, outbreak_pop_x,
             outbreak_x_number, enabled=False, model=bc_model))

outbreak_x_delay_pars = {
    'mean': Parameter('outbreak_x_delay_mean', 7., 0., 50.,
                      'mean delay time for outbreak_x'),
    'sigma': Parameter('outbreak_x_delay_sigma', 1., 0.01, 20.,
                       'standard deviation of outbreak_x times')
}
outbreak_x_delay = Delay('outbreak_x_delay', 'gamma', outbreak_x_delay_pars, bc_model)

bc_model.add_connector(
    Propagator('outbreaks to infected_x', outbreak_pop_x, infected_pop_x,
               outbreak_fraction, outbreak_x_delay))

outbreak_x2_time = Parameter('outbreak_x2_time', 210, 0, 600,
                             'number of days since t0 when outbreak_x2 established',
                             parameter_type='int', hidden=False)

outbreak_x2_number = Parameter('outbreak_x2_number', 10., 0., 50000.,
                               'number of infections in outbreak_x2',
                               hidden=False)

bc_model.add_transition(
    Injector('outbreak_x2', 'rel_days', outbreak_x2_time, outbreak_pop_x,
             outbreak_x2_number, enabled=False, model=bc_model))

# Assume outbreaks are limited to unvaccinated populations for now

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

reported_fractions = []
reported_fractions.append(Parameter('report_frac_0', 0.95, 0., 1.,
                                    'initial reported fraction', hidden=False))

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

bc_model.boot_setup(contagious_pop, 0.1,
                    exclusion_populations=[total_pop, susceptible_pop, susvaccan_pop, vaccan_pop])

bc_model.save_file('ref_model_3_0.pypm')
