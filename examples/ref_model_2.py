# -*- coding: utf-8 -*-
"""
pyPM.ca reference model #2

To make the critical parameters less dependent of each other, break the measured
populations into separate paths.

Reduced the lag time between contagious and reported, (as compared to reference model #1)
to better match the situation in BC and elsewhere - transitions now occur where they should

Add a contact tracing population, to allow another mechanism to remove people from
the contagious population early. Initial value has p=0, so it does not act until it
is switched on.

@author: karlen
"""

from pypmca import Model, Population, Delay, Parameter, Multiplier, Propagator, \
    Splitter, Adder, Subtractor, Chain, Modifier, Injector

# Test by building a population model for BC

bc_model = Model('ref_model_2')
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

initial_contagious_par = Parameter('cont_0', 55., 0., 5000.,
                                   'Number of contagious people at t0',
                                   hidden=False)

contagious_pop = Population('contagious', initial_contagious_par,
                            'number of people that can cause someone to become infected',
                            hidden=False, color='red')

# this value is only used if the transition is removed
trans_rate = Parameter('alpha', 0.390, 0., 2.,
                       'mean number of people that a contagious person infects ' +
                       'per day', hidden=True)
infection_delay = Delay('fast', 'fast', model=bc_model)

bc_model.add_connector(
    Multiplier('infection cycle', [susceptible_pop, contagious_pop, total_pop],
               infected_pop, trans_rate, infection_delay, bc_model))

contagious_frac = Parameter('cont_frac', 0.9, 0., 1.,
                            'fraction of infected people that become contagious',
                            hidden=False)
contagious_delay_pars = {
    'mean': Parameter('cont_delay_mean', 2., 0., 50.,
                      'mean time from being infected to becoming contagious'),
    'sigma': Parameter('cont_delay_sigma', 1., 0.01, 20.,
                       'standard deviation of times from being infected to becoming contagious')
}

contagious_delay = Delay('cont_delay', 'norm', contagious_delay_pars, bc_model)

bc_model.add_connector(
    Propagator('infected to contagious', infected_pop,
               contagious_pop, contagious_frac, contagious_delay))

# The contagious either recover or die
# This split is only used to track the deaths.
# The removal of recoveries/positive tests etc done in a separate path
# oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo

recovered_pop = Population('recovered', 0,
                           'People who have recovered from the illness ' +
                           'and are therefore no longer susceptible', color='limegreen')
deaths_pop = Population('deaths', 0,
                        'people who have died from the illness', hidden=False,
                        color='indigo', show_sim=True)
recover_fraction = Parameter('recover_frac', 0.99, 0., 1.,
                             'fraction of infected people who recover', hidden=False)
recover_delay_pars = {
    'mean': Parameter('recover_delay_mean', 14., 0., 50., 'mean time from infection to recovery'),
    'sigma': Parameter('recover_delay_sigma', 5., 0.01, 20.,
                       'standard deviation of times from infection to recovery')
}
recover_delay = Delay('recover_delay', 'norm', recover_delay_pars, bc_model)

death_delay_pars = {
    'mean': Parameter('death_delay_mean', 21., 0., 50.,
                      'mean time from infection to death', hidden=False),
    'sigma': Parameter('death_delay_sigma', 5., 0.01, 20.,
                       'standard deviation of times from infection to death')
}
death_delay = Delay('death_delay', 'norm', death_delay_pars, bc_model)

bc_model.add_connector(
    Splitter('recovery', contagious_pop, [recovered_pop, deaths_pop],
             [recover_fraction], [recover_delay, death_delay]))

# The newly contagious are split into three groups:
# contact-traced, symptomatic and non-symptomatic.
# contact-traced are removed when then are identified: by adding them to the reported sample
# symptomatic are sampled by a testing path: those receiving a positive test are removed
# When the non-symptomatic recover thay are removed from the contagious population.

contact_traced_pop = Population('contact_traced', 0,
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
contact_traced_delay = Delay('contact_traced_delay', 'norm', contact_traced_delay_pars,
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
symptomatic_delay = Delay('symptomatic_delay', 'norm', symptomatic_delay_pars,
                          bc_model)

asymptomatic_recovered_pop = Population('asymptomatic recovered', 0,
                                        'People who have not shown symptoms', color='silver')

asymptomatic_delay_pars = {
    'mean': Parameter('asymp_rec_delay_mean', 12., 0., 50.,
                      'mean time from becoming contagious to recovery (without symptoms)'),
    'sigma': Parameter('asymp_rec_delay_sigma', 4., 0.01, 20.,
                       'std dev of times from becoming contagious to recovery with no symptoms')
}
asymptomatic_delay = Delay('asymp_rec_delay', 'norm', asymptomatic_delay_pars,
                           bc_model)

bc_model.add_connector(
    Splitter('symptoms', contagious_pop, [contact_traced_pop, symptomatic_pop, asymptomatic_recovered_pop],
             [contact_traced_fraction, symptomatic_fraction],
             [contact_traced_delay, symptomatic_delay, asymptomatic_delay]))

bc_model.add_connector(
    Subtractor('remove asymptomatic recoveries', contagious_pop, asymptomatic_recovered_pop))

# The symptomatic are sampled by two independent paths:  testing and hospitalization

# TESTING - REPORTING
# break into group that receives a positive test result, and those who do not
# (either because they were not tested or had a negative test)
# those receiving a positive test result are removed from contagious immediately
# those not are removed from contagious when they recover

unreported_recovered_pop = Population('unreported', 0,
                                      'Symptomatic people who do not receive a positive test result',
                                      color='deepskyblue')

unreported_delay_pars = {
    'mean': Parameter('unrep_rec_delay_mean', 12., 0., 50.,
                      'mean time from symptoms to recovery for unreported',
                      hidden=False),
    'sigma': Parameter('unrep_rec_delay_sigma', 4., 0.01, 20.,
                       'standard deviation of times from having symptoms to recovery for unreported')
}
unreported_delay = Delay('unreported_delay', 'norm', unreported_delay_pars, bc_model)

reported_pop = Population('reported', 0,
                          'Symptomatic people who received a positive test report',
                          hidden=False, color='forestgreen', show_sim=True)
reported_fraction = Parameter('reported_frac', 0.8, 0., 1.,
                              'fraction of symptomatic people who will ' + \
                              'receive a positive report')
reported_delay_pars = {
    'mean': Parameter('reported_delay_mean', 3., 0., 50.,
                      'mean time from becoming symptomatic to getting positive report', hidden=False),
    'sigma': Parameter('reported_delay_sigma', 1., 0.01, 20.,
                       'standard deviation of times from having symptoms to getting positive report')
}
reported_delay = Delay('reported_delay', 'norm', reported_delay_pars, bc_model)

bc_model.add_connector(
    Splitter('testing', symptomatic_pop, [reported_pop, unreported_recovered_pop],
             [reported_fraction], [reported_delay, unreported_delay]))

# include those being contact_traced as getting a positive report
bc_model.add_connector(
    Adder('report contact_traced', contact_traced_pop, reported_pop))

bc_model.add_connector(
    Subtractor('remove those reported', contagious_pop, reported_pop))

bc_model.add_connector(
    Subtractor('remove unreported recoveries', contagious_pop, unreported_recovered_pop))

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
#
#
# Use the sum of the two streams to be compared to hospitalized data

non_icu_hospitalized_pop = Population('non_icu_hospitalized', 0,
                                      'Total non_icu hospitalization cases', color='dimgrey', show_sim=False)
non_icu_hospitalized_fraction = Parameter('non_icu_hosp_frac', 0.2, 0., 1.,
                                          'fraction of those with symptoms who will ' + \
                                          'be admitted to non_icu hospital', hidden=False)
non_icu_hospitalized_delay_pars = {
    'mean': Parameter('non_icu_hosp_delay_mean', 12., 0., 50.,
                      'mean time from symptoms to non_icu hospitalization', hidden=False),
    'sigma': Parameter('non_icu_hosp_delay_sigma', 3., 0.01, 20.,
                       'standard deviation of times from symptoms to non_icu hospitalization')
}
non_icu_hospitalized_delay = Delay('non_icu_hosp_delay', 'norm', non_icu_hospitalized_delay_pars, bc_model)

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
non_icu_delay = Delay('non_icu_rel_delay', 'norm', non_icu_delay_pars, bc_model)

release_fraction = Parameter('release_frac', 1., 1., 1.,
                             'fraction for all released == 1', hidden=True)

bc_model.add_connector(
    Propagator('non_icu hospital to released', non_icu_hospitalized_pop,
               non_icu_released_pop, release_fraction, non_icu_delay))

# Symptoms -> ICU admission
# Keep track of how many currently in ICU
##########################################

icu_pop = Population('icu admissions', 0,
                     'People admitted to ICU', hidden=True, color='deeppink', show_sim=True)

to_icu_fraction = Parameter('icu_frac', 0.1, 0., 1.,
                            'fraction of symptomatic people who go to ' +
                            'icu', hidden=False)
to_icu_delay_pars = {
    'mean': Parameter('to_icu_delay_mean', 14., 0., 50.,
                      'mean time from symptoms to icu', hidden=False),
    'sigma': Parameter('to_icu_delay_sigma', 1., 0.01, 20.,
                       'standard deviation of times from symptoms to icu')
}
to_icu_delay = Delay('to_icu_delay', 'norm', to_icu_delay_pars, bc_model)

bc_model.add_connector(
    Propagator('symptomatic to icu', symptomatic_pop,
               icu_pop, to_icu_fraction, to_icu_delay))

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
to_vent_delay = Delay('to_vent_delay', 'norm', to_vent_delay_pars, bc_model)

non_ventilated_rel_pop = Population('non_ventilated_rel', 0,
                                    'ICU non-vent released', color='palevioletred')

non_vent_icu_delay_pars = {
    'mean': Parameter('non_vent_icu_delay_mean', 14., 0., 50.,
                      'mean time from non-vent icu admission to release', hidden=False),
    'sigma': Parameter('non_vent_icu_delay_sigma', 5., 0.01, 20.,
                       'standard deviation of times from non-vent icu admission to release')
}
non_vent_icu_delay = Delay('in_icu_delay', 'norm', non_vent_icu_delay_pars, bc_model)

bc_model.add_connector(
    Splitter('ventilator', icu_pop, [ventilated_pop, non_ventilated_rel_pop],
             [to_icu_fraction], [to_icu_delay, non_vent_icu_delay]))

# VENTILATOR -> RELEASED
########################

in_vent_delay_pars = {
    'mean': Parameter('in_vent_delay_mean', 10., 0., 50.,
                      'mean time from ventilator admission to departure', hidden=False),
    'sigma': Parameter('in_vent_delay_sigma', 5., 0.01, 20.,
                       'standard deviation of times from ventilator admission to departure')
}
in_vent_delay = Delay('in_vent_delay', 'norm', in_vent_delay_pars, bc_model)

ventilated_rel_pop = Population('ventilated_rel', 0,
                                'ICU ventilated released', color='aqua')

vent_rel_fraction = Parameter('vent_frac', 1., 1., 1.,
                              'fraction of those on ventialators eventually released == 1')

bc_model.add_connector(
    Propagator('ventilator to released', ventilated_pop,
               ventilated_rel_pop, vent_rel_fraction, in_vent_delay))

# Need new populations to track total number in hospital (non_icu + icu admissions)

hospitalized_pop = Population('hospitalized', 0,
                              'Total hospitalization cases', color='slategrey', show_sim=True)

bc_model.add_connector(
    Adder('include non_icu in hospitalized', non_icu_hospitalized_pop, hospitalized_pop))
bc_model.add_connector(
    Adder('include icu in hospitalized', icu_pop, hospitalized_pop))

# make a copy of hospital admissions to keep track of how many remain in hospital
# ooooooooooooooooooooooooooooooooooooooooooooooooooooooooo

in_hospital_pop = Population('in_hospital', 0,
                             'People currently in hospital',
                             hidden=False, color='darkcyan', show_sim=True)
bc_model.add_connector(
    Adder('copy hospitalizations', hospitalized_pop, in_hospital_pop))

# make a copy of icu admissions to keep track of how many remain in icu
# ooooooooooooooooooooooooooooooooooooooooooooooooooooooooo

in_icu_pop = Population('in_icu', 0,
                        'People currently in ICU', hidden=False, color='deeppink', show_sim=True)

bc_model.add_connector(
    Adder('copy icu admissions', icu_pop, in_icu_pop))

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

# adjust other populations as required
# ooooooooooooooooooooooooooooooooooooo

bc_model.add_connector(
    Subtractor('subtract deaths from total', total_pop, deaths_pop))

bc_model.add_connector(
    Subtractor('subtract infected from susceptible', susceptible_pop, infected_pop))

# transitional aspects of the model
# oooooooooooooooooooooooooooooooooo

trans_rate_1_time = Parameter('trans_rate_1_time', 16, 0, 300,
                              'number of days before 1st transmission rate change',
                              parameter_type='int', hidden=False)

trans_rate_2_time = Parameter('trans_rate_2_time', 50, 0, 300,
                              'number of days before 2nd transmission rate change',
                              parameter_type='int', hidden=False)

trans_rate_3_time = Parameter('trans_rate_3_time', 70, 0, 300,
                              'number of days before 3nd transmission rate change',
                              parameter_type='int', hidden=False)

trans_rate_0 = Parameter('alpha_0', 0.385, 0., 2.,
                         'initial transmission rate', hidden=False)

trans_rate_1 = Parameter('alpha_1', 0.062, 0., 2.,
                         'transmission rate after 1st transition', hidden=False)

trans_rate_2 = Parameter('alpha_2', 0.062, 0., 2.,
                         'transmission rate after 2nd transition', hidden=False)

trans_rate_3 = Parameter('alpha_3', 0.062, 0., 2.,
                         'transmission rate after 3nd transition', hidden=False)

bc_model.add_transition(
    Modifier('trans_rate_1', 'rel_days', trans_rate_1_time, trans_rate,
             trans_rate_0, trans_rate_1, enabled=True, model=bc_model))

bc_model.add_transition(
    Modifier('trans_rate_2', 'rel_days', trans_rate_2_time, trans_rate,
             trans_rate_1, trans_rate_2, enabled=False, model=bc_model))

bc_model.add_transition(
    Modifier('trans_rate_3', 'rel_days', trans_rate_3_time, trans_rate,
             trans_rate_2, trans_rate_3, enabled=False, model=bc_model))

trans_traced_1_time = Parameter('trans_trace_1_time', 100, 0, 300,
                                'number of days before contact traced fraction changes',
                                parameter_type='int', hidden=False)

trans_traced_0 = Parameter('trans_traced_0', 0., 0., 1.,
                           'initial contact traced fraction')

trans_traced_1 = Parameter('trans_traced_1', 0.1, 0., 1.,
                         'contact traced fraction after transition', hidden=False)

bc_model.add_transition(
    Modifier('trans_traced_1', 'rel_days', trans_traced_1_time, contact_traced_fraction,
             trans_traced_0, trans_traced_1, enabled=False, model=bc_model))

outbreak_pop = Population('outbreaks', 0,
                          'Infection outbreaks')

outbreak_1_time = Parameter('outbreak_1_time', 14, 0, 100,
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
                       'standard deviation of loutbreak times',
                       hidden=False)
}
outbreak_delay = Delay('outbreak_delay', 'norm', outbreak_delay_pars, bc_model)

bc_model.add_connector(
    Propagator('outbreaks to infected', outbreak_pop, infected_pop,
               outbreak_fraction, outbreak_delay))

outbreak_2_time = Parameter('outbreak_2_time', 21, 0, 100,
                            'number of days since t0 when outbreak_2 established',
                            parameter_type='int', hidden=False)

outbreak_2_number = Parameter('outbreak_2_number', 10., 0., 50000.,
                              'number of infections in outbreak_2',
                              hidden=False)

bc_model.add_transition(
    Injector('outbreak_2', 'rel_days', outbreak_2_time, outbreak_pop,
             outbreak_2_number, enabled=False, model=bc_model))

outbreak_3_time = Parameter('outbreak_3_time', 41, 0, 100,
                            'number of days since t0 when outbreak_3 established',
                            parameter_type='int', hidden=False)

outbreak_3_number = Parameter('outbreak_3_number', 10., 0., 50000.,
                              'number of infections in outbreak_3',
                              hidden=False)

bc_model.add_transition(
    Injector('outbreak_3', 'rel_days', outbreak_3_time, outbreak_pop,
             outbreak_3_number, enabled=False, model=bc_model))

# define boot parameters
# ooooooooooooooooooooooo

bc_model.boot_setup(contagious_pop, 1,
                    exclusion_populations=[total_pop, susceptible_pop])

bc_model.save_file('ref_model_2.pypm')
