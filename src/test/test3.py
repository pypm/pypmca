# -*- coding: utf-8 -*-
"""
Created on Mon Apr  6 14:29:36 2020

Modified version of the bc model, now with proper treatment of ICU
(allows it to lag from contagious, not after proceeding through the
 long reporting chain)

@author: karlen
"""

from Population import Population
from Model import Model
from Delay import Delay
from Parameter import Parameter
from Multiplier import Multiplier
from Propagator import Propagator
from Splitter import Splitter
from Adder import Adder
from Subtractor import Subtractor
from Chain import Chain
from Modifier import Modifier
from Injector import Injector

# Test by building a population model for BC

bc_model = Model('model v4.x')
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
#oooooooooooooooooooooooooooo

initial_contagious_par = Parameter('cont_0', 55., 0., 5000.,
                                   'Number of contagious people at t0',
                                   hidden=False)

contagious_pop = Population('contagious', initial_contagious_par,
                            'number of people that can cause someone to become infected',
                            hidden=False, color='red')

# this value is only used if the transition is removed
trans_rate = Parameter('alpha', 0.390, 0., 2.,
                       'mean number of people that a contagious person infects '+
                       'per day', hidden=True)
infection_delay = Delay('fast', 'fast', model=bc_model)

bc_model.add_connector(
    Multiplier('infection cycle', [susceptible_pop, contagious_pop, total_pop],
               infected_pop, trans_rate, infection_delay, bc_model))

contagious_frac = Parameter('cont_frac', 0.9, 0., 1.,
                            'fraction of infected people that become contagious',
                            hidden=False)
contagious_delay_pars = {
    'mean': Parameter('cont_delay_mean', 3., 0., 50.,
                      'mean time from being infected to becoming contagious'),
    'sigma': Parameter('cont_delay_sigma', 1., 0.01, 20.,
                       'standard deviation of times from being infected to becoming contagious')
}

contagious_delay = Delay('cont_delay', 'norm', contagious_delay_pars, bc_model)

bc_model.add_connector(
    Propagator('infected to contagious', infected_pop,
               contagious_pop, contagious_frac, contagious_delay))

# The infected either recover or die
#ooooooooooooooooooooooooooooooooooo

recovered_pop = Population('recovered', 0,
                           'People who have recovered from the illness '+
                           'and are therefore no longer susceptible', color='limegreen')
deaths_pop = Population('deaths', 0,
                        'people who have died from the illness', hidden=False,
                        color='indigo', show_sim = True)
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
                      'mean time from infection to death',hidden=False),
    'sigma': Parameter('death_delay_sigma', 5., 0.01, 20.,
                       'standard deviation of times from infection to death')
}
death_delay = Delay('death_delay', 'norm', death_delay_pars, bc_model)

bc_model.add_connector(
    Splitter('recovery', infected_pop, [recovered_pop, deaths_pop],
             [recover_fraction], [recover_delay, death_delay]))

# The newly contagious are split into three groups, with the delays
# reference to the start of being contagious.
# The 3 groups are: ICU patients, non-ICU but symptomatic, non-quarantined,
# the last being a catch all...

icu_pop = Population('icu admissions', 0,
                     'People admitted to ICU', hidden=True, color='deeppink', show_sim=True)

to_icu_fraction = Parameter('icu_frac', 0.05, 0., 1.,
                            'fraction of contagious people who go to '+
                            'icu', hidden=False)
to_icu_delay_pars = {
    'mean': Parameter('to_icu_delay_mean', 12., 0., 50.,
                      'mean time from contagious to icu', hidden=False),
    'sigma': Parameter('to_icu_delay_sigma', 3., 0.01, 20.,
                       'standard deviation of times from contagious to icu')
}
to_icu_delay = Delay('to_icu_delay', 'norm', to_icu_delay_pars, bc_model)


symptomatic_pop = Population('symptomatic', 0,
                             'People who have shown symptoms', color='chocolate')
symptomatic_fraction = Parameter('symptomatic_frac', 0.9, 0., 1.,
                                 'fraction of contagious people who become '+
                                 'symptomatic', hidden=False)
symptomatic_delay_pars = {
    'mean': Parameter('symptomatic_delay_mean', 3., 0., 50.,
                      'mean time from becoming contagious to having symptoms'),
    'sigma': Parameter('symptomatic_delay_sigma', 1., 0.01, 20.,
                       'std dev of times from from becoming contagious to having symptoms')
}
symptomatic_delay = Delay('symptomatic_delay', 'norm', symptomatic_delay_pars,
                          bc_model)

non_quarantined_pop = Population('non_quarantined', 0, 'contagious but not identified')
non_quarantined_delay = Delay('fast', 'fast', model=bc_model)

bc_model.add_connector(
    Splitter('contagious split', contagious_pop, 
             [icu_pop, symptomatic_pop, non_quarantined_pop],
             [to_icu_fraction, symptomatic_fraction],
              [to_icu_delay, symptomatic_delay, non_quarantined_delay]))



# Those symptomatic may get tested, and finally
# with a positive report, become isolated and therefore not contagious.
# Meanwhile there will be others who are contagious who may not follow
# that path, but naturally recover. To keep track of recoveries
# outside of the reporting path, we use a chain.
chain = []


# Most of those with symptoms get tested a few days after onset of symptoms
#oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo

tested_pop = Population('tested', 0,
                        'People who have been tested', color='deepskyblue')
tested_fraction = Parameter('tested_frac', 0.8, 0., 1.,
                            'fraction of symptomatic people who get tested',
                            hidden=False)
tested_delay_pars = {
    'mean': Parameter('tested_delay_mean', 3., 0., 50.,
                      'mean time from having symptoms to getting tested',
                      hidden=False),
    'sigma': Parameter('tested_delay_sigma', 1., 0.01, 20.,
                       'standard deviation of times from from having symptoms to getting tested')
}
tested_delay = Delay('tested_delay', 'norm', tested_delay_pars, bc_model)

chain.append(
    Propagator('symptomatic to tested', symptomatic_pop,
               tested_pop, tested_fraction, tested_delay))

# Once tested, infected symptomatic people must wait for their positive report
#ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo

reported_pop = Population('reported', 0,
                          'Infected people who received a positive test report',
                          hidden=False, color='forestgreen', show_sim=True)
reported_fraction = Parameter('reported_frac', 0.95, 0., 1.,
                              'fraction of tested infected people who will '+\
                                  'receive a positive report')
reported_delay_pars = {
    'mean': Parameter('reported_delay_mean', 5., 0., 50.,
                      'mean time from having test to getting report', hidden=False),
    'sigma': Parameter('reported_delay_sigma', 2., 0.01, 20.,
                       'standard deviation of times from having test to getting report')
}
reported_delay = Delay('reported_delay', 'norm', reported_delay_pars, bc_model)

chain.append(
    Propagator('tested to reported', tested_pop,
               reported_pop, reported_fraction, reported_delay))

# This is the end of the chain. All remainders are not quarantined
#ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo

chain_delay = Delay('fast', 'fast', model=bc_model)
chain_fraction = Parameter('chain_fraction', 1., 0., 1., 'fraction of remainder to non_quarantined')

bc_model.add_connector(
    Chain('reporting chain', symptomatic_pop, non_quarantined_pop,
          chain, chain_fraction, chain_delay, bc_model))

nq_recovered_pop = Population('rec_non_quarantined', 0, 'natural recovery, never isolated')
nq_recover_fraction = Parameter('nq_recover_frac', 1.0, 0., 1.,
                                'fraction of asymptomatic contagious people who recover')
bc_model.add_connector(
    Propagator('non-quarantined to recovered', non_quarantined_pop,
               nq_recovered_pop, nq_recover_fraction, recover_delay))

# Some with positive tests will be admitted to hospital
#oooooooooooooooooooooooooooooooooooooooooooooooooooooo

hospitalized_pop = Population('hospitalized', 0,
                              'Total hospitalization cases', color='slategrey', show_sim=True)
hospitalized_fraction = Parameter('hosp_frac', 0.2, 0., 1.,
                                  'fraction of those with postive tests who will '+\
                                  'be admitted to hospital', hidden=False)
hospitalized_delay_pars = {
    'mean': Parameter('hosp_delay_mean', 3., 0., 50.,
                      'mean time from positive report to hospitalization', hidden=False),
    'sigma': Parameter('hosp_delay_sigma', 2., 0.01, 20.,
                       'standard deviation of times from positive report to hospitalization')
}
hospitalized_delay = Delay('hospitalized_delay', 'norm', hospitalized_delay_pars, bc_model)

bc_model.add_connector(
    Propagator('reported to hospitalized', reported_pop,
               hospitalized_pop, hospitalized_fraction, hospitalized_delay))

# make a copy to keep track of how many remain in hospital
#ooooooooooooooooooooooooooooooooooooooooooooooooooooooooo

in_hospital_pop = Population('in_hospital', 0,
                             'People currently in hospital',
                             hidden=False, color='darkcyan', show_sim=True)
bc_model.add_connector(
    Adder('copy hospitalizations', hospitalized_pop, in_hospital_pop))

# people are released from hospital (ICU tracked separately)
#oooooooooooooooooooooooooooooooooo

released_pop = Population('released', 0,
                          'People released from hospital', color='darkolivegreen')

released_fraction = Parameter('released_frac', 1., 0., 1.,
                              'fraction of those released from hospital == 1', hidden=False)
released_delay_pars = {
    'mean': Parameter('released_delay_mean', 5., 0., 50.,
                      'mean time from hospital admission to release', hidden=False),
    'sigma': Parameter('released_delay_sigma', 2., 0.01, 20.,
                       'standard deviation of times from hospital admission to release')
}
released_delay = Delay('released_delay', 'norm', released_delay_pars, bc_model)

bc_model.add_connector(
    Propagator('hospitalized to released', hospitalized_pop,
               released_pop, released_fraction, released_delay))


# people stay a while in the ICU and then leave
#oooooooooooooooooooooooooooooooooooooooooooooo

depart_icu_pop = Population('departed ICU', 0,
                            'People who left ICU', color='plum')

in_icu_pop = Population('in_icu', 0,
                     'People remaining ICU', hidden=False, color='deeppink', show_sim=True)

bc_model.add_connector(
    Adder('copy icu admissions', icu_pop, in_icu_pop))

in_icu_delay_pars = {
    'mean': Parameter('in_icu_delay_mean', 14., 0., 50.,
                      'mean time from icu admission to departure', hidden=False),
    'sigma': Parameter('in_icu_delay_sigma', 5., 0.01, 20.,
                       'standard deviation of times from icu admission to departure')
}
in_icu_delay = Delay('in_icu_delay', 'norm', in_icu_delay_pars, bc_model)

depart_icu_fraction = Parameter('depart_icu_frac', 1., 0., 1.,
                              'fraction of those released who eventually leave ICU==1')

bc_model.add_connector(
    Propagator('ICU to departure', in_icu_pop,
               depart_icu_pop, depart_icu_fraction, in_icu_delay))

# keep track of how many remain in hospital and ICU
#oooooooooooooooooooooooooooooooooooooooooooooooooo

bc_model.add_connector(
    Subtractor('remove released patients', in_hospital_pop, released_pop))

bc_model.add_connector(
    Subtractor('remove ICU departures from ICU', in_icu_pop, depart_icu_pop))

#bc_model.add_connector(
#    Subtractor('remove ICU departures from hospital', in_hospital_pop, depart_icu_pop))

# adjust other populations as required
#ooooooooooooooooooooooooooooooooooooo

bc_model.add_connector(
    Subtractor('subtract deaths from total', total_pop, deaths_pop))

bc_model.add_connector(
    Subtractor('remove those tested positive '+\
               'from contagious (isolation)',
               contagious_pop, reported_pop))

bc_model.add_connector(
    Subtractor('removed natural recoveries', contagious_pop, nq_recovered_pop))

bc_model.add_connector(
    Subtractor('removed icu admissions', contagious_pop, icu_pop))

bc_model.add_connector(
    Subtractor('subtract infected from susceptible', susceptible_pop, infected_pop))

# transitional aspects of the model
#oooooooooooooooooooooooooooooooooo

trans_rate_time1 = Parameter('trans_rate_time1', 12, 0., 300.,
                            'number of days before 1st transmission rate change',
                            parameter_type='int', hidden=False)

trans_rate_time2 = Parameter('trans_rate_time2', 50, 0., 300.,
                            'number of days before 2nd transmission rate change',
                            parameter_type='int', hidden=False)

trans_rate_before = Parameter('alpha_before', 0.385, 0., 2.,
                             'initial transmission rate', hidden=False)

trans_rate_after1 = Parameter('alpha_after1', 0.062, 0., 2.,
                             'transmission rate after 1st transition', hidden=False)

trans_rate_after2 = Parameter('alpha_after2', 0.062, 0., 2.,
                              'transmission rate after 2nd transition', hidden=False)

bc_model.add_transition(
    Modifier('transition_rate1', 'rel_days', trans_rate_time1, trans_rate,
             trans_rate_before, trans_rate_after1, enabled=True, model=bc_model))

bc_model.add_transition(
    Modifier('transition_rate2', 'rel_days', trans_rate_time2, trans_rate,
             trans_rate_after1, trans_rate_after2, enabled=False, model=bc_model))

traveller_pop = Population('travellers', 0,
                           'Infected travellers')

traveller_1_time = Parameter('traveller_1_time', 14, 0, 100,
                           'number of days since t0 when travellers_1 depart',
                           parameter_type='int', hidden=False)

traveller_1_number = Parameter('traveller_1_number', 10., 0., 50000.,
                             'number of infected travellers_1',
                             hidden=False)

bc_model.add_transition(
    Injector('infected_travellers_1', 'rel_days', traveller_1_time, traveller_pop,
             traveller_1_number, enabled=False, model=bc_model))

traveller_fraction = Parameter('traveller_frac', 1., 0., 1.,
                               'fraction of infected travellers that enter')
traveller_delay_pars = {
    'mean': Parameter('traveller_delay_mean', 7., 0., 50.,
                      'mean travel time for travellers', hidden=False),
    'sigma': Parameter('traveller_delay_sigma', 1., 0.01, 20.,
                       'standard deviation of travel times',
                       hidden=False)
}
traveller_delay = Delay('traveller_delay', 'norm', traveller_delay_pars, bc_model)

bc_model.add_connector(
    Propagator('travellers to infected', traveller_pop, infected_pop,
               traveller_fraction, traveller_delay))

traveller_2_time = Parameter('traveller_2_time', 21, 0, 100,
                           'number of days since t0 when travellers_2 depart',
                           parameter_type='int', hidden=False)

traveller_2_number = Parameter('traveller_2_number', 10., 0., 50000.,
                             'number of infected travellers_2',
                             hidden=False)

bc_model.add_transition(
    Injector('infected_travellers_2', 'rel_days', traveller_2_time, traveller_pop,
             traveller_2_number, enabled=False, model=bc_model))


# define boot parameters
#ooooooooooooooooooooooo

bc_model.boot_setup(contagious_pop, 1,
                    exclusion_populations=[total_pop, susceptible_pop])

# bootstrap the model:
# - start from just 1 contagious person and run until goal is reached
# (in this case, the number of contagious people reaches a set value)
# - histories (beyond current) are removed and
# histories and futures are scaled to correspond to the goal

#goal_cont_0 = initial_contagious_par.get_value()
#contagious_pop.history[-1] = 1

###
### DEBUG
###
#infected_pop.future=[1]
#contagious_pop.history[-1] = 0

#last_cont = -1
#while contagious_pop.history[-1] < goal_cont_0 and \
#    contagious_pop.history[-1] >= last_cont:
#    last_cont = contagious_pop.history[-1]
#    bc_model.evolve_expectations(1)

#scale = goal_cont_0/contagious_pop.history[-1]

#bc_model.reboot(scale, exclusions=[total_pop, susceptible_pop])

#'symptomatic_frac'

#bc_model.parameters['symptomatic_frac'].set_value(0.9)
#bc_model.boot()
#bc_model.evolve_expectations(200)
#print(bc_model.populations['non_quarantined'].history[-1])

#bc_model.parameters['symptomatic_frac'].set_value(0.09)
#bc_model.boot()
#bc_model.evolve_expectations(200)
#print(bc_model.populations['non_quarantined'].history[-1])

#bc_model.boot()
#bc_model.evolve_expectations(200)
#bc_model.generate_data(200)

#recover_delay_pars['mean'].set_value(4.)

i=1
#bc_model.update_lists()
i=2

#with open('model.pickle', 'wb') as f:
#    pickle.dump(bc_model, f, pickle.HIGHEST_PROTOCOL)#

bc_model.save_file('model_v4_0.pypm')
