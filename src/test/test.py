# -*- coding: utf-8 -*-
"""
Created on Mon Apr  6 14:29:36 2020

@author: karlen
"""
import pickle

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

bc_model = Model('BC v1')
bc_model.set_t0(2020, 3, 1)

# Initialization

total_pop = Population('total', 5000000,
                       'total population of the province')
susceptible_pop = Population('susceptible', total_pop.history[0],
                             'number of people who could become infected')
infected_pop = Population('infected', 0,
                          'total number of people ever infected')

# Define the infection cycle
#oooooooooooooooooooooooooooo

initial_contagious_par = Parameter('cont_0', 55., 'Number of contagious people at t0')

contagious_pop = Population('contagious', initial_contagious_par,
                            'number of people that can cause someone to become infected')

trans_rate = Parameter('alpha', 0.385,
                       'mean number of people that a contagious person infects per day')
infection_delay = Delay('fast', 'fast', model=bc_model)

bc_model.add_connector(
    Multiplier('infection_cycle', [susceptible_pop, contagious_pop, total_pop],
               infected_pop, trans_rate, infection_delay, bc_model))

contagious_frac = Parameter('cont_frac', 0.9,
                            'fraction of infected people that become contagious')
contagious_delay_pars = {
    'mean': Parameter('cont_delay_mean', 3.,
                      'mean time from being infected to becoming contagious'),
    'sigma': Parameter('cont_delay_sigma', 1.,
                       'standard deviation of times from being infected to becoming contagious')
}

contagious_delay = Delay('cont_delay', 'norm', contagious_delay_pars, bc_model)

bc_model.add_connector(
    Propagator('infect_to_contagious', infected_pop,
               contagious_pop, contagious_frac, contagious_delay))

# The infected either recover or die
#ooooooooooooooooooooooooooooooooooo

recovered_pop = Population('recovered', 0,
                           'People who have recovered from the illness '+
                           'and are therefore no longer susceptible')
deaths_pop = Population('deaths', 0,
                        'people who have died from the illness')
recover_fraction = Parameter('recover_frac', 0.99,
                             'fraction of infected people who recover')
recover_delay_pars = {
    'mean': Parameter('recover_delay_mean', 14., 'mean time from infection to recovery'),
    'sigma': Parameter('recover_delay_sigma', 5.,
                       'standard deviation of times from infection to recovery')
}
recover_delay = Delay('recover_delay', 'norm', recover_delay_pars, bc_model)

death_delay_pars = {
    'mean': Parameter('death_delay_mean', 14., 'mean time from infection to death'),
    'sigma': Parameter('death_delay_sigma', 5.,
                       'standard deviation of times from infection to death')
}
death_delay = Delay('death_delay', 'norm', death_delay_pars, bc_model)

bc_model.add_connector(
    Splitter('recovery', infected_pop, [recovered_pop, deaths_pop],
             recover_fraction, [recover_delay, death_delay]))

# Those contagious may become symptomatic, then get tested, and finally
# with a positive report, become isolated and therefore not contagious.
# Meanwhile there will be others who are contagious who may not follow
# that path, but naturally recover. To keep track of recoveries
# outside of the reporting path, we use a chain.
chain = []

# Most of those contagious become symptomatic
#oooooooooooooooooooooooooooooooooooooooooooo

symptomatic_pop = Population('symptomatic', 0,
                             'People who have shown symptoms')
symptomatic_fraction = Parameter('symptomatic_frac', 0.9,
                                 'fraction of contagious people who become symptomatic')
symptomatic_delay_pars = {
    'mean': Parameter('symptomatic_delay_mean', 3.,
                      'mean time from becoming contagious to having symptoms'),
    'sigma': Parameter('symptomatic_delay_sigma', 1.,
                       'std dev of times from from becoming contagious to having symptoms')
}
symptomatic_delay = Delay('symptomatic_delay', 'norm', symptomatic_delay_pars,
                          bc_model)

chain.append(
    Propagator('contagious to symptomatic', contagious_pop,
               symptomatic_pop, symptomatic_fraction, symptomatic_delay))

# Most of those with symptoms get tested a few days after onset of symptoms
#oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo

tested_pop = Population('tested', 0,
                        'People who have been tested')
tested_fraction = Parameter('tested_frac', 0.8,
                            'fraction of symptomatic people who become tested')
tested_delay_pars = {
    'mean': Parameter('tested_delay_mean', 3.,
                      'mean time from having symptoms to getting tested'),
    'sigma': Parameter('tested_delay_sigma', 1.,
                       'standard deviation of times from from having symptoms to getting tested')
}
tested_delay = Delay('tested_delay', 'norm', tested_delay_pars, bc_model)

chain.append(
    Propagator('symptomatic to tested', symptomatic_pop,
               tested_pop, tested_fraction, tested_delay))

# Once tested, infected symptomatic people must wait for their positive report
#ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo

reported_pop = Population('reported', 0,
                          'Infected people who received a positive test report')
reported_fraction = Parameter('reported_frac', 0.95,
                              'fraction of tested infected people who will '+\
                                  'receive a positive report')
reported_delay_pars = {
    'mean': Parameter('reported_delay_mean', 5., 'mean time from having test to getting report'),
    'sigma': Parameter('reported_delay_sigma', 2.,
                       'standard deviation of times from having test to getting report')
}
reported_delay = Delay('reported_delay', 'norm', reported_delay_pars, bc_model)

chain.append(
    Propagator('tested to reported', tested_pop,
               reported_pop, reported_fraction, reported_delay))

# This is the end of the chain. All remainders are not quarantined
#ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo

non_quarantined_pop = Population('non_quarantined', 0, 'contagious but not identified')
chain_delay = Delay('fast', 'fast', model=bc_model)
chain_fraction = Parameter('chain_fraction', 1.)

bc_model.add_connector(
    Chain('reporting chain', contagious_pop, non_quarantined_pop,
          chain, chain_fraction, chain_delay, bc_model))

nq_recovered_pop = Population('non_quarantined_recovered', 0, 'natural recovery, never isolated')
nq_recover_fraction = Parameter('recover_frac', 1.0,
                             'fraction of asymptomatic contagious people who recover')
bc_model.add_connector(
    Propagator('non-quarantined to recovered', non_quarantined_pop,
               nq_recovered_pop, nq_recover_fraction, recover_delay))

# Some with positive tests will be admitted to hospital
#oooooooooooooooooooooooooooooooooooooooooooooooooooooo

hospitalized_pop = Population('hospitalized', 0,
                              'Total hospitalization cases')
hospitalized_fraction = Parameter('hospitalized_frac', 0.2,
                                  'fraction of those with postive tests who will '+\
                                  'be admitted to hospital')
hospitalized_delay_pars = {
    'mean': Parameter('hospitalized_delay_mean', 3.,
                      'mean time from positive report to hospitalization'),
    'sigma': Parameter('hospitalized_delay_sigma', 2.,
                       'standard deviation of times from positive report to hospitalization')
}
hospitalized_delay = Delay('hospitalized_delay', 'norm', hospitalized_delay_pars, bc_model)

bc_model.add_connector(
    Propagator('reported to hospitalized', reported_pop,
               hospitalized_pop, hospitalized_fraction, hospitalized_delay))

# make a copy to keep track of how many remain in hospital
#ooooooooooooooooooooooooooooooooooooooooooooooooooooooooo

in_hospital_pop = Population('in_hospital', 0,
                             'People currently in hospital')
bc_model.add_connector(
    Adder('copy hospitalizations', hospitalized_pop, in_hospital_pop))

# people are released from hospital
#oooooooooooooooooooooooooooooooooo

released_pop = Population('released', 0,
                          'People released from hospital')
released_fraction = Parameter('released_frac', 1.,
                              'fraction of those eventually released from hospital')
released_delay_pars = {
    'mean': Parameter('released_delay_mean', 14.,
                      'mean time from hospital admission to release'),
    'sigma': Parameter('released_delay_sigma', 5.,
                       'standard deviation of times from hospital admission to release')
}
released_delay = Delay('released_delay', 'norm', released_delay_pars, bc_model)

bc_model.add_connector(
    Propagator('hospitalized to released', hospitalized_pop,
               released_pop, released_fraction, released_delay))

# keep track of how many remain in hospital
#oooooooooooooooooooooooooooooooooooooooooo

bc_model.add_connector(
    Subtractor('remove released patients', in_hospital_pop, released_pop))

# adjust other populations as required
#ooooooooooooooooooooooooooooooooooooo

bc_model.add_connector(
    Subtractor('subtract deaths from total', total_pop, deaths_pop))

bc_model.add_connector(
    Subtractor('remove those who get a positive report '+\
               'from the contagious population (isolation)',
               contagious_pop, reported_pop))

bc_model.add_connector(
    Subtractor('removed natural recoveries', contagious_pop, nq_recovered_pop))

bc_model.add_connector(
    Subtractor('subtract infected from susceptible', susceptible_pop, infected_pop))

# transitional aspects of the model
#oooooooooooooooooooooooooooooooooo

trans_rate_time = Parameter('trans_rate_time', 12.,
                            'number of days before transmission rate changes')
trans_rate_after = Parameter('trans_rate_after', 0.062,
                             'transmission rate after lockdown')

bc_model.add_transition(
    Modifier('transition_rate', 'rel_days', trans_rate_time, trans_rate,
             trans_rate_after, bc_model))

traveller_pop = Population('traveller_pop', 0,
                           'Infected travellers returning home')

traveller_time = Parameter('traveller_time', 20.,
                           'number of days before travellers start to return')

traveller_number = Parameter('traveller_number', 50.,
                             'number of infected travellers returning')

bc_model.add_transition(
    Injector('infected_travellers', 'rel_days', traveller_time, traveller_pop,
             traveller_number, bc_model))

traveller_fraction = Parameter('traveller_frac', 1.,
                              'fraction of infected travellers that enter')
traveller_delay_pars = {
    'mean': Parameter('traveller_delay_mean', 14.,
                      'mean time for all travellers to enter'),
    'sigma': Parameter('released_delay_sigma', 7.,
                       'standard deviation of times of travellers entering')
}
traveller_delay = Delay('traveller_delay', 'norm', traveller_delay_pars, bc_model)

bc_model.add_connector(
    Propagator('traveller_con', traveller_pop, infected_pop,
               traveller_fraction, traveller_delay))

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




bc_model.evolve_expectations(200)
#bc_model.generate_data(200)

#recover_delay_pars['mean'].set_value(4.)

i=1

#with open('model.pickle', 'wb') as f:
#    pickle.dump(bc_model, f, pickle.HIGHEST_PROTOCOL)
