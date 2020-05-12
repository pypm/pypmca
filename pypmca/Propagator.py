# -*- coding: utf-8 -*-
"""
Propagator: A Connector class that propagates a fraction of its incoming
population to others, either in the next time step or distributed in time.

Note that propagation is done independently to the other populations.
To divide the incoming population into different sub-populations (dependently),
use the Split class instead.

There is no benefit in having a list of to_populations, since it is equivalent
to having two independent propagators.

@author: karlen
"""
from scipy import stats

from pypmca.Connector import Connector
from pypmca.Delay import Delay
from pypmca.Parameter import Parameter


class Propagator(Connector):
    """
    Propagator sends a fraction of an incoming population to other populations,
    either in the next time step or distributed in time:
        - connector_name: string, short descriptor
        - from_population: Population object, source of population to be
        propagated
        - to_population: Population object or list of Population objects,
        destination population(s).
        - fraction: Parameter object or list of Parameter objects, expected
        fraction of population to be propagated. If list provided,
        its length must match to_population, and the fractions are applied
        independently. The sum of all fractions can exceed 1.
        - delay: Delay object or list of Delay objects that define how
        propagation is spread over time. If list provided, its length
        must match length of to_population.
    """

    def __init__(self, connector_name: str, from_population, to_population,
                 fraction, delay: Delay):
        """Constructor
        """
        super().__init__(connector_name, from_population, to_population)

        if isinstance(from_population, list):
            raise TypeError('Propagator (' + self.name +
                            ') from_population cannot be a list')

        if isinstance(fraction, list):
            if not isinstance(to_population, list):
                raise TypeError('Propagator (' + self.name +
                                ') fraction cannot be a list when' +
                                ' to_population is not a list')
            if len(fraction) != len(to_population):
                raise ValueError('Propagator (' + self.name +
                                 ') fraction list length does not match' +
                                 ' to_population list length')
            for frac in fraction:
                if not isinstance(frac, Parameter):
                    raise TypeError('Propagator (' + self.name +
                                    ') fraction list must contain Parameter objects')
        elif not isinstance(fraction, Parameter):
            raise TypeError('Propagator (' + self.name +
                            ') fraction must be a Parameter object')
        self.fraction = fraction

        if not isinstance(delay, Delay):
            raise TypeError('Propagator (' + self.name +
                            ') delay must be a Delay object')
        self.delay = delay

        if isinstance(fraction, list):
            for i in range(len(fraction)):
                name = 'fraction_' + str(to_population[i])
                self.parameters[name] = self.fraction[i]
        else:
            self.parameters['fraction'] = self.fraction
        if isinstance(delay, list):
            for i in range(len(delay)):
                name = 'delay_' + str(to_population[i])
                if delay[i].delay_parameters is not None:
                    for key in delay[i].delay_parameters:
                        self.parameters[name + '_' + key] = delay[i].delay_parameters[key]
        else:
            if delay.delay_parameters is not None:
                for key in delay.delay_parameters:
                    self.parameters['delay_' + key] = delay.delay_parameters[key]

    def update_expectation(self):
        """
        Calculate contributions to other populations by updating their
        future_expectations
        """
        if len(self.from_population.future) > 0:
            incoming = self.from_population.future[0]
            if isinstance(self.to_population, list):
                if isinstance(self.fraction, list):
                    if isinstance(self.delay, list):
                        for i in range(len(self.to_population)):
                            to_pop = self.to_population[i]
                            scale = self.fraction[i].get_value() * incoming
                            to_pop.update_future_expectation(scale, self.delay[i])
                    else:
                        for i in range(len(self.to_population)):
                            to_pop = self.to_population[i]
                            scale = self.fraction[i].get_value() * incoming
                            to_pop.update_future_expectation(scale, self.delay)
                else:
                    scale = self.fraction.get_value() * incoming
                    if isinstance(self.delay, list):
                        for i in range(len(self.to_population)):
                            to_pop = self.to_population[i]
                            to_pop.update_future_expectation(scale, self.delay[i])
                    else:
                        for to_pop in self.to_population:
                            to_pop.update_future_expectation(scale, self.delay)
            else:
                scale = self.fraction.get_value() * incoming
                self.to_population.update_future_expectation(scale, self.delay)

    def update_data(self):
        """
        Simulate data for other populations by updating their
        future_expectations
        """
        if len(self.from_population.future) > 0:
            incoming = self.from_population.future[0]
            if isinstance(self.to_population, list):
                if isinstance(self.fraction, list):
                    if isinstance(self.delay, list):
                        for i in range(len(self.to_population)):
                            to_pop = self.to_population[i]
                            scale = stats.binom.rvs(incoming, self.fraction[i].get_value())
                            to_pop.update_future_data(scale, self.delay[i])
                    else:
                        for i in range(len(self.to_population)):
                            to_pop = self.to_population[i]
                            scale = stats.binom.rvs(incoming, self.fraction[i].get_value())
                            to_pop.update_future_data(scale, self.delay)
                else:
                    if isinstance(self.delay, list):
                        for i in range(len(self.to_population)):
                            to_pop = self.to_population[i]
                            scale = stats.binom.rvs(incoming, self.fraction.get_value())
                            to_pop.update_future_data(scale, self.delay[i])
                    else:
                        for to_pop in self.to_population:
                            scale = stats.binom.rvs(incoming, self.fraction.get_value())
                            to_pop.update_future_data(scale, self.delay)
            else:
                scale = stats.binom.rvs(incoming, self.fraction.get_value())
                self.to_population.update_future_data(scale, self.delay)
