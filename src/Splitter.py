# -*- coding: utf-8 -*-
"""
Splitter: A Connector class that distributes its incoming
population to two others, either in the next time step or distributed in time.

@author: karlen
"""
from scipy import stats

from Connector import Connector
from Delay import Delay
from Parameter import Parameter

class Splitter(Connector):
    """
    Splitter distributes its incoming population to other populations,
    either in the next time step or distributed in time:
        - connector_name: string, short descriptor
        - from_population: Population object, source of population to be
        propagated
        - to_population: list of two Population objects,
        the destination populations.
        - fraction: Parameter object with expected
        fraction of population to be propagated to the first to_population,
        with the remainder going to the other population
        The fractions must be in the range (0.,1.)
        - delay: Delay object or list of two Delay objects that define how
        propagation is spread over time.
    """

    def __init__(self, connector_name, from_population, to_population,
                 fraction, delay):
        """Constructor
        """
        super().__init__(connector_name, from_population, to_population)

        if isinstance(from_population, list):
            raise TypeError('Splitter ('+self.name+
                            ') from_population cannot be a list')

        if not isinstance(to_population, list):
            raise TypeError('Splitter ('+self.name+
                            ') to_population must be a list')

        if len(to_population) != 2:
            raise ValueError('Splitter ('+self.name+
                             ') to_population must be a list of length 2')

        if not isinstance(fraction, Parameter):
            raise TypeError('Splitter ('+self.name+
                            ') fraction must be a Parameter object')

        if fraction.get_value() < 0.:
            raise ValueError('Splitter ('+self.name+
                             ') fraction cannot be negative')
        if fraction.get_value() > 1.:
            raise ValueError('Splitter ('+self.name+
                             ') fraction cannot exceed 1')

        self.fraction = fraction

        if isinstance(delay, list):
            for d in delay:
                if not isinstance(d, Delay):
                    raise TypeError('Splitter ('+self.name+
                                    ') delay list must only have Delay objects')
        elif not isinstance(delay, Delay):
            raise TypeError('Splitter ('+self.name+
                            ') delay must be a Delay object')
        self.delay = delay
        
        self.parameters['fraction'] = self.fraction
        if isinstance(delay, list):
            for i in range(len(delay)):
                name = 'delay_'+str(to_population[i])
                if delay[i].delay_parameters is not None:
                    for key in delay[i].delay_parameters:
                        self.parameters[name+key] = delay[i].delay_parameters[key]
        else:
            if delay.delay_parameters is not None:
                for key in delay.delay_parameters:
                    self.parameters['delay_'+key] = delay.delay_parameters[key]

    def update_expectation(self):
        """
        Calculate contributions to other populations by updating their
        future_expectations
        """
        incoming = self.from_population.future[0]
        fractions = [self.fraction.get_value(), 1. - self.fraction.get_value()]
        if isinstance(self.delay, list):
            for i in range(len(self.to_population)):
                to_pop = self.to_population[i]
                scale = fractions[i] * incoming
                to_pop.update_future_expectation(scale, self.delay[i])
        else:
            for i in range(len(self.to_population)):
                to_pop = self.to_population[i]
                scale = fractions[i] * incoming
                to_pop.update_future_expectation(scale, self.delay)

    def update_data(self):
        """
        Simulate data for other populations by updating their
        future_expectations
        """
        incoming = self.from_population.future[0]
        number = stats.binom.rvs(incoming, self.fraction.get_value())
        remainder = incoming - number
        scales = [number, remainder]
        if isinstance(self.delay, list):
            for i in range(len(self.to_population)):
                to_pop = self.to_population[i]
                to_pop.update_future_data(scales[i], self.delay[i])
        else:
            for i in range(len(self.to_population)):
                to_pop = self.to_population[i]
                to_pop.update_future_data(scales[i], self.delay)
