# -*- coding: utf-8 -*-
"""
Splitter: A Connector class that distributes its incoming
population to n others, either in the next time step or distributed in time.

@author: karlen
"""
from scipy import stats

from pypmca.Connector import Connector
from pypmca.Delay import Delay
from pypmca.Parameter import Parameter


class Splitter(Connector):
    """
    Splitter distributes its incoming population to other populations,
    either in the next time step or distributed in time:
        - connector_name: string, short descriptor
        - from_population: Population object, source of population to be
        propagated
        - to_population: list of two or more Population objects,
        the destination populations.
        - fractions: list of Parameter object with expected
        fraction of population to be propagated to the first to_population,
        the second is the fraction of the remaining to go to the next to_population
        and so on, with the remainder going to the other population
        The fractions must be in the range (0.,1.)
        - delay: Delay object or list of the Delay objects that define how
        propagation is spread over time.
    """

    def __init__(self, connector_name: str, from_population, to_population: list,
                 fractions: list, delay: Delay):
        """Constructor
        """
        super().__init__(connector_name, from_population, to_population)

        if isinstance(from_population, list):
            raise TypeError('Splitter (' + self.name +
                            ') from_population cannot be a list')

        if not isinstance(to_population, list):
            raise TypeError('Splitter (' + self.name +
                            ') to_population must be a list')

        if len(to_population) < 2:
            raise ValueError('Splitter (' + self.name +
                             ') to_population must be a list of length > 1')

        if not isinstance(fractions, list):
            raise TypeError('Splitter (' + self.name +
                            ') fractions must be a list object')

        if len(fractions) != len(to_population) - 1:
            raise ValueError('Splitter (' + self.name +
                             ') fractions length must be len(to_population)-1')

        for fraction in fractions:
            if not isinstance(fraction, Parameter):
                raise TypeError('Splitter (' + self.name +
                                ') fractions list must contain Parameter object')

            if fraction.get_value() < 0.:
                raise ValueError('Splitter (' + self.name +
                                 ') fraction cannot be negative')
            if fraction.get_value() > 1.:
                raise ValueError('Splitter (' + self.name +
                                 ') fraction cannot exceed 1')

        self.fractions = fractions

        if isinstance(delay, list):
            for d in delay:
                if not isinstance(d, Delay):
                    raise TypeError('Splitter (' + self.name +
                                    ') delay list must only have Delay objects')

            if len(delay) != len(to_population):
                raise ValueError('Splitter (' + self.name +
                                 ') delay length must match to_population length' +
                                 ' or be a single (common) delay')

        elif not isinstance(delay, Delay):
            raise TypeError('Splitter (' + self.name +
                            ') delay must be a Delay object')
        self.delay = delay

        for i in range(len(self.fractions)):
            self.parameters['fraction' + str(i)] = self.fractions[i]

        if isinstance(delay, list):
            for i in range(len(delay)):
                name = 'delay_' + str(to_population[i])
                if delay[i].delay_parameters is not None:
                    for key in delay[i].delay_parameters:
                        self.parameters[name + key] = delay[i].delay_parameters[key]
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

            remaining = 1.
            fractions = []
            for fraction in self.fractions:
                frac = remaining * fraction.get_value()
                fractions.append(frac)
                remaining -= frac
            fractions.append(remaining)

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
        if len(self.from_population.future) > 0:
            incoming = self.from_population.future[0]
            if incoming > 0.:

                remaining = 1.
                fractions = []
                for fraction in self.fractions:
                    frac = remaining * fraction.get_value()
                    fractions.append(frac)
                    remaining -= frac
                fractions.append(remaining)

                scales = stats.multinomial.rvs(incoming, fractions)
                if isinstance(self.delay, list):
                    for i in range(len(self.to_population)):
                        to_pop = self.to_population[i]
                        to_pop.update_future_data(scales[i], self.delay[i])
                else:
                    for i in range(len(self.to_population)):
                        to_pop = self.to_population[i]
                        to_pop.update_future_data(scales[i], self.delay)
