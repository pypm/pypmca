# -*- coding: utf-8 -*-
"""
Multiplier: A Connector class that propagates a population determined by
the product of the two "from_populations" to the "to_population" either in
the next time step or distributed in time:

@author: karlen
"""

from scipy import stats

from Connector import Connector
from Delay import Delay
from Parameter import Parameter

class Multiplier(Connector):
    """
    Multiplier send a population determined by the product of
    two "from_populations" to the "to_population",
    either in the next time step or distributed in time:
        - connector_name: string, short descriptor
        - from_population: list of two Population objects to be multiplied
        If a third from_population is included, it enters in the denominator
        - to_population: Population object that is the destination population
        - scale_parameter: Parameter object that multiplies the product
        of the from_population totals
        - delay: Delay object that defines how propagation is spread over time.
        - model: necessary to access the time_step
    """

    def __init__(self, connector_name, from_population, to_population,
                 scale_parameter, delay, model):
        """Constructor
        """
        super().__init__(connector_name, from_population, to_population)

        if not isinstance(from_population, list):
            raise TypeError('Multiplier ('+self.name+
                            ') from_population must be a list')

        if len(from_population) not in [2, 3]:
            raise ValueError('Multiplier ('+self.name+
                             ') from_population list length must be 2 or 3')

        if not isinstance(scale_parameter, Parameter):
            raise TypeError('Multiplier ('+self.name+
                            ') scale_parameter must be a Parameter object')
        self.scale_parameter = scale_parameter

        if not isinstance(delay, Delay):
            raise TypeError('Multiplier ('+self.name+
                            ') delay must be a Delay object')
        self.delay = delay
        
        if model is None:
            raise TypeError('Multiplier ('+self.name+') model cannot be None')
        if not hasattr(model, 'get_time_step'):
            raise TypeError('Multiplier ('+self.name+
                            ') model must be a model object')
        self.model = model
        
        self.parameters[str(self.scale_parameter)] = self.scale_parameter
        if delay.delay_parameters is not None:
            for key in delay.delay_parameters:
                self.parameters['delay_'+key] = delay.delay_parameters[key]

    def update_expectation(self):
        """
        Calculate contributions to other populations by updating their
        future_expectations
        """
        scale = self.__get_scale()
        self.to_population.update_future_expectation(scale, self.delay)

    def update_data(self):
        """
        Simulate data for other populations by updating their
        future_expectations
        """
        scale = self.__get_scale()
        n = stats.poisson.rvs(scale)
        self.to_population.update_future_data(n, self.delay)

    def __get_scale(self):
        """
        Calculate expected number to be sent to "to_population"
        """
        denom = 1.
        if len(self.from_population) == 3 and \
            self.from_population[2].history[-1] > 0:
            denom = self.from_population[2].history[-1]

        ratio = 1. * self.from_population[0].history[-1] / denom
        scale = self.scale_parameter.get_value() * ratio * self.from_population[1].history[-1]
        scale *= self.model.get_time_step()
        return scale
