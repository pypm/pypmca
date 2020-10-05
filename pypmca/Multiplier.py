# -*- coding: utf-8 -*-
"""
Multiplier: A Connector class that propagates a population determined by
the product of the two "from_populations" to the "to_population" either in
the next time step or distributed in time:

@author: karlen
"""

from scipy import stats

from pypmca.Connector import Connector
from pypmca.Delay import Delay
from pypmca.Parameter import Parameter


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
        - distribution: newly infected production distribution. 
            poisson (indendent infections)
            nbinom (negative binomial - in popular use)
            nbinom_par = p in the scipy nbinom. limit 0.001<p<0.999
    """

    def __init__(self, connector_name: str, from_population: list, to_population,
                 scale_parameter: Parameter, delay: Delay, model=None, distribution: str = 'poisson',
                 nbinom_par: Parameter = None):
        """Constructor
        """
        super().__init__(connector_name, from_population, to_population)

        if not isinstance(from_population, list):
            raise TypeError('Multiplier (' + self.name +
                            ') from_population must be a list')

        if len(from_population) not in [2, 3]:
            raise ValueError('Multiplier (' + self.name +
                             ') from_population list length must be 2 or 3')

        if not isinstance(scale_parameter, Parameter):
            raise TypeError('Multiplier (' + self.name +
                            ') scale_parameter must be a Parameter object')
        self.scale_parameter = scale_parameter

        if not isinstance(delay, Delay):
            raise TypeError('Multiplier (' + self.name +
                            ') delay must be a Delay object')
        self.delay = delay

        if model is None:
            raise TypeError('Multiplier (' + self.name + ') model cannot be None')
        if not hasattr(model, 'get_time_step'):
            raise TypeError('Multiplier (' + self.name +
                            ') model must be a model object')
        self.model = model

        self.parameters[str(self.scale_parameter)] = self.scale_parameter
        if delay.delay_parameters is not None:
            for key in delay.delay_parameters:
                self.parameters['delay_' + key] = delay.delay_parameters[key]

        self.__distribution = None
        self.__nbinom_par = None
        self.set_distribution(distribution, nbinom_par)

    def set_distribution(self, distribution, nbinom_par):
        if distribution not in ['poisson', 'nbinom']:
            raise ValueError('Multiplier (' + self.name +
                             ') distribution must be poisson or nbinom')
        self.__distribution = distribution

        if distribution == 'nbinom' and \
                (nbinom_par is None or not isinstance(nbinom_par, Parameter)):
            raise TypeError('Multiplier (' + self.name +
                            ') nbinom_par must be a Parameter object')
        self.__nbinom_par = nbinom_par

        if distribution == 'nbinom':
            self.parameters['nbinom_par'] = self.__nbinom_par
        # in case parameter changed, update the model list of parameters
        self.model.update_lists()

    def get_distribution(self):
        return self.__distribution, self.__nbinom_par

    def update_expectation(self):
        """
        Calculate contributions to other populations by updating their
        future_expectations
        """
        scale = self.__get_scale()
        if scale > 0:
            self.to_population.update_future_expectation(scale, self.delay)

    def update_data(self):
        """
        Simulate data for other populations by updating their
        future_expectations
        """
        scale = self.__get_scale()
        if scale > 0:
            if self.__distribution == 'poisson':
                n = stats.poisson.rvs(scale)
                self.to_population.update_future_data(n, self.delay)
            else:
                p = self.__nbinom_par.get_value()
                if p < 0.001:
                    p = 0.001
                if p > 0.999:
                    p = 0.999
                r = scale * p / (1. - p)
                n = 0
                if r > 0.:
                    n = stats.nbinom.rvs(r, p)
                self.to_population.update_future_data(n, self.delay)

    def __get_scale(self):
        """
        Calculate expected number to be sent to "to_population"
        """
        scale = 0.
        if self.from_population[0].history[-1] > 0 and self.from_population[1].history[-1] > 0:
            denom = 1.
            if len(self.from_population) == 3 and \
                    self.from_population[2].history[-1] > 0:
                denom = self.from_population[2].history[-1]

            ratio = 1. * self.from_population[0].history[-1] / denom
            scale = self.scale_parameter.get_value() * ratio * self.from_population[1].history[-1]
            scale *= self.model.get_time_step()
        else:
            iii=1
        return scale
