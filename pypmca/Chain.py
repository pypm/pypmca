# -*- coding: utf-8 -*-
"""
Chain: A Connector class processes a chain of one-to-one propagators, and
accumulating the remainder for a new population

@author: karlen
"""
from scipy import stats
import numpy as np

from pypmca.Connector import Connector
from pypmca.Delay import Delay
from pypmca.Parameter import Parameter
from pypmca.Population import Population
from pypmca.Propagator import Propagator


class Chain(Connector):
    """
    Chain: processes multiple one-to-one propagators and
    sends the remainder of the population to the to_population
        - connector_name: string, short descriptor
        - from_population: Population object at the start of the chain
        - to_population: Population object that receives the remainder
        the destination populations.
        - chain: list of one-to-one Propagator objects with 'norm' delays
        - fraction: Parameter object with expected
        fraction of remaining population to be propagated to the to_population,
        - delay: Delay object that defines how the remainder
        propagation is spread over time.
    """

    def __init__(self, connector_name: str, from_population: Population, to_population: Population,
                 chain: list, fraction: Parameter, delay: Parameter, model):
        """Constructor
        """
        super().__init__(connector_name, from_population, to_population)

        if not isinstance(from_population, Population):
            raise TypeError('Chain (' + self.name +
                            ') from_population must be a Population object')

        if not isinstance(to_population, Population):
            raise TypeError('Chain (' + self.name +
                            ') to_population must be a Population object')

        if not isinstance(chain, list):
            raise TypeError('Chain (' + self.name +
                            ') chain must be a list of Propagator objects')

        for prop in chain:
            if not isinstance(prop, Propagator):
                raise TypeError('Chain (' + self.name +
                                ') chain must be a list of only Propagator objects')
            if not isinstance(prop.from_population, Population):
                raise ValueError('Chain (' + self.name +
                                 ') chain must be a list of one-to-one Propagator objects')
            if not isinstance(prop.to_population, Population):
                raise ValueError('Chain (' + self.name +
                                 ') chain must be a list of one-to-one Propagator objects')
            if prop.delay.delay_type != 'norm':
                raise ValueError('Chain (' + self.name +
                                 ') chain must be a list of Propagator objects with "norm" delays')

            for key in prop.delay.delay_parameters:
                # notify parameter that this object's update method must be call if changed
                prop.delay.delay_parameters[key].set_must_update(self)

            prop.fraction.set_must_update(self)

        self.chain = chain

        if not isinstance(fraction, Parameter):
            raise TypeError('Chain (' + self.name +
                            ') fraction must be a Parameter object')
        if fraction.get_value() < 0.:
            raise ValueError('Chain (' + self.name +
                             ') fraction cannot be negative')
        if fraction.get_value() > 1.:
            raise ValueError('Chain (' + self.name +
                             ') fraction cannot exceed 1')
        self.fraction = fraction

        if not isinstance(delay, Delay):
            raise TypeError('Chain (' + self.name +
                            ') delay must be a Delay object')
        self.delay = delay

        if model is None:
            raise TypeError('Chain (' + self.name + ') model cannot be None')
        if not hasattr(model, 'get_time_step'):
            raise TypeError('Chain (' + self.name +
                            ') model must be a model object')
        self.model = model

        self.parameters['fraction'] = self.fraction
        if delay.delay_parameters is not None:
            for key in delay.delay_parameters:
                self.parameters['delay_' + key] = delay.delay_parameters[key]

        self.propagators = None
        self.remainder_propagator = None
        self.__setup_compound_propagation()

    # the same issue for delay parameters (update required after a change)
    def update(self):
        self.__setup_compound_propagation()

    def __setup_compound_propagation(self):
        """ For each link in chain, work out its effective combined delay
        and the combined fraction, to define a set of propagators """
        remainder_fraction = 0.
        current_fraction = 1.
        mean = 0.
        sigma = 0.
        self.propagators = []
        for prop in self.chain:
            prop_frac = prop.fraction.get_value()
            remainder_fraction += current_fraction * (1. - prop_frac)
            current_fraction *= prop_frac
            mean += prop.delay.delay_parameters['mean'].get_value()
            prop_sig = prop.delay.delay_parameters['sigma'].get_value()
            sigma2 = sigma * sigma + prop_sig * prop_sig
            sigma = np.sqrt(sigma2)

            frac = Parameter('frac', prop_frac)
            delay_pars = {
                'mean': Parameter('mean', mean, parameter_min=0., parameter_max=2.*mean),
                'sigma': Parameter('sigma', sigma, parameter_min= 0., parameter_max=2.*mean)
            }
            delay = Delay('delay_'+str(prop), 'norm', delay_pars, self.model)
            self.propagators.append(
                Propagator(prop.name + '_chain', self.from_population,
                           prop.to_population, frac, delay))
        frac = Parameter('frac', remainder_fraction * self.fraction.get_value())
        self.remainder_propagator = \
            Propagator(self.name + '_remainder', self.from_population,
                       self.to_population, frac, self.delay)

    def update_expectation(self):
        """
        Calculate contributions to other populations by updating their
        future_expectations
        """
        if len(self.from_population.future) > 0:
            incoming = self.from_population.future[0]

            current_fraction = 1.
            for prop in self.propagators:
                current_fraction *= prop.fraction.get_value()
                scale = current_fraction * incoming
                prop.to_population.update_future_expectation(scale, prop.delay)

            prop = self.remainder_propagator
            scale = prop.fraction.get_value() * incoming
            prop.to_population.update_future_expectation(scale, prop.delay)

    def update_data(self):
        """
        Simulate data for other populations by updating their
        future_expectations
        """
        if len(self.from_population.future) > 0:
            incoming = self.from_population.future[0]

            remaining = 0
            in_chain = incoming
            for prop in self.propagators:
                frac = prop.fraction.get_value()
                selected = stats.binom.rvs(in_chain, frac)
                remaining += (in_chain - selected)
                in_chain = selected
                prop.to_population.update_future_data(selected, prop.delay)

            prop = self.remainder_propagator
            # note that there is no need to apply the remainder fraction here
            selected = stats.binom.rvs(remaining, self.fraction.get_value())
            prop.to_population.update_future_data(selected, prop.delay)
