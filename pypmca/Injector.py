# -*- coding: utf-8 -*-
"""
Injector: A Transition class that allows a injection of new members to a populations at
a specific point in the evolution

@author:karlen
"""

from pypmca.Parameter import Parameter
from pypmca.Population import Population
from pypmca.Transition import Transition


class Injector(Transition):
    """
    An injection object will cause the number in a population to change at specified point
    in the evolution.

    The time_step for which the action is taken can be specified as one of 'time_spec':
    - 'rel_days': int, number of days after start (t0)
    - 'rel_steps': int, number of time_steps after start (t0)

    - to_population: will have its numbers increased by an amount defined by
    - injection: a float type Parameter object

    In all cases, the time is specified by a Parameter object

    """

    def __init__(self, transition_name: str, time_spec: str, transition_time: Parameter,
                 to_population: Population, injection: Parameter, enabled: bool = True, model=None):
        """Constructor
        """
        super().__init__(transition_name, time_spec, transition_time, enabled, model)

        if not isinstance(to_population, Population):
            raise TypeError('Injector (' + self.name +
                            ') to_population must be a Population object')
        self.to_population = to_population

        if injection is None:
            raise TypeError('Error in constructing transition (' + self.name +
                            '): injection cannot be None')

        if not isinstance(injection, Parameter):
            raise TypeError('Error in constructing transition (' + self.name +
                            '): injection must be a Parameter object')

        if injection.parameter_type != 'float':
            raise TypeError('Error in constructing transition (' + self.name +
                            '): injection must be a float type Parameter object')
        self.injection = injection

        self.parameters[str(injection)] = injection

    def take_action(self, expectations=True):
        """
        Inject the population
        """
        value = self.injection.get_value()
        if expectations:
            self.to_population.update_future_fast(value)
        else:
            self.to_population.update_future_fast(int(round(value)))

    def reset(self):
        """
        Nothing to reset
        """
