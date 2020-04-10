# -*- coding: utf-8 -*-
"""
Modifier: A Transition class that allows a parameter value to be modified at
a specific point in the evolution

@author:karlen
"""

from Parameter import Parameter
from Transition import Transition

class Modifier(Transition):
    """
    A modifier object will cause a parameter value to change at specified point
    in the evolution.

    The time_step for which the action is taken can be specified as one of 'time_spec':
    - 'rel_days': float, number of days after start (t0)
    - 'rel_steps': int, number of time_steps after start (t0)

    The value of parameter_before is change to the value contained in parameter_after.

    In all cases, the time is specified by a Parameter object

    """


    def __init__(self, transition_name, time_spec, transition_time,
                 parameter_before, parameter_after, model=None):
        """Constructor
        """
        super().__init__(transition_name, time_spec, transition_time, model)

        if parameter_before is None:
            raise TypeError('Error in constructing transition ('+self.name+
                            '): parameter_before cannot be None')

        if parameter_after is None:
            raise TypeError('Error in constructing transition ('+self.name+
                            '): parameter_after cannot be None')

        if not isinstance(parameter_before, Parameter):
            raise TypeError('Error in constructing transition ('+self.name+
                            '): parameter_before must be a Parameter object')

        if not isinstance(parameter_after, Parameter):
            raise TypeError('Error in constructing transition ('+self.name+
                            '): parameter_after must be a Parameter object')

        self.parameter_before = parameter_before
        self.parameter_after = parameter_after

        self.initial_value = parameter_before.get_value()


    def take_action(self, expectations=True):
        """
        Modify the value of the parameter
        """
        self.parameter_before.set_value(self.parameter_after.get_value())

    def reset(self):
        """
        Revert to the oriinal value for the parameter
        """
        self.parameter_before.set_value(self.initial_value)
