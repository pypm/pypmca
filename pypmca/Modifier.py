# -*- coding: utf-8 -*-
"""
Modifier: A Transition class that allows a parameter value to be modified at
a specific point in the evolution

@author:karlen
"""

from pypmca.Parameter import Parameter
from pypmca.Transition import Transition


class Modifier(Transition):
    """
    A modifier object will cause a parameter value to change at specified point
    in the evolution.

    The time_step for which the action is taken can be specified as one of 'time_spec':
    - 'rel_days': int, number of days after start (t0)
    - 'rel_steps': int, number of time_steps after start (t0)
    
    After the time_step is reached in the model, the value for parameter is
    changed to the value contained in parameter_after.
    
    The _before value of the parameter is specified so that
    the transition can be reversed. To restart the model the transitions are
    reset in reverse order, so that the _before value for the first transition
    is used at the beginning.

    The transition time is specified by a Parameter object.

    Transitions are not applied during the model "boot"

    """

    def __init__(self, transition_name: str, time_spec: str, transition_time: Parameter, parameter: Parameter,
                 parameter_before: Parameter, parameter_after: Parameter, enabled: bool = True, model=None):
        """Constructor
        """
        super().__init__(transition_name, time_spec, transition_time, enabled, model)

        if parameter is None:
            raise TypeError('Error in constructing transition (' + self.name +
                            '): parameter_before cannot be None')

        if parameter_before is None:
            raise TypeError('Error in constructing transition (' + self.name +
                            '): parameter_before cannot be None')

        if parameter_after is None:
            raise TypeError('Error in constructing transition (' + self.name +
                            '): parameter_after cannot be None')

        if not isinstance(parameter, Parameter):
            raise TypeError('Error in constructing transition (' + self.name +
                            '): parameter must be a Parameter object')

        if not isinstance(parameter_before, Parameter):
            raise TypeError('Error in constructing transition (' + self.name +
                            '): parameter_before must be a Parameter object')

        if not isinstance(parameter_after, Parameter):
            raise TypeError('Error in constructing transition (' + self.name +
                            '): parameter_after must be a Parameter object')

        self.parameter = parameter
        self.parameter_before = parameter_before
        self.parameter_after = parameter_after

        self.initial_value = parameter_before.get_value()

        self.parameters[str(parameter_after)] = parameter_after
        self.parameters[str(parameter_before)] = parameter_before

    def take_action(self, expectations=True):
        """
        Modify the value of the parameter
        """
        self.parameter.set_value(self.parameter_after.get_value())

    def reset(self):
        """
        Revert to the value of the before parameter
        """
        self.parameter.set_value(self.parameter_before.get_value())
