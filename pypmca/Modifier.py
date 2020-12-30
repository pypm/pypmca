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

    If linear = True - the parameter is modified in a linear fashion: parameter_after defines the slope (change per step)
    If n_step is not None and has a value >0, the modification stops after n_steps

    Transitions are not applied during the model "boot"

    """

    def __init__(self, transition_name: str, time_spec: str, transition_time: Parameter, parameter: Parameter,
                 parameter_before: Parameter, parameter_after: Parameter, enabled: bool = True, model=None,
                 linear: bool=False, n_step = None):
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

        self.linear = linear
        if n_step is not None:
            if not isinstance(n_step, Parameter):
                raise TypeError('Error in constructing transition (' + self.name +
                                '): n_step must be a Parameter object')
        self.n_step = n_step

    def take_action(self, expectations=True):
        """
        Modify the value of the parameter
        """
        linear = getattr(self,'linear',False)

        if linear:
            current_value = self.parameter.get_value()
            slope = (self.parameter_after.get_value() - self.parameter_before.get_value())/self.n_step.get_value()
            new_value = current_value+slope
            if new_value > self.parameter.get_max():
                new_value = self.parameter.get_max()
            if new_value < self.parameter.get_min():
                new_value = self.parameter.get_min()
            self.parameter.set_value(new_value)
        else:
            self.parameter.set_value(self.parameter_after.get_value())

    def reset(self):
        """
        Revert to the value of the before parameter
        """
        self.parameter.set_value(self.parameter_before.get_value())
