# -*- coding: utf-8 -*-
"""
Transition: A base'virtual' class for defining transitions during the evolution

The derived classes must define the action methods


@author:karlen
"""

from pypmca.Parameter import Parameter


class Transition:
    """
    Base class for defining transitions during the evolution

    The time_step for which the action is taken can be specified as one of 'time_spec':
    - 'rel_days': float, number of days after start (t0)
    - 'rel_steps': int, number of time_steps after start (t0)
    
    The transition can be disabled by setting enabled=False
    
    enabled: specifies default state for the transition

    In all cases, the time is specified by a Parameter object
    """

    TIME_SPECS = {'rel_days': 'int', 'rel_steps': 'int'}

    def __init__(self, transition_name, time_spec, transition_time,
                 enabled=True, model=None):
        """Constructor
        """
        self.name = str(transition_name)

        if not isinstance(time_spec, str):
            raise TypeError('Error in constructing transition (' + self.name +
                            '): time_spec argument must be a str')

        if time_spec not in self.TIME_SPECS:
            specs = []
            for key in self.TIME_SPECS:
                specs.append(key)
            buff = '/'.join(specs)
            raise TypeError('Error in constructing transition (' + self.name +
                            '): time_spec argument must be one of' + buff)
        self.time_spec = time_spec

        if not isinstance(transition_time, Parameter):
            raise TypeError('Error in constructing transition (' + self.name +
                            '): transition_time argument must be a Parameter object')

        if transition_time.parameter_type != self.TIME_SPECS[self.time_spec]:
            raise TypeError('Transition (' + self.name +
                            ') transition time (' +
                            transition_time.parameter_type +
                            ') does not match time_spec (' +
                            self.time_spec + ')')

        if transition_time.get_value() <= 0:
            raise ValueError('Transition (' + self.name +
                             ') transition time must be larger than zero.')

        self.transition_time = transition_time
        self.transition_time.set_must_update(self)

        if model is None:
            raise TypeError('Transition (' + self.name + ') model cannot be None')
        if not hasattr(model, 'get_time_step'):
            raise TypeError('Transition (' + self.name +
                            ') model must be a model object')
        self.model = model

        self.trigger_step = None
        self.__calculate_trigger_step()

        self.parameters = {str(self.transition_time): self.transition_time}

        self.enabled = enabled

    def update(self):
        """ Method called if transition_time parameter changed
        """
        self.__calculate_trigger_step()

    def __calculate_trigger_step(self):
        """ calculate the step when this transition should take place
        """
        time_step = self.model.get_time_step()
        day = 1. / time_step

        if self.time_spec == 'rel_days':
            self.trigger_step = int(round(self.transition_time.get_value() * day))
        elif self.time_spec == 'rel_steps':
            self.trigger_step = self.transition_time.get_value()

    def __str__(self):
        return self.name

    def take_action(self, expectations=True):
        """
        This method must be defined in derived classes

        Raises
        ------
        NotImplementedError if not defined in the derived class

        Returns
        -------
        None.

        """
        raise NotImplementedError()

    def reset(self):
        """
        This method must be defined in derived classes

        Raises
        ------
        NotImplementedError if not defined in the derived class

        Returns
        -------
        None.

        """
        raise NotImplementedError()
