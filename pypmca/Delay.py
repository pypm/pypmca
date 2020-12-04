# -*- coding: utf-8 -*-
"""
Delay: A class that holds information about delayed propagation

@author: karlen
"""

from scipy import stats

from pypmca.Parameter import Parameter


class Delay:
    """
    A class to define delays in propagation
        - delay_name: string, short descriptor
        - delay_type: string:
            - fast: all passed in next time step
            - norm: delay_parameters are mean and standard deviation
            - uniform: delay_parameters are mean and half width
            - erlang: delay_parameters are mean and 'k' (number of ODE stages)
            - gamma: delay_parameters are mean and standard deviation
        - delay_parameters: dictionary of Parameter objects. The units are days.
            - Can be None if delay_type is 'fast'
        - model: Model object that this delay is to be used with
            - Can be None if delay_type is 'fast', but otherwise must be
              specified to access the model time_step
    """

    DELAY_TYPES = ['fast', 'norm', 'uniform', 'erlang', 'gamma']
    DELAY_PAR_KEYS = {'fast': [], 'norm': ['mean', 'sigma'],
                      'uniform': ['mean', 'half_width'],
                      'erlang': ['mean', 'k'], 'gamma': ['mean', 'sigma']}
    EPSILON = 0.0001  # futures assigned until cdf reaches (1.-EPS)

    def __init__(self, delay_name: str, delay_type: str,
                 delay_parameters: dict = None, model = None):
        """Constructor
        """

        self.name = str(delay_name)
        self.delay_type = None
        self.delay_parameters = None
        self.future_expectations = None

        if delay_type != 'fast':
            if model is None:
                raise TypeError('Delay (' + self.name + ') model cannot be None')
            if not hasattr(model, 'get_time_step'):
                raise TypeError('Delay (' + self.name +
                                ') model must be a model object')
        self.model = model

        self.setup_delay(delay_type, delay_parameters)

    def __str__(self):
        return self.name

    def setup_delay(self, delay_type, delay_parameters):

        if str(delay_type) not in self.DELAY_TYPES:
            buff = '/'.join(self.DELAY_TYPES)
            raise ValueError('Delay (' + self.name + ') type "' + str(delay_type) +
                             '" invalid: not in: ' + buff)

        if delay_type != 'fast':
            if delay_parameters is None:
                raise TypeError('Delay (' + self.name + ') parameters cannot be None')
            if len(delay_parameters) < 2:
                raise ValueError('Delay (' + self.name +
                                 ') parameters dictionary must have length > 1')
            for key in delay_parameters:
                if not isinstance(delay_parameters[key], Parameter):
                    raise TypeError('Delay (' + self.name +
                                    ') parameters must be Parameter objects')
                if key == 'k':
                    if delay_parameters[key].parameter_type != 'int':
                        raise TypeError('Delay (' + self.name +
                                        ') k erlang parameter must be an integer parameter')
                elif delay_parameters[key].parameter_type != 'float':
                    raise TypeError('Delay (' + self.name +
                                    ') parameter ' + key + ' must be a float parameter')
                # flag the parameter that if changed, the parent
                # update method must be called
                delay_parameters[key].set_must_update(self)
            keys = self.DELAY_PAR_KEYS[delay_type]
            for key in keys:
                if key not in delay_parameters:
                    raise ValueError('Delay (' + self.name +
                                     ') parameters missing parameter: ' + key)

        self.delay_type = delay_type
        self.delay_parameters = delay_parameters

        self.__set_future_expectations()

    def update(self):
        self.__set_future_expectations()

    def __set_future_expectations(self):
        """
        initializes the probability for future time steps
        this needs to be redone when delay parameters are changed
        or if the model time_step changes
        """
        if self.delay_type == 'fast':
            self.future_expectations = [1.]

        else:
            time_step = self.model.get_time_step()
            day = 1. / time_step
            mean = self.delay_parameters['mean'].get_value()
            loc = 0.
            scale = 0.
            a = 0
            dist = None
            if self.delay_type == 'norm':
                loc = mean * day
                sigma = self.delay_parameters['sigma'].get_value()
                scale = sigma * day
                dist = stats.norm
            elif self.delay_type == 'gamma':
                loc = 0.
                sigma = self.delay_parameters['sigma'].get_value()
                mu = max(0.001,mean)
                a = (mu/sigma)**2
                scale = sigma**2/mu * day
                dist = stats.gamma
            elif self.delay_type == 'uniform':
                half_width = self.delay_parameters['half_width'].get_value()
                loc = mean * day - half_width * day
                scale = 2. * half_width * day
                dist = stats.uniform
            elif self.delay_type == 'erlang':
                # this mimics an ODE type delay with # of stages = a 
                loc = 0.
                a = self.delay_parameters['k'].get_value()
                scale = mean * day / a
                dist = stats.erlang

            # use cdf with 0.5 time_step offsets
            cdf = 0.
            if self.delay_type in ['gamma','erlang']:
                cdf = dist.cdf(0.5, a, loc=loc, scale=scale)
            else:
                cdf = dist.cdf(0.5, loc=loc, scale=scale)

            future = [cdf]
            step = 0
            while cdf < (1 - self.EPSILON):
                last_cdf = cdf
                step += 1
                if self.delay_type in ['gamma','erlang']:
                    cdf = dist.cdf(step + 0.5, a, loc=loc, scale=scale)
                else:
                    cdf = dist.cdf(step + 0.5, loc=loc, scale=scale)

                increment = cdf - last_cdf
                future.append(increment)
            missing = 1. - cdf
            future[-1] += missing
            self.future_expectations = future
