# -*- coding: utf-8 -*-
"""
Parameter: A class for defining parameters of the model:

Parameters with parameter_status = 'variable' identifies parameters that are to
be adjusted within optimization or MCMC exploration methods. For such parameters,
priors must be supplied using parameter_function and prior_parameters.

Parameters with parameter_status = 'fixed' can be adjusted manually in an
interactive exploration of model parameters. 

The automatic transition of a parameter at defined point of the evolution is
done by adding a Modifier object to the Model.

The pre-calculated delay distributions need to be recalculated when its
parameters are changed. Such parameters have must_update = True

@author: karlen
"""

class Parameter:
    """
    Model parameter:
    - parameter_name: short descriptor
    - description: long form description for documentation
    - value: value to be set initially or when reset requested
    - parameter_min, parameter_max, allowed range for int or float parameter
    - parameter_type: 'int', 'float', 'bool'
    - parameter_status: 'fixed', 'variable'
    - prior_function: if (real or integer) and variable, the functional form
      of prior. Allowed: 'norm' or 'uniform'
    - prior_parameters for (list of floats, not parameter objects):
        - for 'norm': dict [mean, sigma]
        - for 'uniform': dict [mean, half_width]
        - for boolean: truth probability
    """

    PARAMETER_TYPES = {'int':int, 'float':float, 'bool':bool}
    PARAMETER_STATUSES = ['fixed', 'variable']
    PRIOR_FUNCTIONS = ['norm', 'uniform']
    PRIOR_PARS = {'norm': ['mean','sigma'], 'uniform': ['mean','half_width']}

    def __init__(self, parameter_name, value,
                 parameter_min = 0, parameter_max= 1, description='', 
                 parameter_type='float', parameter_status='fixed',
                 prior_function=None, prior_parameters=None):

        self.name = str(parameter_name)
        self.must_update = False
        self.parent = None

        if parameter_type not in self.PARAMETER_TYPES:
            buff = '/'.join(self.PARAMETER_TYPES)
            raise ValueError('Parameter ('+self.name+') type "'+
                             str(parameter_type)+
                             '" invalid: not in: '+buff)
        self.parameter_type = parameter_type

        self.set_value(value)
        self.initial_value = value
        self.description = description
        self.parameter_min = parameter_min
        self.parameter_max = parameter_max

        if parameter_status not in self.PARAMETER_STATUSES:
            buff = '/'.join(self.PARAMETER_STATUSES)
            raise ValueError('Parameter ('+self.name+') status "'+
                             str(parameter_type)+
                             '" invalid: not in: '+buff)
        self.__status = parameter_status
        self.__status_changed = False

        self.prior_function = None
        self.prior_parameters = None
        if parameter_status == 'variable':
            self.set_variable(prior_function, prior_parameters)

    def __str__(self):
        return self.name

    def get_min(self):
        return self.parameter_min
    
    def get_max(self):
        return self.parameter_max

    def get_value(self):
        """
        Return current value of parameter

        """
        return self.__value

    def set_value(self, new_value):
        """
        Change the value of the parameter

        Parameters
        ----------
        new_value : TYPE must match parameter_type
            new value of the parameter

        Returns
        -------
        None.

        """
        if not isinstance(new_value, self.PARAMETER_TYPES[self.parameter_type]):
            raise TypeError('Parameter ('+self.name+
                            ') value type ('+type(new_value).__name__+
                            ') does not match parameter_type ('+
                            self.parameter_type+')')
        self.__value = new_value

        # if this parameter is used for delay - do a recalculation of the distribution
        if self.must_update:
            self.parent.update()

    def set_must_update(self, parent):
        """ identify this parameter as one that changes to its value requires 
        that a calculation in the parent needs to be updated
        """
        self.must_update = True
        self.parent = parent

    def get_status(self):
        return self.__status

    def set_fixed(self):
        """
        Change the parameter status from 'variable' to 'fixed'.

        Returns
        -------
        None.

        """
        self.__status = 'fixed'
        self.__status_changed = True

    def set_variable(self, prior_function=None, prior_parameters=None):
        """
        Change the parameter status from 'fixed' to 'variable'. If the
        parameter had previously been 'variable', and the prior is not
        specified, the prior previously used is reapplied by default.

        Parameters
        ----------
        prior_function : TYPE, optional
            DESCRIPTION. if (real or integer) and variable, the functional form
            of prior. Allowed: 'norm' or 'uniform'
        prior_parameters : TYPE, optional
            DESCRIPTION. (list of floats, not parameter objects):
                - for 'norm': dict  having keys ['mean', 'sigma']
                - for 'uniform': dict having keys ['mean', 'half_width']
                - for boolean: float, truth probability

        if prior not supplied, the previously used prior will be used

        Returns
        -------
        None.

        """
        
        ## by default use the previously used prior
        if prior_function is None and prior_parameters is None:
            prior_function = self.prior_function
            prior_parameters = self.prior_parameters

        if self.parameter_type != 'bool':
            if prior_function is None:
                raise TypeError('Variable parameter ('+self.name+
                                ') prior_function cannot be None')
            if prior_function not in self.PRIOR_FUNCTIONS:
                buff = '/'.join(self.PRIOR_FUNCTIONS)
                raise ValueError('Parameter ('+self.name+') prior_function "'+
                                 str(self.parameter_type)+
                                 '" invalid: not in: '+buff)
            self.prior_function = prior_function

            if prior_parameters is None:
                raise TypeError('Variable parameter ('+self.name+
                                ') prior_parameters cannot be None')
            if not isinstance(prior_parameters, dict):
                raise TypeError('Variable parameter ('+self.name+
                                ') prior_parameters must be a dictionary')
            for key in self.PRIOR_PARS[prior_function]:
                if key not in prior_parameters:
                    raise ValueError('Variable parameter ('+self.name+
                                     ') prior_parameters must include '+key)
                if not isinstance(prior_parameters[key], float):
                    raise TypeError('Variable parameter ('+self.name+
                                    ' prior_parameters must be floats')

        elif self.parameter_type == 'bool':
            if self.prior_parameters is None:
                raise TypeError('Variable parameter ('+self.name+
                                ') prior_parameters cannot be None')
            if not isinstance(prior_parameters, float):
                raise TypeError('Variable parameter ('+self.name+
                                ' prior_parameters must be a float')
        self.prior_parameters = prior_parameters

        self.__status = 'variable'
        self.__status_changed = True
