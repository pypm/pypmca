# -*- coding: utf-8 -*-
"""
A class: A class for performing mathematical operations on parameters

The operator object returns the result using the .getValue() method,
so that and Operator object can be used in place of a parameter object

@author: karlen
"""

from pypmca.Parameter import Parameter

class Operator:
    """
    - parameters: list of parameters
    - operations: string indicating operation to be performed in sequence. Length should be 1 less than length of list
    """

    def __init__(self, parameters: list, operations: str):

        if not isinstance(parameters, list):
            raise TypeError('Operator parameters must be a list')
        else:
            for parameter in parameters:
                if not isinstance(parameter, Parameter):
                    raise TypeError('Operator parameters must be a list of parameters')

        self.parameters = parameters

        if not isinstance(operations, str):
            raise TypeError('Operator operations must be a string')
        if len(operations) != len(parameters)-1:
            raise ValueError('Operator operations must be one less than number of parameters')

        self.operations = operations

    def get_value(self):
        """
        Return the evaluation of the operations
        """
        buff = str(self.parameters[0].get_value())
        for i in range(len(self.operations)):
            buff += self.operations[i:i+1]
            buff += str(self.parameters[i+1].get_value())

        return eval(buff)