# -*- coding: utf-8 -*-
"""
Adder: A Connector class that 'immediately' adds population.

The purpose of this connector is to copy incoming population into the
from_population to the to_population.

@author: karlen
"""

from pypmca.Connector import Connector
from pypmca.Population import Population

class Adder(Connector):
    """
    Adder copies the new incoming population into the from_population, 
    to the to_population:
        - connector_name: string, short descriptor
        - from_population: Population object which has newcomers. Its future
        is not affected by this connector.
        - to_population: Population object that receives the newcomers
        entering from_population.
    """

    def __init__(self, connector_name: str, from_population: Population, to_population: Population):
        """Constructor
        """
        super().__init__(connector_name, from_population, to_population)

        if not isinstance(from_population, Population):
            raise TypeError('Adder('+self.name+
                            ') from_population must be a Population object')

        if not isinstance(to_population, Population):
            raise TypeError('Adder ('+self.name+
                            ') to_population must be a Population object')

    def update_expectation(self):
        """
        do addition
        """
        self.__do_addition()

    def update_data(self):
        """
        do subtraction
        """
        self.__do_addition()

    def __do_addition(self):
        """
        Add new population going into from_population into to_population
        """

        newcomers = 0
        if len( self.from_population.future) > 0:
            newcomers = self.from_population.future[0]
        self.to_population.update_future_fast(newcomers)
