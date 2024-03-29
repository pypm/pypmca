# -*- coding: utf-8 -*-
"""
Adder: A Connector class that 'immediately' combines changes to several independent populations and applies to
a population that represents the combined population

The connector adds all the incoming populations in the from_population list and sends to the to_population.

For the collector population to represent the sum, the initial size of the population must be set to
be the sum of the initial sizes of the other populations.

@author: karlen
"""

from pypmca.Connector import Connector
from pypmca.Population import Population


class Collector(Connector):
    """
    Collector copies the changes in the populations in the list from_population to the to_population:
        - connector_name: string, short descriptor
        - from_population: List of Population objects which has newcomers
        - to_population: Population object that receives changes
    """

    def __init__(self, connector_name: str, from_population: list, to_population: Population):
        """Constructor
        """
        super().__init__(connector_name, from_population, to_population)

        if not isinstance(from_population, list):
            raise TypeError('Collector(' + self.name +
                            ') from_population must be a list')

        # Identify the to_population as one not to derive daily numbers if any of the from_populations
        # have a subtraction
        monotonic = True
        for from_pop in from_population:
            if not from_pop.monotonic:
                monotonic = False
        self.to_population.monotonic = monotonic

    def update_expectation(self):
        """
        add all immediate futures
        """
        newcomers = 0

        for from_pop in self.from_population:
            if len(from_pop.future) > 0:
                # all populations are forced to remain positive
                # do not decrease by more than existing population size
                incoming = from_pop.future[0]
                if incoming < 0:
                    if incoming < -1*from_pop.history[-1]:
                        incoming = -1*from_pop.history[-1]
                newcomers += incoming

        self.to_population.update_future_fast(newcomers)

    def update_data(self):
        """
        same treatment as for expectation
        """

        self.update_expectation()
