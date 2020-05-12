# -*- coding: utf-8 -*-
"""
Subtractor: A Connector class that 'immediately' removes population.

The purpose of this connector is to include the effect of removal from
populations. For example it is used to track those currently in hospital, instead
of just monitoring the total hospitalized.

@author: karlen
"""

from pypmca.Connector import Connector
from pypmca.Population import Population


class Subtractor(Connector):
    """
    Subtractor removes from the from_population, the newcomers to the
    to_population:
        - connector_name: string, short descriptor
        - from_population: Population object which will be reduced
        by the newcomers to the to_population
        - to_population: Population object having received the newcomers
        of the from_population totals. Its future is not affected by
        this connector.
    """

    def __init__(self, connector_name: str, from_population: Population, to_population: Population):
        """Constructor
        """
        super().__init__(connector_name, from_population, to_population)

        if not isinstance(from_population, Population):
            raise TypeError('Subtractor (' + self.name +
                            ') from_population must be a Population object')

        if not isinstance(to_population, Population):
            raise TypeError('Subtractor (' + self.name +
                            ') to_population must be a Population object')

        if hasattr(from_population, 'report_noise') and from_population.report_noise:
            raise TypeError('Subtractor (' + self.name +
                            ') from_population cannot have report_noise set to True')

        # Identify the from_population as one not to derive daily numbers
        self.from_population.monotonic = False

    def update_expectation(self):
        """
        do subtraction
        """
        self.__do_subtraction()

    def update_data(self):
        """
        do subtraction
        """
        self.__do_subtraction()

    def __do_subtraction(self):
        """
        Subtract new population going into to_population from from_population
        """
        if len(self.to_population.future) > 0:
            reduction = -1 * self.to_population.future[0]
            self.from_population.update_future_fast(reduction)
