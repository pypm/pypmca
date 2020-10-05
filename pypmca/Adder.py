# -*- coding: utf-8 -*-
"""
Adder: A Connector class that 'immediately' adds population.

The purpose of this connector is to copy incoming population into the
from_population to the to_population.

A scale factor can be applied. The scale factor should be positive.

@author: karlen
"""
from scipy import stats

from pypmca.Connector import Connector
from pypmca.Population import Population
from pypmca.Parameter import Parameter

class Adder(Connector):
    """
    Adder copies the new incoming population into the from_population, 
    to the to_population:
        - connector_name: string, short descriptor
        - from_population: Population object which has newcomers. Its future
        is not affected by this connector.
        - to_population: Population object that receives the newcomers
        entering from_population.
        - scale_factor: Parameter object that multiplies expectation. Data treats fraction as binomial
        - ratio_populations: List of 2 populations, the ratio of those populations is applied as
        a scale factor
    """

    def __init__(self, connector_name: str, from_population: Population, to_population: Population,
                 scale_factor: Parameter = None, ratio_populations: list = None):
        """Constructor
        """
        super().__init__(connector_name, from_population, to_population)

        if not isinstance(from_population, Population):
            raise TypeError('Adder('+self.name+
                            ') from_population must be a Population object')

        if not isinstance(to_population, Population):
            raise TypeError('Adder ('+self.name+
                            ') to_population must be a Population object')

        if scale_factor is not None:
            if not isinstance(scale_factor, Parameter):
                raise TypeError('Adder('+self.name+
                            ') scale_factor must be a Parameter object')
        self.scale_factor = scale_factor

        if ratio_populations is not None:
            if not isinstance(ratio_populations, list) or len(ratio_populations) != 2:
                raise TypeError('Adder('+self.name+
                            ') ratio_populations must be a list of two Population objects')
            for i in range(2):
                if not isinstance(ratio_populations[i], Population):
                    raise TypeError('Adder(' + self.name +
                                    ') ratio_populations must be a list of two Population objects')
        self.ratio_populations = ratio_populations

    def update_expectation(self):
        """
        do addition for expectation
        """
        newcomers = 0
        if len(self.from_population.future) > 0:
            newcomers = self.from_population.future[0]

        if newcomers > 0:
            if getattr(self, "scale_factor", None) is not None:
                newcomers = self.scale_factor.get_value() * newcomers
            if getattr(self, "ratio_populations", None) is not None:
                if self.ratio_populations[1].history[-1] > 0.:
                    ratio = self.ratio_populations[0].history[-1]/self.ratio_populations[1].history[-1]
                    newcomers = ratio * newcomers

            self.to_population.update_future_fast(newcomers)

    def update_data(self):
        """
        do addition for data
        """
        newcomers = 0
        if len(self.from_population.future) > 0:
            newcomers = self.from_population.future[0]

        if newcomers > 0:
            scale = 1
            if getattr(self,"scale_factor",None) is not None:
                scale = self.scale_factor.get_value() * scale
            if getattr(self, "ratio_populations", None) is not None:
                if self.ratio_populations[1].history[-1] > 0.:
                    ratio = 1. * self.ratio_populations[0].history[-1] / self.ratio_populations[1].history[-1]
                    scale = ratio * scale

            if scale != 1:
                iscale = int(scale)
                i_newcomers = iscale*newcomers
                fscale = scale - iscale
                f_newcomers = stats.binom.rvs(newcomers, fscale)
                newcomers = i_newcomers + f_newcomers

            self.to_population.update_future_fast(newcomers)

