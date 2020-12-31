# -*- coding: utf-8 -*-
"""
Subtractor: A Connector class that 'immediately' removes population.

The purpose of this connector is to include the effect of removal from
populations. For example it is used to track those currently in hospital, instead
of just monitoring the total hospitalized.

A scale factor can be applied. The scale factor should be positive.

@author: karlen
"""
from scipy import stats

from pypmca.Connector import Connector
from pypmca.Population import Population
from pypmca.Parameter import Parameter

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
            raise TypeError('Subtractor (' + self.name +
                            ') from_population must be a Population object')

        if not isinstance(to_population, Population):
            raise TypeError('Subtractor (' + self.name +
                            ') to_population must be a Population object')

        if hasattr(from_population, 'report_noise') and from_population.report_noise:
            raise TypeError('Subtractor (' + self.name +
                            ') from_population cannot have report_noise set to True')

        if scale_factor is not None:
            if not isinstance(scale_factor, Parameter):
                raise TypeError('Subtractor('+self.name+
                            ') scale_factor must be a Parameter object')
        self.scale_factor = scale_factor

        if ratio_populations is not None:
            if not isinstance(ratio_populations, list) or len(ratio_populations) != 2:
                raise TypeError('Subtractor('+self.name+
                            ') ratio_populations must be a list of two Population objects')
            for i in range(2):
                if not isinstance(ratio_populations[i], Population):
                    raise TypeError('Subtractor(' + self.name +
                                    ') ratio_populations must be a list of two Population objects')
        self.ratio_populations = ratio_populations

        # Identify the from_population as one not to derive daily numbers
        self.from_population.monotonic = False

    def update_expectation(self):
        """
        do subtraction for expectation
        """
        reduction = 0
        if len(self.to_population.future) > 0:
            reduction = -1 * self.to_population.future[0]

        if reduction < 0:
            if getattr(self, "scale_factor", None) is not None:
                reduction = self.scale_factor.get_value() * reduction
            if getattr(self, "ratio_populations", None) is not None:
                if self.ratio_populations[1].history[-1] > 0.:
                    ratio = self.ratio_populations[0].history[-1] / self.ratio_populations[1].history[-1]
                    reduction = ratio * reduction

            self.from_population.update_future_fast(reduction)

    def update_data(self):
        """
        do subtraction
        """
        reduction = 0
        if len(self.to_population.future) > 0:
            reduction = -1 * self.to_population.future[0]

        if reduction < 0:
            scale = 1
            if getattr(self, "scale_factor", None) is not None:
                scale = self.scale_factor.get_value() * scale
            if getattr(self, "ratio_populations", None) is not None:
                if self.ratio_populations[1].history[-1] > 0.:
                    ratio = 1. * self.ratio_populations[0].history[-1] / self.ratio_populations[1].history[-1]
                    scale = ratio * scale

            if scale != 1:
                iscale = int(scale)
                i_reduction = iscale * reduction
                fscale = scale - iscale
                f_reduction = -stats.binom.rvs(-reduction, fscale)
                reduction = i_reduction + f_reduction

            self.from_population.update_future_fast(reduction)
