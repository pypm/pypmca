# -*- coding: utf-8 -*-
"""
Connector: A base'virtual' class for defining connections between populations.

The derived classes must define the methods for updating the future values of populations.

Depending on the connector type, from_population and to_population can be a
single population or a list of populations.

@author:karlen
"""

from Population import Population

class Connector:
    """
    Base class for defining connections between populations
    """

    def __init__(self, connector_name, from_population, to_population):
        """Constructor
        """
        self.name = str(connector_name)

        pops = {'from_population':from_population,
                'to_population':to_population}
        for key in pops:
            if isinstance(pops[key], list):
                for pop in pops[key]:
                    if not isinstance(pop, Population):
                        raise TypeError('Error in constructing '+self.name+
                                        ': '+key+' list must contain Population objects')
            else:
                if not isinstance(pops[key], Population):
                    raise TypeError('Error in constructing '+self.name+
                                    ': '+key+' argument must be a Population object')

        self.from_population = from_population
        self.to_population = to_population
        self.parameters = {}

    def __str__(self):
        return self.name

    def update_expectation(self):
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

    def update_data(self):
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
