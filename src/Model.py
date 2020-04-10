# -*- coding: utf-8 -*-
"""
Model: A class consisting of a list of populations and connectors.

With this information, the evolution of each population is performed,
either by updating expectations or simulated data.

@author: karlen
"""
from datetime import date

from Connector import Connector
from Transition import Transition

class Model:
    """
    Class that defines populations and their evolution through connectors
    """

    def __init__(self, model_name):

        self.name = str(model_name)
        self.t0 = date.today()
        self.__time_step = 1. # default time step is one day. (hidden)
        self.__time_step_changed = False
        self.populations = {}
        self.connectors = {}
        self.connector_list = []
        self.transitions = {}

        self.boot_pars = {}
        self.boot_needed = True

    def set_t0(self, year=2020, month=1, day=1):
        """Define the date of the first element in the history lists
        """
        self.t0 = date(year, month, day)

    def set_time_step(self, time_step):
        """
        Define the length of time corresponding to one step in the evolution

        Parameters
        ----------
        time_step : size of time step in units of days. So a 1 hour time step
                    is specified as 1./24.

        Returns
        -------
        None.

        """
        if not isinstance(time_step, float):
            raise TypeError('Model ('+self.name+
                            ' time_step must be a float')
        if not time_step > 0.:
            raise ValueError('Model ('+self.name+
                             ' time_step must be positive')
        self.__time_step = time_step
        self.__time_step_changed = True

    def get_time_step(self):
        """
        Get the length of time corresponding to one step in the evolution

        Returns
        -------
        real value for one time step

        """
        return self.__time_step

    def boot_setup(self, boot_population, boot_value, exclusion_populations=None):
        """
        Typically most populations are initialized to zero. Propogation and
        growith requires at least one population to not be zero. However,
        starting a model with large values in a source population can lead to
        unusual behaviour for initial time steps, as that would be inconsistent
        with the zeros in most other populations. On the other hand, starting
        with small values in a source population can lead to extreme sensitivity to
        time. To solve this issue, we boot the model, by starting with
        small populations and evolve the system until a population exceeds
        a target value. History is removed and the current populations and futures
        are used for the initial state of the model. To allow for
        continuous dependence, the current value and futures are scaled by the ratio of
        the target value to actual value, except for the exclusion_populations.
            - boot_population: str, source population
            - boot_value: float or int, small value to begin boot sequence.
            - exclusion_populations: str or str-list of populations to exclude
            from current value scaling. The history is set to the initial value for
            the exclusion populations, the future contributions are scaled.
        The boot sequence ends when the boot_population reaches or exceeds
        the value set for the initial_value of that population.

        One downside of the rescaling is that the long term of contagious
        population can be non-zero (a few).

        It can be useful to save a copy of a booted model, and perform restarts
        of that model.
        """

        name_boot_pop = str(boot_population)
        if name_boot_pop not in self.populations:
            raise ValueError('Error in boot_setup for model ('+self.name+
                             '). '+name_boot_pop+'is not present in model.')
        exc_pop_list = []
        if exclusion_populations is not None:
            if isinstance(exclusion_populations, list):
                for exc_pop in exclusion_populations:
                    name_exc_pop = str(exc_pop)
                    if name_exc_pop not in self.populations:
                        raise ValueError('Error in boot_setup for model ('+self.name+
                                         '). '+name_exc_pop+'is not present in model.')
                    exc_pop_list.append(name_exc_pop)
            else:
                name_exc_pop = str(exclusion_populations)
                if name_exc_pop not in self.populations:
                    raise ValueError('Error in boot_setup for model ('+self.name+
                                     '). '+name_exc_pop+'is not present in model.')
                exc_pop_list.append(name_exc_pop)

        self.boot_pars['boot_population'] = name_boot_pop
        self.boot_pars['boot_value'] = boot_value
        self.boot_pars['exclusion_populations'] = exc_pop_list

    def reset(self):
        """
        Reinitialize all populations, eliminating history and future, and
        setting initial values. Also reset transitioned parameters
        """
        for key in self.populations:
            self.populations[key].reset()

        for key in self.transitions:
            self.transitions[key].reset()

        self.boot_needed = True

    def boot(self, expectations=True):
        """ Using information in boot_pars, evolve the system to a point where
        the target population is reached or exceeded:
            - expectations: if True, history begins with expectations,
            otherwise, a Poisson number is drawn using that expectation.
        """
        if len(self.boot_pars) == 0:
            raise RuntimeError('Model boot parameters not set.')

        # start from a clean slate of populations set to their initial values
        self.reset()

        # need to set this to False before evolving expectations
        self.boot_needed = False

        boot_pop = self.populations[self.boot_pars['boot_population']]
        goal_value = boot_pop.history[-1]
        boot_pop.history[-1] = self.boot_pars['boot_value']

        # exclusion_populations are reset to their original values after the
        # boot. Save them here.
        exclusions = self.boot_pars['exclusion_populations']
        exc_vals = {}
        for exc_pop_name in exclusions:
            exc_pop = self.populations[exc_pop_name]
            exc_vals[exc_pop_name] = exc_pop.history[-1]

        last_value = -1
        while boot_pop.history[-1] < goal_value and boot_pop.history[-1] >= last_value:
            last_value = boot_pop.history[-1]
            self.evolve_expectations(1)

        scale = goal_value/boot_pop.history[-1]

        #Erase history until current. Scale all histories and futures.
        for key in self.populations:
            pop = self.populations[key]
            pop.remove_history()
            if exclusions is None or str(pop) not in exclusions:
                pop.scale_history(scale, expectations)
            elif exclusions is not None and str(pop) in exclusions:
                pop.history[-1] = exc_vals[str(pop)]
            pop.scale_future(scale, expectations)

        # restore goal_value for data production
        if not expectations:
            boot_pop.history[-1] = int(goal_value)

    def evolve_expectations(self, n_step):
        """
        First check if a boot is needed: this is required if
        reset_populations is called

        Take n_steps of
        0) check if it is time to apply a transition
        1) updating futures
        2) appending new value to history array
        3) remove current time step in future array
        """
        if self.boot_needed:
            self.boot(expectations=True)

        for step in range(n_step):
            for transition_name in self.transitions:
                transition = self.transitions[transition_name]
                if step == transition.trigger_step:
                    transition.take_action(expectations=True)
            # calculate future expectations
            for connector_name in self.connector_list:
                connector = self.connectors[connector_name]
                connector.update_expectation()
            # make one time step
            for key in self.populations:
                self.populations[key].do_time_step()

    def generate_data(self, n_step):
        """
        Produce data
        """
        if self.boot_needed:
            self.boot(expectations=False)

        for step in range(n_step):
            for transition_name in self.transitions:
                transition = self.transitions[transition_name]
                if step == transition.trigger_step:
                    transition.take_action(expectations=False)
            # calculate future data
            for connector_name in self.connector_list:
                connector = self.connectors[connector_name]
                connector.update_data()
            # make one time step
            for key in self.populations:
                self.populations[key].do_time_step()

    def add_transition(self, transition):
        """
        Add a transition to the model. Special action is taken in the specied
        time step. Transitions do not apply to the boot phase.
        """
        if not isinstance(transition, Transition):
            raise TypeError('Error in adding transition to '+self.name+
                            ': add_transition argument must be a Transition object')
        if str(transition) in self.transitions:
            raise ValueError('Error in adding transition ('+
                             str(transition)+') to model ('+self.name+
                             '). Transition with that name already present in model.')
        self.transitions[str(transition)] = transition

    def add_connector(self, connector, before_connector=None, after_connector=None):
        """
        Add a connector to the model. The model is evolved by processing
        the connectors in the order that they appear in connector_list.
        Connectors must have unique names.

        Parameters
        ----------
        connector : Connector object
        """
        if not isinstance(connector, Connector):
            raise TypeError('Error in adding connector to '+self.name+
                            ': add_connector argument must be a Connector object')
        if str(connector) in self.connectors:
            raise ValueError('Error in adding connector ('+
                             str(connector)+') to model ('+self.name+
                             '). Connector with that name already present in model.')
        if after_connector is not None and before_connector is not None:
            raise ValueError('Error in adding connector('+
                             str(connector)+') to model ('+self.name+
                             '). Cannot specify both before and after locations.')
        index_loc = -1
        if before_connector is not None:
            name = str(before_connector)
            if name not in self.connector_list:
                raise ValueError('Error in adding connector('+
                                 str(connector)+') to model ('+self.name+
                                 '): '+str(before_connector)+' not found.')
            index_loc = self.connector_list.index(name)
        elif after_connector is not None:
            name = str(after_connector)
            if name not in self.connector_list:
                raise ValueError('Error in adding connector('+
                                 str(connector)+') to model ('+self.name+
                                 '): '+str(after_connector)+' not found.')
            index_loc = self.connector_list.index(name) + 1
        else:
            index_loc = len(self.connector_list)

        self.connector_list.insert(index_loc, str(connector))
        self.connectors[str(connector)] = connector

        self.__update_population_list(connector)

    def remove_connector(self, connector):
        """
        Remove a connector from the model.

        Parameters
        ----------
        connector : Connector or str
            Connector to be removed. Can be specified by its name.

        Returns
        -------
        removed_connector : Connector
            Returned so that it can be easily added back.
        """

        name = str(connector)
        if name not in self.connector_list:
            raise ValueError('Error in removing connector('+
                             str(connector)+') from model ('+self.name+
                             '): '+str(connector)+' not found.')

        self.connector_list.remove(name)
        removed_connector = self.connectors.pop(name, None)

        # remake list of active populations
        self.populations = {}
        for con_name in self.connector_list:
            con = self.connectors[con_name]
            self.__update_population_list(con)

        return removed_connector

    def __update_population_list(self, connector):

        # make a list of populations involved in this connector
        pops = []
        args = [connector.from_population, connector.to_population]
        for arg in args:
            if isinstance(arg, list):
                for pop in arg:
                    pops.append(pop)
            else:
                pops.append(arg)
        if hasattr(connector, 'chain'):
            chain = connector.chain
            for link in chain:
                pops.append(link.from_population)
                pops.append(link.to_population)

        # If necessary add to population dictionary.
        # Check there are no duplicate keys
        for pop in pops:
            key = pop.name
            if key in self.populations:
                dict_pop = self.populations[key]
                if dict_pop != pop:
                    raise ValueError('Two populations share the same name: '+
                                     pop.name)
            else:
                self.populations[key] = pop
