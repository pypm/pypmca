# -*- coding: utf-8 -*-
"""
Model: A class consisting of a list of populations and connectors.

With this information, the evolution of each population is performed,
either by updating expectations or simulated data.

@author: karlen
"""
from datetime import date
import pickle
import copy
from pathlib import Path
from pypmca.Connector import Connector
from pypmca.Population import Population
from pypmca.Transition import Transition
from pypmca.Parameter import Parameter

class Model:
    """
    Class that defines populations and their evolution through connectors
    """

    def __init__(self, model_name):

        self.name = str(model_name)
        self.description = ''
        self.t0 = date(2020, 3, 1)
        # default time step is one day. (hidden)
        self.__time_step = 1.
        self.populations = {}
        self.connectors = {}
        self.connector_list = []
        self.transitions = {}
        self.parameters = {}
        self.delays = {}
        self.lnl = None
        self.population_mcmc_time_step_range = {}
        self.user_dict = {}

        self.boot_pars = {}
        self.boot_needed = True

    def set_t0(self, year=2020, month=1, day=1):
        """Define the date of the first element in the history lists
        """
        self.t0 = date(year, month, day)

    def set_time_step(self, time_step):
        """
        Define the length of time corresponding to one step in the evolution.
        Update all delays and transitions if this is changed.

        Parameters
        ----------
        time_step : size of time step in units of days. So a 1 hour time step
                    is specified as 1./24.

        Returns
        -------
        None.

        """
        if not isinstance(time_step, float):
            raise TypeError('Model (' + self.name +
                            ' time_step must be a float')
        if not time_step > 0.:
            raise ValueError('Model (' + self.name +
                             ' time_step must be positive')
        self.__time_step = time_step

        # update the objects that are sensitive to time_step
        for delay_name in self.delays:
            delay = self.delays[delay_name]
            delay.update()
        for transition_name in self.transitions:
            transition = self.transitions[transition_name]
            transition.update()
        for connector_name in self.connectors:
            connector = self.connectors[connector_name]
            if type(connector).__name__ == 'Chain':
                connector.update()

    def get_time_step(self):
        """
        Get the length of time corresponding to one step in the evolution

        Returns
        -------
        real value for one time step

        """
        return self.__time_step

    def get_history(self, population_name):
        if not isinstance(population_name, str):
            raise TypeError('model.get_history argument must be a str')

        pop = self.populations[population_name]
        return pop.history

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

        To fit data to a model it will be necessary to set the initial size of
        a critical population (the boot population) to match the scale for data.
        For that reason, we insist that the boot population be initialized by
        a parameter below.

        One downside of the rescaling is that the long term of contagious
        population can be non-zero (a few).

        It can be useful to save a copy of a booted model, and perform restarts
        of that model.
        """

        name_boot_pop = str(boot_population)
        if name_boot_pop not in self.populations:
            raise ValueError('Error in boot_setup for model (' + self.name +
                             '). ' + name_boot_pop + 'is not present in model.')
        initial_value = self.populations[name_boot_pop].initial_value
        if not isinstance(initial_value, Parameter):
            raise ValueError('Error in boot_setup for model (' + self.name +
                             '). ' + name_boot_pop + 'is not initialized by a parameter.')
        exc_pop_list = []
        if exclusion_populations is not None:
            if isinstance(exclusion_populations, list):
                for exc_pop in exclusion_populations:
                    name_exc_pop = str(exc_pop)
                    if name_exc_pop not in self.populations:
                        raise ValueError('Error in boot_setup for model (' + self.name +
                                         '). ' + name_exc_pop + 'is not present in model.')
                    exc_pop_list.append(name_exc_pop)
            else:
                name_exc_pop = str(exclusion_populations)
                if name_exc_pop not in self.populations:
                    raise ValueError('Error in boot_setup for model (' + self.name +
                                     '). ' + name_exc_pop + 'is not present in model.')
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

        # do this in reverse order, so earlier transitions get precidence
        trans_dict = {}
        rc = 0.
        for key in self.transitions:
            rc += 0.01
            trans = self.transitions[key]
            trig_step = trans.trigger_step + rc
            trans_dict[trig_step] = trans
        trans_list = list(trans_dict.keys())
        trans_list.sort()
        trans_list.reverse()
        for trans_step in trans_list:
            trans = trans_dict[trans_step]
            if trans.enabled:
                trans.reset()

        self.boot_needed = True

    def boot(self, expectations=True):
        """ Using information in boot_pars, evolve the system to a point where
        the target population is reached or exceeded:
            - expectations: if True, after the boot completes, the
            history will begin with an expectation,
            otherwise, a Poisson number is drawn using that expectation.
        Note: check is added so that if boot requirement is not made within
        10k steps, it fails (rather than stay in infinite loop)
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
        i_step = 0
        while boot_pop.history[-1] < goal_value and boot_pop.history[-1] >= last_value:
            last_value = boot_pop.history[-1]
            self.evolve_expectations(1)
            i_step += 1
            if i_step % 1000 == 0:
                if boot_pop.history[-1] < 0.1 * goal_value:
                    self.boot_needed = True
                    raise ValueError('The boot process is taking too long to complete. ' +
                                     'Adjust parameters to allow the goal to be reached.')

        if boot_pop.history[-1] < goal_value:
            raise ValueError('The boot process did not reach its target. Consider reduce the target.')

        scale = goal_value / boot_pop.history[-1]

        # Erase history until current. Scale all histories and futures.
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

    def evolve_expectations(self, n_step, from_step=0):
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

        for step in range(from_step, n_step+from_step):
            self.do_transitions(step, expectations=True)

            # calculate future expectations
            self.calculate_future(expectations=True)

            # make one time step
            self.do_time_step(expectations=True)

    def generate_data(self, n_step, from_step=0, data_start=0):
        """
        Produce data:
        - n_step: number of steps to take
        - from_step: allows a continuation of a previous call to generate_data: no check that this was done
        - data_start: allows that model is first evolved (expectation values) to a date, then generate data
        """
        if data_start >= from_step and data_start>0:
            n_exp_step = data_start - from_step
            if n_exp_step > 0:
                self.evolve_expectations(n_exp_step, from_step=from_step)

            # Convert current value and future values from float to int, for all populations
            # in order to produce data going forward
            for key in self.populations:
                pop = self.populations[key]
                nu = pop.history[-1]
                pop.history[-1] = int(round(nu))
                pop.scale_future(1., expectations=False)

            data_steps = n_step + from_step - data_start
            self.generate_data(data_steps, from_step=data_start, data_start=0)
            return

        if self.boot_needed:
            self.boot(expectations=False)

        for step in range(from_step, n_step+from_step):
            self.do_transitions(step, expectations=False)

            # calculate future data
            self.calculate_future(expectations=False)

            # make one time step
            self.do_time_step(expectations=False)

    def do_transitions(self, step, expectations=True):
        for transition_name in self.transitions:
            transition = self.transitions[transition_name]
            if transition.enabled:
                linear = getattr(transition,'linear',False)
                if linear:
                    if step >= transition.trigger_step:
                        n_step_par = getattr(transition,'n_step',None)
                        n_step = -1
                        if n_step_par is not None:
                            n_step = n_step_par.get_value()
                        if n_step <0 or n_step > step - transition.trigger_step:
                            transition.take_action(expectations)
                else:
                    if step == transition.trigger_step:
                        transition.take_action(expectations)

    def calculate_future(self, expectations=True):
        for connector_name in self.connector_list:
            connector = self.connectors[connector_name]
            if expectations:
                connector.update_expectation()
            else:
                connector.update_data()

    def do_time_step(self, expectations=True):
        for key in self.populations:
            self.populations[key].do_time_step(expectations)

    def add_transition(self, transition):
        """
        Add a transition to the model. Special action is taken in the specified
        time step. Transitions do not apply to the boot phase.
        """
        if not isinstance(transition, Transition):
            raise TypeError('Error in adding transition to ' + self.name +
                            ': add_transition argument must be a Transition object')
        if str(transition) in self.transitions:
            raise ValueError('Error in adding transition (' +
                             str(transition) + ') to model (' + self.name +
                             '). Transition with that name already present in model.')
        self.transitions[str(transition)] = transition

        self.__update_parameter_list(transition)

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
            raise TypeError('Error in adding connector to ' + self.name +
                            ': add_connector argument must be a Connector object')
        if str(connector) in self.connectors:
            raise ValueError('Error in adding connector (' +
                             str(connector) + ') to model (' + self.name +
                             '). Connector with that name already present in model.')
        if after_connector is not None and before_connector is not None:
            raise ValueError('Error in adding connector(' +
                             str(connector) + ') to model (' + self.name +
                             '). Cannot specify both before and after locations.')
        index_loc = -1
        if before_connector is not None:
            name = str(before_connector)
            if name not in self.connector_list:
                raise ValueError('Error in adding connector(' +
                                 str(connector) + ') to model (' + self.name +
                                 '): ' + str(before_connector) + ' not found.')
            index_loc = self.connector_list.index(name)
        elif after_connector is not None:
            name = str(after_connector)
            if name not in self.connector_list:
                raise ValueError('Error in adding connector(' +
                                 str(connector) + ') to model (' + self.name +
                                 '): ' + str(after_connector) + ' not found.')
            index_loc = self.connector_list.index(name) + 1
        else:
            index_loc = len(self.connector_list)

        self.connector_list.insert(index_loc, str(connector))
        self.connectors[str(connector)] = connector

        self.__update_population_list(connector)
        self.__update_parameter_list(connector)
        self.__update_delay_list(connector)

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
            raise ValueError('Error in removing connector(' +
                             str(connector) + ') from model (' + self.name +
                             '): ' + str(connector) + ' not found.')

        self.connector_list.remove(name)
        removed_connector = self.connectors.pop(name, None)

        self.update_lists()
        return removed_connector

    def remove_transition(self, transition):
        name = str(transition)
        if name not in self.transitions:
            raise ValueError('Error in removing transition(' +
                             str(transition) + ') from model (' + self.name +
                             '): ' + str(transition) + ' not found.')

        removed_transition = self.transitions.pop(name, None)

        self.update_lists()
        return removed_transition

    def update_lists(self):
        """ After modifying any element of the model that involves
            populations or parameters, call this to be sure that
            the list of active populations and parameters is current
        """
        # remake list of active populations and their parameters
        self.populations = {}
        self.parameters = {}
        self.delays = {}

        for con_name in self.connector_list:
            con = self.connectors[con_name]
            self.__update_population_list(con)
            self.__update_parameter_list(con)
            self.__update_delay_list(con)

        # add to list of active parameters
        for trans_name in self.transitions:
            trans = self.transitions[trans_name]
            self.__update_parameter_list(trans)

    def __update_parameter_list(self, obj):

        if type(obj).__name__ == 'Injector':
            self.__add_to_parameter_list(obj.transition_time)
            self.__add_to_parameter_list(obj.injection)

        elif type(obj).__name__ == 'Modifier':
            self.__add_to_parameter_list(obj.transition_time)
            self.__add_to_parameter_list(obj.parameter_before)
            self.__add_to_parameter_list(obj.parameter_after)
            if getattr(obj,'n_step',None) is not None:
                self.__add_to_parameter_list(obj.n_step)

        elif isinstance(obj, Connector):
            con = obj
            for par_key in con.parameters:
                par = con.parameters[par_key]
                self.__add_to_parameter_list(par)
            if type(con).__name__ == 'Chain':
                for link_con in con.chain:
                    for link_par_key in link_con.parameters:
                        link_par = link_con.parameters[link_par_key]
                        self.__add_to_parameter_list(link_par)

        elif isinstance(obj, Population):
            pop = obj
            for par_key in pop.parameters:
                par = pop.parameters[par_key]
                self.__add_to_parameter_list(par)

    def __add_to_parameter_list(self, obj):
        # check that the name of the parameter is unique
        if str(obj) in self.parameters:
            par = self.parameters[str(obj)]
            if obj != par:
                raise ValueError('Two parameters share the same name: ' + str(obj))
        else:
            self.parameters[str(obj)] = obj

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
                    raise ValueError('Two populations share the same name: ' +
                                     pop.name)
            else:
                self.populations[key] = pop
                self.__update_parameter_list(pop)
                # inform the population of its parent model - to get t0 and timestep information
                pop.set_model(self)

    def __update_delay_list(self, connector):
        # make a list of delays involved in this connector
        delays = []
        if hasattr(connector, 'delay'):
            if isinstance(connector.delay, list):
                for delay in connector.delay:
                    delays.append(delay)
                    if type(connector).__name__ == 'Chain':
                        for chain_con in connector.chain:
                            delays.append(chain_con.delay)

            else:
                delays.append(connector.delay)
                if type(connector).__name__ == 'Chain':
                    for chain_con in connector.chain:
                        delays.append(chain_con.delay)

        # If necessary add to delay dictionary.
        # Check there are no duplicate keys
        for delay in delays:
            key = delay.name
            if key in self.delays:
                dict_delay = self.delays[key]
                if dict_delay != delay:
                    raise ValueError('Two delays share the same name: ' +
                                     delay.name)
            else:
                self.delays[key] = delay

    def copy_values_from(self, from_model):
        """Copy the parameter values and transition states from another model to this one
           Also copy model.name and description. Also copy parameter/population hidden info"""
        if not isinstance(from_model, Model):
            raise ValueError('Error in copy_values_from. Argument must be a model object')

        for par_name in from_model.parameters:
            from_par = from_model.parameters[par_name]
            if par_name in self.parameters:
                par = self.parameters[par_name]
                par.set_value(from_par.get_value())
                par.initial_value = from_par.initial_value
                par.hidden = from_par.hidden
                par.mcmc_step = from_par.mcmc_step
                try:
                    par.std_estimator = from_par.std_estimator
                except:
                    pass
                par.parameter_min = from_par.parameter_min
                par.parameter_max = from_par.parameter_max
                if from_par.prior_function is not None:
                    prior_function = copy.copy(from_par.prior_function)
                if from_par.prior_parameters is not None:
                    prior_parameters = copy.deepcopy(from_par.prior_parameters)
                if from_par.get_status() == 'variable':
                    par.set_variable(prior_function, prior_parameters)
                else:
                    if from_par.prior_function is not None:
                        par.prior_function = prior_function
                    if from_par.prior_parameters is not None:
                        par.prior_parameters = prior_parameters

        for tran_name in from_model.transitions:
            from_tran = from_model.transitions[tran_name]
            if tran_name in self.transitions:
                tran = self.transitions[tran_name]
                tran.enabled = from_tran.enabled

        for pop_name in from_model.populations:
            from_pop = from_model.populations[pop_name]
            if pop_name in self.populations:
                pop = self.populations[pop_name]
                pop.hidden = from_pop.hidden
                # assign initial population value as necessary
                pop.set_initial_value(from_pop.initial_value)

        self.name = 'copy: '+from_model.name
        self.description = from_model.description

        self.t0 = from_model.t0
        self.set_time_step(from_model.get_time_step())
        if getattr(from_model, "user_dict", None) is not None:
            self.user_dict = from_model.user_dict

        self.boot_pars['boot_population'] = from_model.boot_pars['boot_population']
        self.boot_pars['boot_value'] = from_model.boot_pars['boot_value']
        # self.boot_pars['exclusion_populations'] = exc_pop_list
        self.boot_needed = True

    def save_file(self, filename):
        """
        Save a copy of the current model to a file. The model can be restored
        at a late time using Model.open_file(filename).

        Before saving, the model futures are removed, and the initial_values
        are set to the current values for all parameters.

        If no extension is provided, the default extension, .pypm is added.

        Parameters
        ----------
        filename : Path or str
            name of file to save model

        Returns
        -------
        None.

        """

        try:
            filepath = Path(filename).resolve()
        except:
            raise TypeError('Input arg could not be converted to a valid path: {}' +
                            '\n It must be a str or Path-like.'.format(filename))

        if len(filepath.suffix) < 2:
            filepath = filepath.with_suffix('.pypm')

        #       eliminate histories, undo transitioned parameter adjustments

        self.reset()

        #       set initial parameter values to their current values to define a new
        #       reference

        for par_name in self.parameters:
            par = self.parameters[par_name]
            par.new_initial_value()

        with open(filepath, 'wb') as f:
            pickle.dump(self, f, protocol=4)


    @classmethod
    def open_file(cls, filepath):
        """
        Restore a model that was save to a file using Model.save_file(filename)

        Parameters
        ----------
        filepath : Path or str
            name of existing model file to open

        Returns
        -------
        Model
            The model object saved in the file

        """

        try:
            filepath = Path(filepath).resolve()
        except:
            raise TypeError('Input arg could not be converted to a valid path: {}' +
                            '\n It must be a str or Path-like.'.format(filepath))

        if not filepath.exists():
            raise ValueError('Filepath does not exist: {}'.format(filepath))

        if len(filepath.suffix) < 2:
            filepath = filepath.with_suffix('.pypm')

        with open(filepath, 'rb') as f:
            return pickle.load(f)
