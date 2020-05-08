# -*- coding: utf-8 -*-
"""
Ensemble: A Model class whose behaviour is defined by an ensemble of models

An ensemble allows for categorization of populations by age, risk, or other factors. It can also
be used to combine many models that are nearly independent (such as separate provinces to make a
Canada wide model).

Each category has its own model, and one model can influence the others.
The ensemble sums the histories of all its models to represent the evolution of
the entire system.

@author: karlen
"""
from datetime import date
import numpy as np
from scipy import stats

import pickle
from pypm import Model, Parameter, Population


class Ensemble(Model):
    """
    Ensemble overwrites some model methods, so that actions are applied to all models.
    It's own connectors are not used for evolution: each model is evolved and ensemble
    histories are just the sum of the models histories.
    """

    def __init__(self, name: str, model):
        """Constructor
         - Model specifies the reference model: populations etc should match those in
         the ensemble of models. It can be specified as a model, or as a filename
        """
        super().__init__(name)

        self.model = None
        if isinstance(model, Model):
            self.model = model
        elif isinstance(model, str):
            self.model = Model.open_file(model)
        else:
            raise TypeError('Error creating ensemble. The model argument needs to be a model object ' +
                            'or a .pypm filename')

        # use the population objects in self.model to save the histories
        # this allows for consistent access to history between models and ensembles
        for pop_name in self.model.populations:
            pop = self.model.populations[pop_name]
            self.populations[pop_name] = pop

        # The list of models (either provided directly or as a list of .pypm files)
        # They are required to have unique names for user interaction
        # Ordering is preserved to match the transmission matrix
        self.model_list = []
        self.models = {}
        self.infected_name = None
        self.susceptible_name = None
        self.total_name = None
        self.contagious_name = None
        self.alpha_name = None
        self.infection_cycle_name = None
        self.contact_matrix = None
        self.null_pop = Population('null', 0.)

        self.__distribution = None
        self.__nbinom_par = None
        self.set_distribution('poisson', None)

    def set_distribution(self, distribution, nbinom_par):
        if distribution not in ['poisson', 'nbinom']:
            raise ValueError('Ensemble (' + self.name +
                             ') distribution must be poisson or nbinom')
        self.__distribution = distribution

        if distribution == 'nbinom' and \
                (nbinom_par is None or not isinstance(nbinom_par, Parameter)):
            raise TypeError('Ensemble (' + self.name +
                            ') nbinom_par must be a Parameter object')
        self.__nbinom_par = nbinom_par

    def get_distribution(self):
        return self.__distribution, self.__nbinom_par

    def upload_models(self, models: list):
        if not isinstance(models, list):
            raise TypeError('Error in upload_models: argument must be a list of Model objects')
        for model in models:
            if not isinstance(model, Model):
                raise TypeError('Error in upload_models: argument must be a list of Model objects')
        for model in models:
            self.__add_model(model)

    def upload_model_files(self, model_filenames: list):
        if not isinstance(model_filenames, list):
            raise TypeError('Error in upload_model_files: argument must be a list of .pypm filenames')
        for model_file in model_filenames:
            if not isinstance(model_file, str):
                raise TypeError('Error in upload_model_files: argument must be a list of .pypm filenames')

        for model_filename in model_filenames:
            self.__add_model(Model.open_file(model_filename))

    def __add_model(self, model):
        if model.name in self.models:
            raise ValueError('Error: Two models have same name: ' + model.name)
        for population_name in self.populations:
            if population_name not in model.populations:
                raise ValueError('Error: Model ' + model.name + ' does not contain population:' + population_name)

        self.model_list.append(model.name)
        self.models[model.name] = model
        self.reset_cross_transmission()

    def define_cross_transmission(self, infection_cycle_name: str, infected_name: str,
                                  susceptible_name: str, total_name: str,
                                  contagious_name: str, alpha_name: str, contact: list):
        """
        Define how the populations mix to produce new infections:
        For example new infections in model A arising from interactions with group B
        are calculated by Susceptible_A / M * f[A][B] * Contagious_B * alpha_AB
            - M is the total effective total population = sum f[A][B] * N_B.
            - alpha_AB is the geometric mean = sqrt(alpha_A * alpha_B)
            - the contact matrix, f, needs to be n x n, where n is the number of models
            and the order of the matrix rows matches the order of the models in model_list
        In a homogeneous society, all alphas are the same, and all terms of the matrix are 1.
        f represents the relative probability for a contagious members of B group to infect a
        random member of the A group (relative to infecting a random member of the B group).
        The off diagonal elements are typically less than 1. The diagonal elements are 1 by
        definition. A diagonal matrix describes a set of independent populations.
        Given the definintion of M, both homogenous and independent populations yield their
        original infection rate without adjusting alpha.
        """
        if infection_cycle_name not in self.model.connectors:
            raise ValueError('Infection cycle (' + infection_cycle_name + ') not found in connectors')
        self.infection_cycle_name = infection_cycle_name
        # these could have been derived from belpw, but for non-symmetric contact matricies
        # the order of the two multipliers matter.
        if infected_name not in self.populations:
            raise ValueError('Population ' + infected_name + ' not found in reference model.')
        self.infected_name = infected_name
        if susceptible_name not in self.populations:
            raise ValueError('Population ' + susceptible_name + ' not found in reference model.')
        self.susceptible_name = susceptible_name
        if total_name not in self.populations:
            raise ValueError('Population ' + total_name + ' not found in reference model.')
        self.total_name = total_name
        if contagious_name not in self.populations:
            raise ValueError('Population ' + contagious_name + ' not found in reference model.')
        self.contagious_name = contagious_name
        if alpha_name not in self.model.parameters:
            raise ValueError('Parameter ' + alpha_name + ' not found in reference model.')
        self.alpha_name = alpha_name

        n_model = len(self.model_list)
        if not isinstance(contact, list):
            raise TypeError('The contact matrix must be a list of lists')
        if len(contact) != n_model:
            raise ValueError('The contact matrix must be an nxn matrix, where n=number of models')
        for i in range(n_model):
            row = contact[i]
            if not isinstance(row, list):
                raise TypeError('The contact matrix must be a list of lists')
            if len(row) != n_model:
                raise ValueError('The contact matrix must be an nxn matrix, where n=number of models')
            if np.abs(contact[i][i] - 1.) > 0.001:
                raise ValueError('The diagonal terms in the contact matrix should be 1.')

        self.contact_matrix = contact

    def reset_cross_transmission(self):
        self.infected_name = None
        self.susceptible_name = None
        self.total_name = None
        self.contagious_name = None
        self.alpha_name = None
        self.contact_matrix = None

    def set_t0(self, year=2020, month=1, day=1):
        """Define the date of the first element in the history lists
        """
        self.t0 = date(year, month, day)
        for model_name in self.models:
            model = self.models[model_name]
            model.set_t0(year, month, day)

    def set_time_step(self, time_step):

        if not isinstance(time_step, float):
            raise TypeError('Model (' + self.name +
                            ' time_step must be a float')
        if not time_step > 0.:
            raise ValueError('Model (' + self.name +
                             ' time_step must be positive')
        self.__time_step = time_step

    def get_time_step(self):
        """
        Get the length of time corresponding to one step in the evolution

        Returns
        -------
        real value for one time step

        """
        return self.__time_step

    def reset(self):
        """
        Reinitialize all populations, eliminating history and future, and
        setting initial values. Also reset transitioned parameters
        """
        for model_name in self.models:
            model = self.models[model_name]
            model.reset()
        self.__reset()

    def __reset(self):
        """
        Reinitialize all populations, eliminating history and future, and
        setting initial values. Also reset transitioned parameters
        """
        for key in self.populations:
            self.populations[key].reset()

        self.boot_needed = True

    def boot(self, expectations=True):

        # start from a clean slate of populations set to their initial values
        self.__reset()

        # need to set this to False before evolving expectations
        self.boot_needed = False

        for model_name in self.models:
            model = self.models[model_name]
            model.boot(expectations)

    def evolve_expectations(self, n_step):
        """
        First check if a boot is needed: this is required if
        reset_populations is called

        Take n_steps of
        0) check if it is time to apply a transition in any of the models
        1) updating internal futures in all the models
        2) update future due to infection from other models
        3) appending new value to history and remove current time step in future array
        """
        if self.boot_needed:
            self.boot(expectations=True)

        # turn off the infection cycle in each model, as the ensemble does this.
        self.__disable_infections()

        for step in range(n_step):
            self.do_transitions(step, expectations=False)

            # calculate influence between models
            self.__cross_model_transmission(expectations=True)

            # calculate future expectations
            self.calculate_future(expectations=True)

            # make one time step
            self.do_time_step(expectations=True)

        # add up the histories and put in local population objects
        self.__combine_histories()

        # turn back on the infection cycles in each model, since they need them to boot
        self.__enable_infections()

    def generate_data(self, n_step):
        """
        Produce data
        """
        if self.boot_needed:
            self.boot(expectations=False)

        # turn off the infection cycle in each model, as the ensemble does this.
        self.__disable_infections()

        for step in range(n_step):
            self.do_transitions(step, expectations=False)

            # calculate influence between models
            self.__cross_model_transmission(expectations=False)

            # calculate future data
            self.calculate_future(expectations=False)

            # make one time step
            self.do_time_step(expectations=False)

        # add up the histories and put in local population objects
        self.__combine_histories()

        # turn back on the infection cycles in each model, since they need them to boot
        self.__enable_infections()

    def do_transitions(self, step, expectations=True):
        for model_name in self.models:
            model = self.models[model_name]
            model.do_transitions(step, expectations)

    def calculate_future(self, expectations=True):
        for model_name in self.models:
            model = self.models[model_name]
            model.calculate_future(expectations)

    def __disable_infections(self):
        for model_name in self.models:
            model = self.models[model_name]
            model.connectors[self.infection_cycle_name].to_population = self.null_pop

    def __enable_infections(self):
        for model_name in self.models:
            model = self.models[model_name]
            model.connectors[self.infection_cycle_name].to_population = \
                model.populations[self.infected_name]

    def __cross_model_transmission(self, expectations=True):
        n = len(self.model_list)

        for ia in range(n):
            model_name_A = self.model_list[ia]
            model_A = self.models[model_name_A]
            sum_denom = 0.
            sum_numer = 0.
            for ib in range(n):
                model_name_B = self.model_list[ib]
                model_B = self.models[model_name_B]

                sum_denom += self.contact_matrix[ia][ib] * model_B.populations[self.total_name].history[-1]
                term1 = np.sqrt(model_A.parameters[self.alpha_name].get_value() *
                                model_B.parameters[self.alpha_name].get_value())
                term2 = model_A.populations[self.susceptible_name].history[-1]
                term3 = model_B.populations[self.contagious_name].history[-1]
                sum_numer += term1 * term2 * term3 * self.contact_matrix[ia][ib]
            new_infections = sum_numer / sum_denom

            if expectations:
                model_A.populations[self.infected_name].update_future_fast(new_infections)
            else:
                if self.__distribution == 'poisson':
                    n = stats.poisson.rvs(new_infections)
                    model_A.populations[self.infected_name].update_future_fast(n)
                else:
                    p = self.__nbinom_par.get_value()
                    if p < 0.001:
                        p = 0.001
                    if p > 0.999:
                        p = 0.999
                    r = new_infections * p / (1. - p)
                    n = stats.nbinom.rvs(r, p)
                    model_A.populations[self.infected_name].update_future_fast(n)

    def do_time_step(self, expectations=True):
        for model_name in self.models:
            model = self.models[model_name]
            model.do_time_step(expectations)

    def __combine_histories(self):
        # add the histories of all models together
        the_models = []
        for model_name in self.models:
            model = self.models[model_name]
            the_models.append(model)

        for pop_name in self.populations:
            ens_pop = self.populations[pop_name]
            ens_history = [0] * len(the_models[0].populations[pop_name].history)
            for the_model in the_models:
                mod_history = the_model.populations[pop_name].history
                for i in range(len(ens_history)):
                    ens_history[i] += mod_history[i]
            ens_pop.history = ens_history

    def add_transition(self, transition):
        raise NotImplementedError()

    def add_connector(self, connector, before_connector=None, after_connector=None):
        """
        Overwrite this method since it is important to not change the populations in the reference model
        """
        raise NotImplementedError()

    def remove_connector(self, connector):
        raise NotImplementedError()

    def update_lists(self):
        raise NotImplementedError()

    def save_file(self, filename):
        """
        Save a copy of the current ensemble to a file. The ensemble can be restored
        at a late time using Ensemble.open_file(filename).

        Before saving, the model futures are removed, and the initial_values
        are set to the current values for all parameters.

        If no extension is provided, the default extension, .pypm is added.

        Parameters
        ----------
        filename : str
            name of file to save model

        Returns
        -------
        None.

        """
        if not isinstance(filename, str):
            raise TypeError('Error saving file. ' +
                            ': filename argument must be a str')

        fullname = filename
        if '.' not in filename:
            fullname = filename + '.pypm_e'

        #       eliminate histories, undo transitioned parameter adjustments

        self.reset()

        #       set initial parameter values to their current values to define a new
        #       reference

        for model_name in self.models:
            model = self.models[model_name]
            for par_name in model.parameters:
                par = model.parameters[par_name]
                par.new_initial_value()

        with open(fullname, 'wb') as f:
            pickle.dump(self, f, pickle.HIGHEST_PROTOCOL)

    @classmethod
    def open_file(cls, filename):
        """
        Restore an ensemble that was saved to a file using Ensemble.save_file(filename)

        Parameters
        ----------
        filename : str
            name of existing model file to open

        Returns
        -------
        Ensemble
            The ensemble object saved in the file

        """
        if not isinstance(filename, str):
            raise TypeError('Error opening file. ' +
                            ': filename argument must be a str')

        fullname = filename
        if '.' not in filename:
            fullname = filename + '.pypm_e'

        with open(fullname, 'rb') as f:
            return pickle.load(f)
