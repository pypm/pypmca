# -*- coding: utf-8 -*-
"""
Ensemble: A Model class whose behaviour is defined by an ensemble of models

An ensemble allows for categorization of populations by age, risk, or other factors. It can also
be used to combine many models that may or may not benearly independent
(such as separate provinces to make a Canada wide model).

Each category has its own model, and one model can influence the growth within other models.
The ensemble sums the histories of all its models to represent the evolution of
the entire system.

When mixing models having different growth rates, achieving the desired initial condition at t_0 through the
boot process is challenging. The boot process starts with a small number in each model's boot_population and
the boot ends when the ensemble exceeds its goal. It is possible that the relative sizes for each model
would differ significantly from the desired proportions, at t=0, and after scaling each model for the t_0 condition,
the state would be far from a steady state. A second boot can be done (by adjusting the small numbers in each model's
boot_population) according to the outcome from the first boot - ie. doing an iterative boot.
This issue is not as serious if, at least initially, the growth behaviours of the different groups are similar.

Owing to the complexity of ensembles made from mixtures of very different growth rate sub groups,
such situations should be treated with care until such behaviour has been thoroughly tested!

@author: karlen
"""
from datetime import date
import numpy as np
from scipy import stats
import copy
from pathlib import Path
import pickle
from pypmca import Model, Parameter, Population


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
        in_model = None
        if isinstance(model, Model):
            in_model = model
        elif isinstance(model, str):
            in_model = Model.open_file(model)
        else:
            raise TypeError('Error creating ensemble. The model argument needs to be a model object ' +
                            'or have a .pypm filename')

        # to avoid the possibility that the reference model is subsequently modified, make a copy
        # that will not be as exposed to tampering
        self.model = copy.deepcopy(in_model)

        # point to the population objects in self.model to save the histories
        # this allows for consistent access to history between models and ensembles
        for pop_name in self.model.populations:
            pop = self.model.populations[pop_name]
            self.populations[pop_name] = pop

        # the parameters will be set once the list of models are read in
        self.parameters = {}

        # the boot parameters point to the ones copied from the reference
        self.boot_pars = self.model.boot_pars

        # the important scaling parameter for the ensemble gets defined here
        # this is forced to be a parameter when each model is created and
        # checked when the models are included in the ensemble

        boot_pop = self.populations[self.boot_pars['boot_population']]
        boot_pop_ip = boot_pop.initial_value
        self.parameters[boot_pop_ip.name] = boot_pop_ip
        # When this scaling parameter is changed, all the model scaling parameters
        # must also be scaled accordingly.
        boot_pop_ip.set_must_update(self)

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
        self.contact = None
        self.contact_matrix = None
        self.contact_type = ''
        self.null_pop = Population('null', 0.)
        self.total_goal = None
        self.boot_achieved = {}
        self.max_scale_factor = None

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

    def update(self):
        """
        This is called when the ensemble scaling parameter is changed. It is necessary
        to update the model scaling parameters accordingly.
        """
        ens_scale = 0.
        for model_name in self.models:
            model = self.models[model_name]
            boot_pop = model.populations[model.boot_pars['boot_population']]
            ens_scale += boot_pop.initial_value.get_value()
        if ens_scale > 0.:
            ens_boot_pop = self.populations[self.boot_pars['boot_population']]
            ratio = ens_boot_pop.initial_value.get_value()/ens_scale
            for model_name in self.models:
                model = self.models[model_name]
                boot_pop = model.populations[model.boot_pars['boot_population']]
                new_scale = boot_pop.initial_value.get_value() * ratio
                boot_pop.initial_value.set_value(new_scale)

            self.boot_needed = True

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
        if model.boot_pars['boot_population'] != self.model.boot_pars['boot_population']:
            raise ValueError('Error: Model ' + model.name + ' boot population name (' + model.boot_pars['boot_population'] +
            ') does not match name in reference boot population (' + self.model.boot_pars['boot_population'])
        boot_pop_ip = model.populations[model.boot_pars['boot_population']].initial_value
        if not isinstance(boot_pop_ip, Parameter):
            raise ValueError(
                'Error: Model ' + model.name + ' boot population (' + model.boot_pars['boot_population'] +
                ') is not initialized by a parameter')

        self.model_list.append(model.name)
        self.models[model.name] = model
        self.reset_cross_transmission()

        # include parameters
        model_app = '#' + str(len(self.model_list) - 1)
        for par_name in model.parameters:
            new_name = par_name + model_app
            self.parameters[new_name] = model.parameters[par_name]

    def define_cross_transmission(self, infection_cycle_name: str, infected_name: str,
                                  susceptible_name: str, total_name: str,
                                  contagious_name: str, alpha_name: str,
                                  contact_type: str = 'independent', contact: list = None):
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

        The contact matrix can be specified in a number of ways, defined by contact_type.
            - 'diagonal','independent': all diagonal terms are 1, non-diagonal are 0.
            This corresponds to the situation where all models develop independently
            - 'equality': all matrix elements are 1.
            This corresponds to the situation where members are blind to the catagorizations
            - 'fixed': an arbitrary contact_matrix is provided as a list of list of floats
            Diagonal elements are inforced to be 1. An ValueError is raised if not.
            - 'simple': all off-diagonal elements are equal and specified by a single
            Parameter object provided in list passed by the contact argument
            - 'symmetric': off-diagonal elements passed by n(n-1)/2 Parameter objects
            provided in the list passed by the contact argument - f12,f13,...,f1n,f23,f24...
        Parameters are included in the ensemble parameter list for ease in adjusting the mixing.

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

        self.contact_type = contact_type
        n_model = len(self.model_list)
        if contact_type == 'fixed':
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
        elif contact_type == 'simple':
            if not isinstance(contact, list):
                raise TypeError('The contact argument must be a list of Parameters')
            if len(contact) != 1:
                raise ValueError('The contact argument must be a list of length 1, ')
            for cont in contact:
                if not isinstance(cont, Parameter):
                    raise TypeError('The contact argument must be a list containing a Parameter object')
            self.contact = contact[0]
            self.__add_to_parameter_list(contact[0])
        elif contact_type == 'symmetric':
            if not isinstance(contact, list):
                raise TypeError('The contact argument must be a list of Parameters')
            if len(contact) != n_model * (n_model - 1) / 2:
                raise ValueError('The contact argument must be a list with ' +
                                 str(n_model * (n_model - 1) / 2) + ' parameters')
            for cont in contact:
                if not isinstance(cont, Parameter):
                    raise TypeError('The contact argument must be a list of Parameter objects')
                self.__add_to_parameter_list(cont)
            self.contact = contact
        elif contact_type not in ['diagonal', 'independent', 'equality']:
            raise ValueError('Contact matrix type not recognized.')

    def __add_to_parameter_list(self, obj):
        # check that the name of the parameter is unique
        if str(obj) in self.parameters:
            par = self.parameters[str(obj)]
            if obj != par:
                raise ValueError('Two parameters share the same name: ' + str(obj))
        else:
            self.parameters[str(obj)] = obj

    def reset_cross_transmission(self):
        self.infected_name = None
        self.susceptible_name = None
        self.total_name = None
        self.contagious_name = None
        self.alpha_name = None
        self.contact_matrix = None
        self.contact_type = ''

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

    def do_boot(self, expectations=True):
        if self.contact_type in ['diagonal', 'independent']:
            # in this case, boot each model separately
            # each will boot to its own goal

            # turn on the infection cycle in each model, in case it was disabled previously
            self.__enable_infections()

            # start from a clean slate of populations set to their initial values
            self.__reset()

            # need to set this to False before evolving expectations
            self.boot_needed = False

            for model_name in self.models:
                model = self.models[model_name]
                model.boot(expectations)

            # prior to full calculation, setup contact matrix etc.

            # turn off the infection cycle in each model, as the ensemble does this.
            self.__disable_infections()

            # setup the contact matrix
            self.__setup_contact_matrix()

        else:

            # try booting once, then a second time in case the scaling factors are too large
            self.boot_achieved = {}
            self.boot(expectations)

            # second boot in case scale factors are large
            # always do this to avoid discontinuities
            self.boot(expectations)

    def boot(self, expectations=True):
        """ Using information in boot_pars, evolve the system to a point where
        the target population in the ensemble is reached or exceeded:
            - expectations: if True, after the boot completes, the
            history will begin with an expectation,
            otherwise, a Poisson number is drawn using that expectation.
        Note: check is in place so that if boot requirement is not made within
        10k steps, it fails (rather than stay in infinite loop)

        The boot process itself may need to be iterated. If the relative sizes of the
        target population differs by a lot from the intended values, then the scaled
        populations after the boot may not be consistent with a steady state solution.
        """
        if len(self.boot_pars) == 0:
            raise RuntimeError('Model boot parameters not set.')

        # start from a clean slate of populations set to their initial values
        self.reset()

        # need to set this to False before evolving expectations
        self.boot_needed = False

        # save the goal_values for each model
        goal_values = {}
        goal_value = 0
        boot_value = 0
        for model_name in self.models:
            model = self.models[model_name]
            model_boot_pop = model.populations[model.boot_pars['boot_population']]
            goal_values[model_name] = model_boot_pop.history[-1]
            goal_value += model_boot_pop.history[-1]
            boot_value += model.boot_pars['boot_value']

        # the ensemble starts with a small population size in the boot_pop (defined by reference model)
        # the goal is to have this population grow to the sum of the sizes in each model
        boot_pop = self.populations[self.boot_pars['boot_population']]
        boot_pop.history[-1] = self.boot_pars['boot_value']
        # this small number needs to be distributed accordingly to the models
        # If this is the first boot, the best guess is to distribute the small number in proportion to
        # the goal numbers for each model. On a subsequent iteration, the scaling from the
        # results of the previous boot will provide a better solution.
        if len(self.boot_achieved) == 0:
            # First boot:
            self.total_goal = 0.
            for model_name in self.models:
                model = self.models[model_name]
                model_boot_pop = model.populations[model.boot_pars['boot_population']]
                self.total_goal += model_boot_pop.history[-1]
            if self.total_goal == 0.:
                raise ValueError('Targets in models are zero. Unable to boot.')
            boot_scale = 1.*boot_value/self.total_goal
            for model_name in self.models:
                model = self.models[model_name]
                model_boot_pop = model.populations[model.boot_pars['boot_population']]
                model_goal = model_boot_pop.history[-1]
                model_boot_pop.history[-1] = model_goal * boot_scale
        else:
            # Second boot:
            # the first boot derived fractions for each model based on goal/total
            # eg. desired population 10 and 30, the fractions are 0.25 and 0.75.
            # If that achieved relative populations after boot significantly different
            # eg.  2 and 38, then the starting fractions need to be adjusted to get
            # closer to the goal propotions after boot.
            # The original proportioning (goal/total) is revised by scaling each by goal/achieved
            # The new scaling is therefore goal^2/(achieved*total). If achieved is close to goal
            # this is a small change to the proportioning.
            total_rev = 0.
            for model_name in self.models:
                model = self.models[model_name]
                model_boot_pop = model.populations[model.boot_pars['boot_population']]
                total_rev += (model_boot_pop.history[-1])**2/self.total_goal/self.boot_achieved[model_name]
            boot_scale = 1. * boot_value / total_rev
            for model_name in self.models:
                model = self.models[model_name]
                model_boot_pop = model.populations[model.boot_pars['boot_population']]
                model_goal = model_boot_pop.history[-1]
                model_boot_pop.history[-1] = model_goal**2/self.total_goal/self.boot_achieved[model_name] * boot_scale

        # exclusion_populations for each model
        # are reset to their original values after the
        # boot. Save them here.

        model_exclusions = {}
        for model_name in self.models:
            model = self.models[model_name]
            exclusions = model.boot_pars['exclusion_populations']
            exc_vals = {}
            for exc_pop_name in exclusions:
                exc_pop = model.populations[exc_pop_name]
                exc_vals[exc_pop_name] = exc_pop.history[-1]
            model_exclusions[model_name] = exc_vals

        # turn off the infection cycle in each model, as the ensemble does this.
        self.__disable_infections()

        # setup the contact matrix
        self.__setup_contact_matrix()

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
            raise ValueError('The boot process did not reach its target. Consider reducing the target or adjusting a parameter.')

        # save the achieved values: needed if there is a second boot
        # quality of boot is good if the scale factors are near 1
        self.max_scale_factor = 0.
        for model_name in self.models:
            model = self.models[model_name]
            model_boot_pop = model.populations[model.boot_pars['boot_population']]
            self.boot_achieved[model_name] = model_boot_pop.history[-1]
            scale_factor = model_boot_pop.history[-1]/goal_values[model_name]
            if scale_factor < 1.:
                scale_factor = 1./scale_factor
            if scale_factor > self.max_scale_factor:
                self.max_scale_factor = scale_factor

        ens_scale = goal_value / boot_pop.history[-1]

        # Erase models history until current. Scale all histories and futures.
        for model_name in self.models:
            model =self.models[model_name]
            exc_vals = model_exclusions[model_name]
            model_boot_pop = model.populations[model.boot_pars['boot_population']]
            scale = goal_values[model_name] / model_boot_pop.history[-1]
            for key in model.populations:
                pop = model.populations[key]
                pop.remove_history()
                if str(pop) not in exc_vals:
                    pop.scale_history(scale, expectations)
                else:
                    pop.history[-1] = exc_vals[str(pop)]
                pop.scale_future(scale, expectations)

        # restore goal_value for data production
        if not expectations:
            for model_name in self.models:
                model = self.models[model_name]
                model_boot_pop = model.populations[model.boot_pars['boot_population']]
                model_boot_pop.history[-1] = int(goal_values[model_name])

    def evolve_expectations(self, n_step, from_step=0):
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
            self.do_boot(expectations=True)

        for step in range(from_step, n_step+from_step):
            self.do_transitions(step, expectations=True)

            # calculate influence between models
            self.__cross_model_transmission(expectations=True)

            # calculate future expectations
            self.calculate_future(expectations=True)

            # make one time step
            self.do_time_step(expectations=True)

        # add up the histories and put in local population objects
        self.__combine_histories()

    def generate_data(self, n_step, from_step):
        """
        Produce data
        """
        if self.boot_needed:
            self.do_boot(expectations=False)

        for step in range(from_step, n_step+from_step):
            self.do_transitions(step, expectations=False)

            # calculate influence between models
            self.__cross_model_transmission(expectations=False)

            # calculate future data
            self.calculate_future(expectations=False)

            # make one time step
            self.do_time_step(expectations=False)

        # add up the histories and put in local population objects
        self.__combine_histories()

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

    def __setup_contact_matrix(self):
        n_model = len(self.model_list)
        if self.contact_type in ['diagonal', 'independent']:
            self.contact_matrix = []
            for i in range(n_model):
                row = [0.] * n_model
                row[i] = 1.
                self.contact_matrix.append(row)
        elif self.contact_type == 'equality':
            self.contact_matrix = []
            for i in range(n_model):
                row = [1.] * n_model
                self.contact_matrix.append(row)
        elif self.contact_type == 'simple':
            self.contact_matrix = []
            off_diag = self.contact.get_value()
            for i in range(n_model):
                row = [off_diag] * n_model
                row[i] = 1.
                self.contact_matrix.append(row)
        elif self.contact_type == 'symmetric':
            self.contact_matrix = [[0. for i in range(n_model) ] for i in range(n_model)]
            i_pnt = 0
            for i in range(n_model):
                self.contact_matrix[i][i] = 1.
                for j in range(i + 1, n_model):
                    self.contact_matrix[i][j] = self.contact[i_pnt].get_value()
                    self.contact_matrix[j][i] = self.contact[i_pnt].get_value()
                    i_pnt += 1

    def __cross_model_transmission(self, expectations=True):
        # check that everything is ready
        if not self.contact_type == 'fixed' and self.contact_matrix is None:
            raise RuntimeError('The cross-transmission information has not been provided yet.')

        n = len(self.model_list)

        for ia in range(n):
            model_name_A = self.model_list[ia]
            model_A = self.models[model_name_A]
            sum_denom = 0.
            sum_numer = 0.
            if self.contact_type in ['diagonal', 'independent']:
                sum_denom += model_A.populations[self.total_name].history[-1]
                term1 = model_A.parameters[self.alpha_name].get_value()
                term2 = model_A.populations[self.susceptible_name].history[-1]
                term3 = model_A.populations[self.contagious_name].history[-1]
                sum_numer += term1 * term2 * term3
            else:
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
                    n_poisson = stats.poisson.rvs(new_infections)
                    model_A.populations[self.infected_name].update_future_fast(n_poisson)
                else:
                    p = self.__nbinom_par.get_value()
                    if p < 0.001:
                        p = 0.001
                    if p > 0.999:
                        p = 0.999
                    r = new_infections * p / (1. - p)
                    n_binom = stats.nbinom.rvs(r, p)
                    model_A.populations[self.infected_name].update_future_fast(n_binom)

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

        If no extension is provided, the default extension, .pypm_e is added.

        Parameters
        ----------
        filename : str
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
            filepath = filepath.with_suffix('.pypm_e')


        #       eliminate histories, undo transitioned parameter adjustments

        self.reset()

        #       set initial parameter values to their current values to define a new
        #       reference

        for model_name in self.models:
            model = self.models[model_name]
            for par_name in model.parameters:
                par = model.parameters[par_name]
                par.new_initial_value()

        with open(filepath, 'wb') as f:
            pickle.dump(self, f, protocol=4)

    @classmethod
    def open_file(cls, filepath):
        """
        Restore an ensemble that was saved to a file using Ensemble.save_file(filename)

        Parameters
        ----------
        filepath : Path or str
            name of existing ensemble file to open

        Returns
        -------
        Ensemble
            The ensemble object saved in the file

        """

        try:
            filepath = Path(filepath).resolve()
        except:
            raise TypeError('Input arg could not be converted to a valid path: {}' +
                            '\n It must be a str or Path-like.'.format(filepath))

        if not filepath.exists():
            raise ValueError('Filepath does not exist: {}'.format(filepath))

        if len(filepath.suffix) < 2:
            filepath = filepath.with_suffix('.pypm_e')

        with open(filepath, 'rb') as f:
            return pickle.load(f)
