"""
AgeModeller:
 - For Germany, it is essential to use age group information to forecast deaths (starting late 2020)
 - The age modeller provides automated ways to build simple models of recent history to be combined to
   form a useful model for each state
 - Also produces interval data so that the generated state model file can be used by the forecast_hub code

"""

from pathlib import Path
import numpy as np
import pandas as pd
import requests
import datetime
from datetime import timedelta
import copy

from pypmca import Model
from pypmca.analysis.Optimizer import Optimizer

class AgeModeller:
    """ Produce age group models automatically
    """

    def __init__(self, t_0, data_folder, hub, forecast_date):
        # t_0: the first day of data to be analyzed
        # data_folder: location of data
        # hub = "USA" or "Germany" for now
        # forecast_date: intended date for forecast (typically upcoming Sunday): datetime.date object
        self.t_0 = t_0
        self.data_folder = data_folder
        self.hub = hub
        self.forecast_date = forecast_date
        self.data_description = None
        self.pd_dict = self.get_data()

    def get_data(self):

        success = True
        pd_dict = {}
        try:
            data_desc_resp = requests.get('http://data.ipypm.ca/get_data_desc/' + self.data_folder)
        except requests.exceptions.RequestException as error:
            print('Error retrieving data description over network:')
            success = False
        if success:
            data_description = data_desc_resp.json()
            data_description['folder'] = self.data_folder
            print(data_description['description'])
            self.data_description = data_description

            # load the data into a panda dictionary
            for filename in data_description['files']:
                path = self.data_folder + '/' + filename
                success = True
                try:
                    csv_resp = requests.get('http://data.ipypm.ca/get_csv/' + path, stream=True)
                except requests.exceptions.RequestException as error:
                    print('Error retrieving data over network:')
                    print()
                    print(error)
                    success = False

                if success:
                    pd_dict[filename] = pd.read_csv(csv_resp.raw)

        return pd_dict

    def get_age_group_data(self,id):
        data_dict = None
        if id in self.data_description['regional_data']:
            data_dict = {}
            for pop in self.data_description['regional_data'][id]:
                filename = self.data_description['regional_data'][id][pop]['total']['filename']
                header = self.data_description['regional_data'][id][pop]['total']['header']
                ds = self.data_description['files'][filename]['date start']
                date_start = datetime.date(ds[0],ds[1],ds[2])
                days_skip = (self.t_0 - date_start).days
                data_dict[pop] = self.pd_dict[filename][header].fillna(0).values[days_skip:]
        return data_dict

    def fit_reported(self, model, data, trans_date_guess, n_rep, verbose):
        model.set_t0(self.t_0.year,self.t_0.month,self.t_0.day)
        trans_day_guess = (trans_date_guess - self.t_0).days

        model.parameters['trans_rate_1_time'].set_value(trans_day_guess)
        # estimate the number who are contagious using the first 8 days of data
        reported_8 = data['reported'][7] - data['reported'][0]
        cont_8 = reported_8
        model.parameters['cont_0'].set_value(cont_8)
        model.parameters['cont_0'].set_max(4*cont_8)
        model.boot_pars['boot_value'] = cont_8/50.

        # do fit of alpha_0, alpha_1, cont_0, trans_rate_1_time
        for par_name in ['alpha_0', 'alpha_1', 'cont_0']:
            par = model.parameters[par_name]
            par.set_variable(None, None)

        # find reasonable values for the parameters
        start_day = 0
        end_day = len(data['reported']) - 1
        optimizer = Optimizer(model, 'total reported', data['reported'], [start_day, end_day], cumul_reset=True)
        popt, pcov = optimizer.fit()

        for par_name in ['alpha_0', 'alpha_1', 'cont_0']:
            par = model.parameters[par_name]
            par.new_initial_value()

        par = model.parameters['trans_rate_1_time']
        par.set_variable(None, None)

        par.set_min(trans_day_guess)
        par.set_max(trans_day_guess+1)

        scan_dict = optimizer.i_fit()

        direction = +1
        min_chi2 = scan_dict['chi2_list'][1]
        trans_day = trans_day_guess + 1
        if scan_dict['chi2_list'][0] < scan_dict['chi2_list'][1]:
            direction = -1
            min_chi2 = scan_dict['chi2_list'][0]
            trans_day = trans_day_guess

        min_found = False
        while not min_found:
            trans_day_try = trans_day + direction
            par.set_min(trans_day_try)
            par.set_max(trans_day_try)

            scan_dict = optimizer.i_fit()
            if scan_dict['chi2_list'][0] < min_chi2:
                min_chi2 = scan_dict['chi2_list'][0]
                trans_day = trans_day_try
            else:
                min_found = True

        model.parameters['trans_rate_1_time'].set_value(trans_day)

        par = model.parameters['trans_rate_1_time']
        par.set_fixed()
        popt, pcov = optimizer.fit()

        # find and assign uncertainties
        self.set_std_estimators(model, optimizer, n_rep, verbose)

        # Include the uncertainty from unknown transition day
        delta_days = [-2, +2]
        mod_alphas = []
        for delta_day in delta_days:
            temp_model = copy.deepcopy(model)
            new_date = temp_model.transitions['trans_rate_1'].transition_time.get_value() + delta_day
            temp_model.transitions['trans_rate_1'].transition_time.set_value(new_date)
            temp_optimizer = Optimizer(temp_model, 'total reported', data['reported'], [start_day, end_day], cumul_reset=True)
            popt, pcov = temp_optimizer.fit()
            mod_alpha = temp_model.parameters['alpha_1'].get_value()
            mod_alphas.append(mod_alpha)

        if verbose:
            print('Transition: trans_rate_1 on day',model.transitions['trans_rate_1'].transition_time.get_value())
            print('alpha values: \n  nom = {0:0.4f} +/- {1:0.4f} \n  {2:+d} days = {3:0.4f} \n  {4:+d} days = {5:0.4f}'.format(
                model.parameters['alpha_1'].get_value(),
                model.parameters['alpha_1'].std_estimator,
                delta_days[0], mod_alphas[0], delta_days[1], mod_alphas[1]
            ))

        # while the following should be divided by 2, leave as is to account for larger delta_day possibility
        mod_alphas_std = np.abs(mod_alphas[0] - mod_alphas[1])
        current_std = model.parameters['alpha_1'].std_estimator
        new_std = np.sqrt(current_std ** 2 + mod_alphas_std ** 2)
        model.parameters['alpha_1'].std_estimator = new_std

        if verbose:
            print('Additional error included in final alpha: {0:0.4f}'.format(mod_alphas_std))

    def fit_deaths(self, model, data, n_rep, verbose, start_deaths):
        # fit death parameters
        for par_name in ['alpha_0', 'alpha_1', 'cont_0']:
            par = model.parameters[par_name]
            par.set_fixed()
        # do fit of recover_frac
        for par_name in ['recover_frac']:
            par = model.parameters[par_name]
            par.set_variable(None, None)
        start_day = start_deaths
        end_day = len(data['deaths']) - 1
        optimizer = Optimizer(model, 'total deaths', data['deaths'], [start_day, end_day], cumul_reset=True)
        popt, pcov = optimizer.fit()

        # find and assign uncertainties
        self.set_std_estimators(model, optimizer, n_rep, verbose)

    def set_std_estimators(self, model, optimizer, n_rep, verbose):
        # Estimate the properties of the estimators, by making many several simulated samples to find the bias and covariance
        # assign the standard deviations to the model parameters
        optimizer.calc_chi2f = True
        optimizer.calc_chi2s = False
        optimizer.calc_sim_gof(n_rep)

        # Save the standard deviations of the estimators
        results = []
        for result in optimizer.opt_lists['opt']:
            results.append(result)
        transposed = np.array(results).T

        if verbose:
            for i, par_name in enumerate(optimizer.variable_names):
                std = np.std(transposed[i])
                model.parameters[par_name].std_estimator = std
                val = model.parameters[par_name].get_value()
                print('  ' + par_name, '{0:0.4f} +/- {1:0.4f}'.format(val,std))

