# -*- coding: utf-8 -*-
"""
Optimizer: Given a model and data, find point estimates for the variable parameters in the model,
along with the covariance (optional)

This is done using the scipy.optimize.curve_fit method

The full_population_name has the prefix daily or total. If daily, then the
differences need to be calculated

@author: karlen
"""
import copy
from scipy.optimize import curve_fit
from scipy import stats
import numpy as np

class Optimizer:
    """ Optimizer: point estimator and covariance
    """

    def __init__(self, model, full_population_name, data, data_range):
        self.model = model
        self.full_population_name = full_population_name
        self.population_name = full_population_name[6:]
        self.population_type = full_population_name[:5]
        self.data = data
        self.data_range = data_range
        self.variable_names = []  # float ype variable parameters
        self.variable_initial_values = None
        self.i_variable_names = []  # float ype variable parameters
        self.i_variable_initial_values = None
        self.chi2d = 0.
        self.chi2m = 0.
        self.chi2s = 0.
        self.accept_fraction = 0.
        self.auto_cov = None

    def func_setup(self):
        # find the variable parameters and return the bounds
        # note: scipy.curve_fit only works for floats
        self.variable_names = []
        self.variable_initial_values = {}
        par_0 = []
        bound_low = []
        bound_high = []
        for par_name in self.model.parameters:
            par = self.model.parameters[par_name]
            if par.get_status() == 'variable':
                if par.parameter_type == 'float':
                    self.variable_names.append(par_name)
                    self.variable_initial_values[par_name] = par.get_value()
                    par_0.append(par.get_value())
                    bound_low.append(par.get_min())
                    bound_high.append(par.get_max())
        return par_0, (bound_low, bound_high)
    
    def i_func_setup(self):
        self.i_variable_names = []
        self.i_variable_initial_values = {}
        par_0 = []
        bound_low = []
        bound_high = []
        for par_name in self.model.parameters:
            par = self.model.parameters[par_name]
            if par.get_status() == 'variable':
                if par.parameter_type == 'int':
                    self.i_variable_names.append(par_name)
                    self.i_variable_initial_values[par_name] = par.get_value()
                    par_0.append(par.get_value())
                    bound_low.append(par.get_min())
                    bound_high.append(par.get_max())
        return par_0, (bound_low, bound_high)

    def reset_variables(self):
        for var_name in self.variable_names:
            self.model.parameters[var_name].set_value(self.variable_initial_values[var_name])
            
    def reset_i_variables(self):
        for var_name in self.i_variable_names:
            self.model.parameters[var_name].set_value(self.i_variable_initial_values[var_name])
            
    def delta(self, cumul):
        diff = []
        for i in range(1,len(cumul)):
            diff.append(cumul[i] - cumul[i-1])
        return diff

    def i_fit(self):
        """ scan the first integer variable parameter over the range provided
            iterating over different integers by hand
            Returns a dict if a scan was done over an integer,
            describing the scan
        """
        
        i_par0, i_bounds = self.i_func_setup()
        if len(i_par0) == 0:
            return None

        par_name = self.i_variable_names[0]
        par = self.model.parameters[par_name]
        scan_dict = {'name': par_name}
        val_list = []
        chi2_list = []
        
        for i in range(par.get_min(), par.get_max()+1):
            par.set_value(i)
            self.fit()
            val_list.append(i)
            chi2_list.append(self.chi2d)
            self.reset_variables()
 
        arg = np.argmin(chi2_list)
        par.set_value(val_list[arg])

        scan_dict['val_list'] = val_list
        scan_dict['chi2_list'] = chi2_list
        return scan_dict        

    def fit(self, with_cov=False):
        """ work out point estimate. Note that this does not seem to work
            well if the auto_covariance is used - perhaps bug in curve_fit
        """

        def func(x, *params):
            i = -1
            for par_val in params:
                i += 1
                par_name = self.variable_names[i]
                self.model.parameters[par_name].set_value(par_val)

            self.model.reset()
            self.model.evolve_expectations(self.data_range[1])
            pop = self.model.populations[self.population_name]
            func_values = []
            if self.population_type == 'total':
                for xi in x:
                    i = int(round(xi))
                    func_values.append(pop.history[i])
            else:
                diff = self.delta(pop.history)
                for xi in x:
                    i = int(round(xi))
                    func_values.append(diff[i])
            return np.array(func_values)

        par_0, bounds = self.func_setup()
        xdata = np.arange(self.data_range[0], self.data_range[1], 1)

        popt = None
        pcov = None
        if not with_cov:
            popt, pcov = curve_fit(func, xdata, self.data[self.data_range[0]:self.data_range[1]],
                                   p0=par_0, bounds=bounds)
        else:
            popt, pcov = curve_fit(func, xdata, self.data[self.data_range[0]:self.data_range[1]],
                                   p0=par_0, bounds=bounds, sigma=self.auto_cov)

        i = -1
        for par_val in popt:
            i += 1
            par_name = self.variable_names[i]
            self.model.parameters[par_name].set_value(par_val)

        self.model.reset()
        self.model.evolve_expectations(self.data_range[1])
        self.chi2d = 0.
        for x in xdata:
            resid = self.data[x] - self.model.populations[self.population_name].history[x]
            self.chi2d += resid**2/self.model.populations[self.population_name].history[x]

        return popt, pcov

    def calc_auto_covariance(self, n_rep=100):
        """ Calculate the autocovariance matrix. Do that by generating simulated datasets.
            At the same time, investigate the chi^2 for cases where the last simulated point
            is close to the expectation (which would be the case for fitting)
            
            STILL TO DO: deal with daily data! Not hard. see above ***
        """
        # Make a copy first

        sim_model = copy.deepcopy(self.model)
        sim_population = sim_model.populations[self.population_name]
        ref_population = self.model.populations[self.population_name]

        sigma_last = np.sqrt(ref_population.history[self.data_range[1]])
        chi2_list = []

        xdata = np.arange(self.data_range[0], self.data_range[1], 1)
        n_p = len(xdata)
        cov = [[0. for i in range(n_p)] for j in range(n_p)]
        for irep in range(n_rep):
            # produce simulated data -> sim_population histories
            sim_model.reset()
            sim_model.generate_data(self.data_range[1])

            # calculate residuals (data-expectation) for cov and chi^2
            resid = []
            chi2 = 0.
            for x in xdata:
                delta = sim_population.history[x]-ref_population.history[x]
                resid.append(delta)
                chi2 += delta**2/ref_population.history[x]

            if abs(resid[-1]) < sigma_last/4.:
                chi2_list.append(chi2)

            for i in range(n_p):
                for j in range(n_p):
                    cov[i][j] += resid[i]*resid[j]/n_rep

        self.chi2m = np.mean(chi2_list)
        self.auto_cov = cov

    def calc_sim_gof(self, n_rep=100):
        """ Calculate the goodness of fit statistic for simulated data.
            The average of this is chi2s
        """

        sim_model = copy.deepcopy(self.model)
        sim_population = sim_model.populations[self.population_name]
        ref_population = self.model.populations[self.population_name]

        xdata = np.arange(self.data_range[0], self.data_range[1], 1)
        n_p = len(xdata)
        zeros = [0. for i in range(n_p)]
        lpdf_zero = stats.multivariate_normal.logpdf(zeros, cov=self.auto_cov)
        chi2 = []
        for i in range(n_rep):
            sim_model.reset()
            sim_model.generate_data(self.data_range[1])
    
            resid = []
            for x in xdata:
                resid.append(sim_population.history[x]-ref_population.history[x])

            lpdf = stats.multivariate_normal.logpdf(resid, cov=self.auto_cov)
            chi2.append(2.*(lpdf_zero-lpdf))

        self.chi2s = np.mean(chi2)

    def mcmc(self, n_dof, n_MCMC):
        """ Make a MCMC chain (n_MCMC points) assuming n_dof defines gof statistic for data. 
            Method returns chain in array of dictionaries
        """

        var_names = []
        var_parameters = {}
        var_values = {}
        prior_function = {}
        hypercube = {}
        mean = {}
        half_width = {}
        sigma = {}

        sim_model = copy.deepcopy(self.model)
        sim_population = sim_model.populations[self.population_name]

        for par_name in sim_model.parameters:
            par = sim_model.parameters[par_name]
            if par.get_status() == 'variable':
                var_names.append(par_name)
                var_parameters[par_name] = par
                var_values[par_name] = par.get_value()
                prior_function[par_name] = par.prior_function
                hypercube[par_name] = par.mcmc_step
                mean[par_name] = par.prior_parameters['mean']
                if par.prior_function == 'uniform':
                    half_width[par_name] = par.prior_parameters['half_width']
                else:
                    sigma[par_name] = par.prior_parameters['sigma']

        xdata = np.arange(self.data_range[0], self.data_range[1], 1)
        zeros = [0. for i in range(len(xdata))]
        lpdf_zero = stats.multivariate_normal.logpdf(zeros, cov=self.auto_cov)
        # calculate log of posterior probability parameters
        def logP():
            # check if selected parameters within bounds (uniform priors)
            for var_name in var_names:
                if prior_function[var_name] == 'uniform':
                    par = var_parameters[var_name]
                    if np.abs(par.get_value() - mean[var_name]) > half_width[var_name]:
                        return -np.inf

            sim_model.reset()
            sim_model.evolve_expectations(self.data_range[1])

            resid = []
            for x in xdata:
                resid.append(self.data[x]-sim_population.history[x])

            lpdf = stats.multivariate_normal.logpdf(resid, cov=self.auto_cov)
            chi2 = 2.*(lpdf_zero-lpdf)

            ln_P = stats.chi2.logpdf(chi2, n_dof)

            #Deal with any normal priors
            for var_name in var_names:
                if prior_function[var_name] == 'normal':
                    par = var_parameters[var_name]
                    dev = (par.get_value() - mean[var_name])/sigma[var_name]
                    ln_P -= dev*dev/2.

            return ln_P

        chain = []

        n_accept = 0
        lp_curr = logP()
        for i in range(n_MCMC):
            new_var_values = var_values.copy()
            for var_name in var_names:
                new_var_values[var_name] += hypercube[var_name]*(1.-2.*stats.uniform.rvs())
                par = var_parameters[var_name]
                par.set_value(new_var_values[var_name])

            lp_new = logP()
            if lp_new != -np.inf:
                del_lp = lp_new - lp_curr
                if del_lp > -30:
                    if del_lp > 0 or (stats.uniform.rvs() < np.exp(del_lp)):
                        var_values = new_var_values
                        lp_curr = lp_new
                        n_accept += 1
            chain.append(var_values)

        self.accept_fraction = 1.*n_accept/n_MCMC
        return chain
