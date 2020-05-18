# -*- coding: utf-8 -*-
"""
Optimizer: Given a model and data, find point estimates for the variable parameters in the model,
along with the covariance (optional)

This is done using the scipy.optimize.curve_fit method

The full_population_name has the prefix daily or total. If daily, then the
differences need to be calculated

The MCMC method relies on calculating a likelihood. The approach was designed around using a cumulative
indicator (reported cases). This was motivated because this indicator suffers from large reporting noise
on a daily basis - with large negative correlations between neighbouring days. The effect is smaller on the
cumulative. Still, reporting noise should be included in simulation samples (it is an optional parameter in pypmca
population models).

Because cumulative data is being used, it is essential to include the auto-covariance in the likelihood.
The shape and normalization of the cumulative indicator are separated so that the auto-covariance
does not include the very large effect from overall normalization. The separation is achieved when
calculating chi^2, by scaling the reference model so that its last point coincides with the simulated/real data.
The contribution from normalization can be evaluated separately by using the values in the last point.
The autocovariance therefore does not include the last point: its rank is one less than the number of points being fit.

It is assumed that the model always has an initial state parameter that scales all final state population sizes.
In the reference models 1&2, this is the cont_0 parameter (number contagious at day 0).
Given that this parameter exists, we can ignore the
contribution to the normalization in evaluating the likelihood.
The simple scaling applied (to separate the shape from normalization) assumes that the initial (t_0)
value for the indicator is near zero.

The following steps are taken to define the likelihood for the MCMC analysis:
(1) The data is fit to find point estimates for the variable parameters. This is done in the analysis tab, and
must be done prior to starting the MCMC analysis. As a check, a naive goodness of fit (shape only) is calculated
and saved as self.chi2d. It is naive, in that it assumes that the data are outcomes of independent Poisson random variables.
The model with the parameters set to these point estimates is referred to as the reference model.

(2) The autocovariance is calculated by producing a set of many simulated datasets from this reference model.
The residuals from the reference model expectations are used to calculate the model auto-covariance. This may differ from the actual
autocovariance if reporting errors are not correctly included or if the model does not contain other important
sources of variance.

(3) An independent set of many simulated datasets from the reference model are produced and the goodness
of fit statistic for each simulated dataset is calculated using different methods to understand the importance of
autocorrelation and to detect any significant difference between noise in data and simulation:
  * chi2_m is the naive goodness of fit (shape only) calculated with respect to the reference model
  * chi2_n is the goodness of fit (normalization only): one degree of freedom - the last data value
  * chi2_f is the naive goodness of fit (shape only) calculated with respect to a reference model fitted to that data
  (this is time consuming and only calculated if self.calc_chi2_f is set to True)
  * chi2_s is the goodness of fit (shape only, with autocorrelation) calculated with respect to the reference model
The means and standard deviations of the gof distributions are saved (in self.chi2m, self.chi2m_sd etc) and the lists
of the values are also available for further study in self.chi2_lists dictionary with keys 'm', 'f', 's'

The chi2_n is calculated with only one degree of freedom (the last data point). This is far from being
well described by a Poisson random variable. Due to fluctuations early in an exponential growth period, the
fluctuations can be very large. Empirically it was found that the distribution for chi2_n/mean(chi2_n) does
follow a chi2 distribution for one degree of freedom. So the likelihood calculation includes, in addition to
the shape contribution, the normalization contribution as calculated by chi2_n/mean(chi2_n), adding one degree
of freedom.

Note: in checking this method with simulated samples, it was found that the chi2_s as calculated by
-2 log(delta likelihood), where likelihood is from the stats.multivariate_normal, follows the chi^2 distribution
much better when the definition is altered (empirically) such that chi2_s = - log(delta likelihood).
Also, the mean of the distribution is much closer to the expected number of degrees of freedom.
It is not understood why this is the case. Since the sampling distribution is better described
and it yields more conservative confidence intervals, the altered statistic is used.

(4) The likelihood is calculated using the probability density of the chi^2 distribution.
The number of degress of freedom can be specified. If the model included all sources of variation correctly,
the best estimate for dof would be the mean(chi2s). This value should be
increased to account for additional degrees of freedom that reflect the unaccounted for noise, coming from the study in part (3).

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
        self.chi2n = 0.
        self.chi2s = 0.
        self.chi2f = 0.
        self.chi2m_sd = 0.
        self.chi2n_sd = 0.
        self.chi2f_sd = 0.
        self.chi2s_sd = 0.
        self.accept_fraction = 0.
        self.auto_cov = None
        self.chi2_lists = None
        self.calc_chi2f = False

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
        for i in range(1, len(cumul)):
            diff.append(cumul[i] - cumul[i - 1])
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

        for i in range(par.get_min(), par.get_max() + 1):
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
        # report the naive chi^2 (no autocovariance) for shape only.
        # the model is scaled so that it matches the data on the final point
        # see description at top of this file
        scale = self.data[xdata[-1]] / self.model.populations[self.population_name].history[xdata[-1]]
        self.chi2d = 0.
        for x in xdata[:-1]:
            resid = self.data[x] - self.model.populations[self.population_name].history[x] * scale
            self.chi2d += resid ** 2 / (self.model.populations[self.population_name].history[x] * scale)

        return popt, pcov

    def calc_auto_covariance(self, n_rep=100):
        """ Calculate the autocovariance matrix. Do that by generating simulated datasets.

            The autocovariance matrix has rank one less than data since last point is used for normalization scale
            STILL TO DO: deal with daily data. Approach is tuned to study cumulative data for now.***
        """
        # Make a copy first

        sim_model = copy.deepcopy(self.model)
        sim_population = sim_model.populations[self.population_name]
        ref_population = self.model.populations[self.population_name]

        xdata = np.arange(self.data_range[0], self.data_range[1], 1)
        # last point not included: since its residual is zero by definition:
        n_p = len(xdata) - 1
        cov = [[0. for i in range(n_p)] for j in range(n_p)]
        for irep in range(n_rep):
            # produce simulated data -> sim_population histories
            sim_model.reset()
            sim_model.generate_data(self.data_range[1])

            # calculate residuals (data-expectation) for cov
            resid = []
            scale = sim_population.history[xdata[-1]] / ref_population.history[xdata[-1]]

            # last point not used for autocovariance xdata[:-1] removes that point
            for x in xdata[:-1]:
                delta = sim_population.history[x] - ref_population.history[x] * scale
                resid.append(delta)

            for i in range(n_p):
                for j in range(n_p):
                    cov[i][j] += resid[i] * resid[j] / n_rep

        self.auto_cov = cov

    def calc_sim_gof(self, n_rep=100):
        """ Calculate goodness of fit statistics for simulated data:
            chi2m: naive unfitted
            chi2f: naive fitted
            chi2s: unfitted (with autocovariance)
        """

        sim_model = copy.deepcopy(self.model)
        sim_population = sim_model.populations[self.population_name]
        ref_population = self.model.populations[self.population_name]

        # Make a copy for fitting simulated data
        sim_fit_model = copy.deepcopy(self.model)

        chi2_s_list = []
        chi2_m_list = []
        chi2_f_list = []
        chi2_n_list = []

        xdata = np.arange(self.data_range[0], self.data_range[1], 1)
        # last point not included: since its residual is zero by definition:
        n_p = len(xdata) - 1
        zeros = [0. for i in range(n_p)]
        lpdf_zero = stats.multivariate_normal.logpdf(zeros, cov=self.auto_cov)
        for i in range(n_rep):
            sim_model.reset()
            sim_model.generate_data(self.data_range[1])

            scale = sim_population.history[xdata[-1]] / ref_population.history[xdata[-1]]
            resid = []
            chi2_m = 0.
            for x in xdata[:-1]:
                delta = sim_population.history[x] - ref_population.history[x] * scale
                resid.append(delta)
                chi2_m += delta ** 2 / (ref_population.history[x] * scale)

            chi2_m_list.append(chi2_m)
            delta = sim_population.history[xdata[-1]] - ref_population.history[xdata[-1]]
            chi2_n = delta ** 2 / (ref_population.history[xdata[-1]])
            chi2_n_list.append(chi2_n)

            lpdf = stats.multivariate_normal.logpdf(resid, cov=self.auto_cov)
            # empirically it is found that the statistic follows a chi^2 distribution if the formula
            # is altered by changing the nominal 2 in this formula to 1:
            chi2_s_list.append(1. * (lpdf_zero - lpdf))

            if self.calc_chi2f:
                # do a fit to this dataset to see the naive gof for a fitted simulation sample
                # this needs to be turned on by the inquisitive user
                # ----------------------------------------------------------------------------
                sim_optimizer = Optimizer(sim_fit_model, self.full_population_name, sim_population.history,
                                          self.data_range)
                sim_optimizer.reset_variables()
                sim_popt, sim_pcov = sim_optimizer.fit()
                sim_ref_population = sim_fit_model.populations[self.population_name]

                # calculate residuals (data-expectation) for cov and chi^2
                resid_f = []
                scale_f = sim_population.history[xdata[-1]] / sim_ref_population.history[xdata[-1]]
                chi2_f = 0.
                # last point not used for autocovariance xdata[:-1] removes that point
                for x in xdata[:-1]:
                    delta = sim_population.history[x] - sim_ref_population.history[x] * scale_f
                    resid_f.append(delta)
                    chi2_f += delta ** 2 / (sim_ref_population.history[x] * scale_f)

                chi2_f_list.append(chi2_f)

        self.chi2m = np.mean(chi2_m_list)
        self.chi2m_sd = np.std(chi2_m_list)
        self.chi2n = np.mean(chi2_n_list)
        self.chi2n_sd = np.std(chi2_n_list)
        if self.calc_chi2f:
            self.chi2f = np.mean(chi2_f_list)
            self.chi2f_sd = np.std(chi2_f_list)
        self.chi2s = np.mean(chi2_s_list)
        self.chi2s_sd = np.std(chi2_s_list)
        # for further study of gof distributions:
        self.chi2_lists = {'m': chi2_m_list, 'f': chi2_f_list, 's': chi2_s_list, 'n': chi2_n_list}

    def mcmc(self, n_dof, chi2n_mean, n_MCMC):
        """ Make a MCMC chain (n_MCMC points) assuming n_dof defines gof statistic for data.
            The gof is calculated by shape (many dof) and normalization (one dof). The latter calculation
            requires a scaling factor: the mean chi2_n.
            Currently only setup to deal with float variables, but a simple extension would allow int/bool
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
        zeros = [0. for i in range(len(xdata)-1)]
        lpdf_zero = stats.multivariate_normal.logpdf(zeros, cov=self.auto_cov)

        # calculate log of posterior probability parameters
        def logP():
            # check if selected parameters within bounds (uniform priors)
            for var_name in var_names:
                if prior_function[var_name] == 'uniform':
                    par = var_parameters[var_name]
                    if np.abs(par.get_value() - mean[var_name]) > half_width[var_name]:
                        return -np.inf
                else:
                    par_val = par.get_value()
                    if par_val < par.get_min() or par_val > par.get_max():
                        return -np.inf

            sim_model.reset()
            sim_model.evolve_expectations(self.data_range[1])

            resid = []
            # do not use the final point
            for x in xdata[:-1]:
                resid.append(self.data[x] - sim_population.history[x])

            lpdf = stats.multivariate_normal.logpdf(resid, cov=self.auto_cov)
            # empirical finding that the statistic follows chi^2 if the
            # nominal 2 in the following form is replaced by 1.
            chi2 = 1. * (lpdf_zero - lpdf)

            # add to the shape chi^2, the normalization chi^2 scaled by chi2n_mean
            # this corresponds to one extra degree of freedom
            delta = self.data[xdata[-1]] - sim_population.history[xdata[-1]]
            chi2_n = delta ** 2 / (sim_population.history[xdata[-1]])
            chi2 += chi2_n/chi2n_mean

            ln_P = stats.chi2.logpdf(chi2, n_dof)

            # Deal with any normal priors
            for var_name in var_names:
                if prior_function[var_name] == 'normal':
                    par = var_parameters[var_name]
                    dev = (par.get_value() - mean[var_name]) / sigma[var_name]
                    ln_P -= dev * dev / 2.

            return ln_P

        chain = []

        n_accept = 0
        lp_curr = logP()
        for i in range(n_MCMC):
            new_var_values = var_values.copy()
            for var_name in var_names:
                new_var_values[var_name] += hypercube[var_name] * (1. - 2. * stats.uniform.rvs())
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

        self.accept_fraction = 1. * n_accept / n_MCMC
        return chain
