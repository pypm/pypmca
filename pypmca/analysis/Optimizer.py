# -*- coding: utf-8 -*-
"""
Optimizer: Given a model and data, find point estimates for the variable parameters in the model,
along with the covariance (optional)

This is done using the scipy.optimize.curve_fit method

The full_population_name has the prefix daily or total. If daily, then the
differences need to be calculated

The description below is for a "global fit". The local fit does not use the scaling factor that sets
the last point in the cumulative to match the data.

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

2021-1-31: Add skip_data (a string field) is added to define dates that should not be used in optimization and chi^2 calculation
Format example '250,263:270' means to remove 250,263,264,...,270

@author: karlen
"""
import copy
from scipy.optimize import curve_fit
from scipy import stats
import numpy as np
from datetime import timedelta


class Optimizer:
    """ Optimizer: point estimator and covariance
    """

    def __init__(self, model, full_population_name, data, data_range, cumul_reset=False, skip_data=None):
        self.model = model
        self.full_population_name = full_population_name
        self.population_name = full_population_name[6:]
        self.population_type = full_population_name[:5]
        self.data = data
        self.data_range = data_range
        self.cumul_reset = cumul_reset  # if True, use the cumulative starting at data_range[0] (a "local" fit)
        self.start_step = 0 # if local fit, this specifies the step to start each iteration (automatic, speeding up fit)
        self.model_ref = None # if local fit, this is the reference model that is copied
        self.variable_names = []  # float type variable parameters
        self.variable_initial_values = None
        self.i_variable_names = []  # integer type variable parameters
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
        self.opt_lists = None
        self.calc_chi2f = False
        self.calc_chi2s = False
        self.fit_statistics = None
        self.skip_data = skip_data
        self.skip_dates = None
        if skip_data is not None and skip_data != '':
            self.skip_dates = []
            blocks = skip_data.split(',')
            for block in blocks:
                if ':' in block:
                    limits = block.split(':')
                    for i in range(int(limits[0]),int(limits[1])+1):
                        self.skip_dates.append(i)
                else:
                    self.skip_dates.append(int(block))

    def remove_data(self, x, y):
        if self.skip_dates is None:
            return x,y
        else:
            x_rem = []
            y_rem = []
            for i,xi in enumerate(x):
                if xi not in self.skip_dates:
                    x_rem.append(xi)
                    y_rem.append(y[i])
            return x_rem,y_rem

    def func_setup(self):
        # if a local fit, setup the start_step
        transition_variable_steps = {}
        if self.cumul_reset:
            self.start_step = self.data_range[0]
            # get list of transitioning parameters and their dates
            for trans_name in self.model.transitions:
                transition = self.model.transitions[trans_name]
                if transition.enabled:
                    transition_variable_name = None
                    # Injector
                    injection_parameter = getattr(transition,'injection',None)
                    if injection_parameter is not None:
                        transition_variable_name = injection_parameter.name
                    else:
                        # Modifier
                        after_parameter = getattr(transition,'parameter_after',None)
                        if after_parameter is not None:
                            transition_variable_name = after_parameter.name
                    if transition_variable_name is not None:
                        transition_variable_steps[transition_variable_name] = transition.trigger_step

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
                    # check to see if this is a transitioning parameter, and if so adjust self.start_step
                    if self.cumul_reset:
                        if par_name in transition_variable_steps:
                            step = transition_variable_steps[par_name]
                            if step < self.start_step:
                                self.start_step = step

        # if local fit, make copy and run to start_step
        if self.cumul_reset:
            self.model_ref = copy.deepcopy(self.model)
            self.model_ref.reset()
            self.model_ref.evolve_expectations(self.start_step)

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
        # first daily value is repeated since val(t0-1) is unknown
        diff.insert(0,diff[0])
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
        """ work out point estimate. Experience shows small bias when auto_cov is used - turn off by default
        """

        def func(x, *params):

            if self.cumul_reset:
                # use a copy of the reference, to reduce calculation time to get to start_step
                model_copy = copy.deepcopy(self.model_ref)
                i = -1
                for par_val in params:
                    i += 1
                    par_name = self.variable_names[i]
                    model_copy.parameters[par_name].set_value(par_val)
                model_copy.evolve_expectations(self.data_range[1]-self.start_step, from_step=self.start_step)
                pop = model_copy.populations[self.population_name]

            else:
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
                cumul_offset = 0
                if self.cumul_reset and pop.monotonic:
                    cumul_offset = pop.history[self.data_range[0]]
                for xi in x:
                    i = int(round(xi))
                    func_values.append(pop.history[i]-cumul_offset)
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
            if self.population_type == 'total' and self.cumul_reset:
                cumul_offset = 0
                if self.model.populations[self.population_name].monotonic:
                    cumul_offset = self.data[self.data_range[0]]
                mod_data = [self.data[i] - cumul_offset for i in range(self.data_range[0],self.data_range[1])]
                x_rem, y_rem = self.remove_data(xdata, mod_data)
                popt, pcov = curve_fit(func, x_rem, y_rem,
                                       p0=par_0, bounds=bounds)
            else:
                x_rem, y_rem = self.remove_data(xdata, self.data[self.data_range[0]:self.data_range[1]])
                popt, pcov = curve_fit(func, x_rem, y_rem,
                                       p0=par_0, bounds=bounds)
        else:
            # with new autocovariance, last point is used for normalization remove from fit
            if self.population_type == 'total' and self.cumul_reset:
                cumul_offset = 0
                if self.model.populations[self.population_name].monotonic:
                    cumul_offset = self.data[self.data_range[0]]
                mod_data = [self.data[i] - cumul_offset for i in range(self.data_range[0],self.data_range[1]-1)]
                x_rem, y_rem = self.remove_data(xdata[:-1], mod_data)
                popt, pcov = curve_fit(func, x_rem, y_rem,
                                       p0=par_0, bounds=bounds, sigma=self.auto_cov)
            else:
                x_rem, y_rem = self.remove_data(xdata[:-1], self.data[self.data_range[0]:self.data_range[1]-1])
                popt, pcov = curve_fit(func, x_rem, y_rem,
                                        p0=par_0, bounds=bounds, sigma=self.auto_cov)

        i = -1
        for par_val in popt:
            i += 1
            par_name = self.variable_names[i]
            self.model.parameters[par_name].set_value(par_val)

        self.model.reset()
        self.model.evolve_expectations(self.data_range[1])
        # report the naive chi^2 (no autocovariance) for shape only.
        # the model is scaled so that it matches the data on the final point (for cumulative fits)
        # see description at top of this file
        scale = 1.

        cumul_offset = 0.
        cumul_offset_model = 0.
        if self.population_type == 'total' and self.cumul_reset:
            cumul_offset = self.data[self.data_range[0]]
            cumul_offset_model = self.model.populations[self.population_name].history[self.data_range[0]]

        if self.population_type == 'total' and not self.cumul_reset:
                scale = self.data[xdata[-1]] / self.model.populations[self.population_name].history[xdata[-1]]

        self.chi2d = 0.
        x_rem, y_rem = self.remove_data(xdata[:-1], xdata[:-1])
        if self.population_type == 'total':
            if self.cumul_reset:
                for x in x_rem:
                    resid = (self.data[x] - cumul_offset) - \
                            (self.model.populations[self.population_name].history[x] - cumul_offset_model)
                    self.chi2d += resid ** 2
            else:
                for x in x_rem:
                    resid = (self.data[x] - cumul_offset) - self.model.populations[self.population_name].history[x] * scale
                    self.chi2d += resid ** 2 / (self.model.populations[self.population_name].history[x] * scale)
        else:
            for x in x_rem:
                diff = self.delta(self.model.populations[self.population_name].history)
                resid = self.data[x] - diff[x]
                self.chi2d += resid ** 2 / (1. + diff[x])

        # report the chi^2 for the fit to daily data and the auto correlation of residuals to the next day
        # if reporting_days is not all days, then do not include non-reporting days or the following day
        # for either statistic
        data = self.data
        expectation = self.model.populations[self.population_name].history
        end_point = self.data_range[1]
        if self.population_type == 'total':
            data = self.delta(self.data)
            expectation = self.delta(np.array(expectation)*scale)

        self.fit_statistics = self.get_fit_statistics(self.model, self.population_name, data, expectation, self.data_range[0], end_point)
        # in addition include the fit statistic for the cumulative
        self.fit_statistics['chi2_c'] = self.chi2d

        return popt, pcov

    def get_fit_statistics(self, model, pop_name, data, expectation, start_point, end_point):
        # work out days of week to include:
        # if reporting_days is not all days, then do not include non-reporting days or the previous or following day
        all_days = True
        include_day = [True]*7
        try:
            report_noise, report_noise_par, report_days = model.populations[pop_name].get_report_noise()
        except:
            report_noise = False
            report_days = 127
        if report_noise and report_days != 127:
            for day_of_week in range(7):
                if report_days.get_value() & 2 ** day_of_week == 0:
                    include_day[day_of_week] = False
                    day_after = (day_of_week+1)%7
                    include_day[day_after] = False
                    day_before = (day_of_week+6)%7
                    include_day[day_before] = False
                    all_days = False

        chi2 = 0.
        sum_0 = 0.
        sum_1 = 0.
        count = 0
        x_rem, y_rem = self.remove_data(range(start_point, end_point), range(start_point, end_point))
        for t in x_rem:
            include_i = True
            if not all_days:
                # check that today should be included
                t0 = model.t0
                days_after = timedelta(days=int((t + 0.5) * model.get_time_step()))
                today = t0 + days_after
                day_of_week = today.weekday()
                include_i = include_day[day_of_week]
            if include_i:
                if len(data) > t + 1 and len(expectation) > t + 1:
                    denom = expectation[t]
                    if denom <= 0.:
                        denom = 0.01
                    chi2 += (data[t] - expectation[t])**2 / denom
                    sum_0 += (data[t] - expectation[t])**2
                    sum_1 += (data[t] - expectation[t]) * (data[t + 1] - expectation[t + 1])
                    count += 1
        cov = 0.
        acor = 0.
        if count > 0:
            cov = sum_0/count
            acor = sum_1/sum_0
        return {'ndof':count, 'chi2':chi2, 'cov':cov, 'acor':acor}

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
            chi2n: normalization

            This also calculates the properties of estimators: bias and variance

            The ensemble of simulations should only include situations which is similar to expectation.
            It is not useful to include a simulation which has the epidemic that dies out, for
            example. A simple condition is applied: the increment to the population for the last 4 days in the simulation
            should be within a factor of 2 of the expectation.

            If a local fit was done, then evolve expectations up until start of fit, then generate data. The conditions
            mentioned above are not required.
        """

        if self.cumul_reset:
            sim_model = None
        else:
            sim_model = copy.deepcopy(self.model)
        ref_population = self.model.populations[self.population_name]

        chi2_s_list = []
        chi2_m_list = []
        chi2_f_list = []
        chi2_n_list = []
        popt_list = []
        pcov_list = []
        fit_stat_list = []

        condition_days = 5
        condition_sum_expectation = 0
        if not self.cumul_reset:
            if ref_population.monotonic:
                condition_sum_expectation = ref_population.history[self.data_range[1]] - \
                                            ref_population.history[self.data_range[1]-condition_days]
            else:
                for day in range(condition_days):
                    condition_sum_expectation += ref_population.history[self.data_range[1]-condition_days]

        xdata = np.arange(self.data_range[0], self.data_range[1], 1)
        if self.cumul_reset:
            n_p = len(xdata)
        else:
            # last point not included: since its residual is zero by definition:
            n_p = len(xdata) - 1
        zeros = [0. for i in range(n_p)]
        if self.calc_chi2s:
            lpdf_zero = stats.multivariate_normal.logpdf(zeros, cov=self.auto_cov)
        for i in range(n_rep):
            if self.cumul_reset:
                sim_model = copy.deepcopy(self.model_ref)
                sim_model.generate_data(self.data_range[1]-self.start_step, from_step=self.start_step, data_start=self.data_range[0])
                sim_population = sim_model.populations[self.population_name]
            else:
                # require the simulation to produce data similar in scale to that observed
                condition_met = False
                while not condition_met:
                    sim_model.reset()
                    sim_model.generate_data(self.data_range[1])
                    sim_population = sim_model.populations[self.population_name]

                    condition_sum_simulation = 0
                    if sim_population.monotonic:
                        condition_sum_simulation = sim_population.history[self.data_range[1]] - \
                                                   sim_population.history[self.data_range[1] - condition_days]
                    else:
                        for day in range(condition_days):
                            condition_sum_simulation += sim_population.history[self.data_range[1] - condition_days]

                    condition_met = condition_sum_simulation < 2. * condition_sum_expectation and \
                                    condition_sum_simulation > 0.5 * condition_sum_expectation

            scale = sim_population.history[xdata[-1]] / ref_population.history[xdata[-1]]
            resid = []
            chi2_m = 0.
            x_rem, y_rem = self.remove_data(xdata[:-1], xdata[:-1])
            for x in x_rem:
                delta = sim_population.history[x] - ref_population.history[x] * scale
                resid.append(delta)
                chi2_m += delta ** 2 / (ref_population.history[x] * scale)

            chi2_m_list.append(chi2_m)
            delta = sim_population.history[xdata[-1]] - ref_population.history[xdata[-1]]
            chi2_n = delta ** 2 / (ref_population.history[xdata[-1]])
            chi2_n_list.append(chi2_n)

            if self.calc_chi2s:
                lpdf = stats.multivariate_normal.logpdf(resid, cov=self.auto_cov)
                # empirically it is found that the statistic follows a chi^2 distribution if the formula
                # is altered by changing the nominal 2 in this formula to 1:
                chi2_s_list.append(1. * (lpdf_zero - lpdf))

            if self.calc_chi2f:
                # do a fit to this dataset to see the naive gof for a fitted simulation sample
                # this needs to be turned on when estimating the properties of the estimators
                # ----------------------------------------------------------------------------
                # Make a copy for fitting simulated data - always start with same initial values...
                sim_fit_model = copy.deepcopy(self.model)
                sim_optimizer = Optimizer(sim_fit_model, self.full_population_name, sim_population.history,
                                          self.data_range, cumul_reset=self.cumul_reset, skip_data=self.skip_data)
                # following is needed only if re-using optimizer for a new fit...
                #sim_optimizer.reset_variables()
                sim_popt, sim_pcov = sim_optimizer.fit()
                popt_list.append(sim_popt)
                pcov_list.append(sim_pcov)

                sim_ref_population = sim_fit_model.populations[self.population_name]
                # calculate residuals (data-expectation) for cov and chi^2
                resid_f = []
                scale_f = sim_population.history[xdata[-1]] / sim_ref_population.history[xdata[-1]]
                chi2_f = 0.
                if self.cumul_reset:
                    for x in x_rem:
                        delta = sim_population.history[x] - sim_ref_population.history[x]
                        resid_f.append(delta)
                        chi2_f += delta ** 2
                else:
                    # last point not used for autocovariance xdata[:-1] removes that point
                    x_rem, y_rem = self.remove_data(xdata[:-1], xdata[:-1])
                    for x in x_rem:
                        delta = sim_population.history[x] - sim_ref_population.history[x] * scale_f
                        resid_f.append(delta)
                        chi2_f += delta ** 2 / (sim_ref_population.history[x] * scale_f)

                chi2_f_list.append(chi2_f)
                fit_stat_list.append(sim_optimizer.fit_statistics)

        self.chi2m = np.mean(chi2_m_list)
        self.chi2m_sd = np.std(chi2_m_list)
        self.chi2n = np.mean(chi2_n_list)
        self.chi2n_sd = np.std(chi2_n_list)
        if self.calc_chi2f:
            self.chi2f = np.mean(chi2_f_list)
            self.chi2f_sd = np.std(chi2_f_list)
        if self.calc_chi2s:
            self.chi2s = np.mean(chi2_s_list)
            self.chi2s_sd = np.std(chi2_s_list)
        # for further study of gof distributions:
        self.chi2_lists = {'m': chi2_m_list, 'f': chi2_f_list, 's': chi2_s_list, 'n': chi2_n_list}
        self.opt_lists = {'opt':popt_list, 'cov':pcov_list}
        self.fit_stat_list = fit_stat_list

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
