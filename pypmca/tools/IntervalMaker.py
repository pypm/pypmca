"""
Interval_maker:
 - produce simulated samples, with a distribution of transmission rates and other
   smeering parameters to calculate quantiles

 - this analysis is required for inclusion into the forecast hub
 - adds necessary data into model file (user_dict['forecast_hub'])

"""

import datetime
import numpy as np
from scipy import stats
import copy


class IntervalMaker:
    """ Run model many times to find intervals
    """

    def __init__(self, hub, forecast_date):
        # hub = "USA" or "Germany" for now
        # forecast_date: intended date for forecast (typically upcoming Sunday): datetime.date object
        self.hub = hub
        self.forecast_date = forecast_date
        self.quantile_dict = self.def_quantile_dict()
        self.period_dict = self.def_period_dict()
        self.point_estimates = None
        self.quantiles = None
        self.inc_periods = None
        self.sim_alphas = None
        self.population_name_dict = {'case':'reported', 'death':'deaths', 'hospitalization':'hospitalized'}

    def def_quantile_dict(self):
        # defines the default quantiles to produce for each category
        dict = {}
        if self.hub == 'USA':
            dict['case'] = [0.025, 0.1, 0.25, 0.5, 0.75, 0.9, 0.975]
            for category in ['death', 'hospitalization']:
                dict[category] = [0.01, 0.025] + [0.05 + 0.05 * i for i in range(19)] + [0.975, 0.99]
        if self.hub == 'Germany':
            for category in ['case','death']:
                dict[category] = [0.01, 0.025] + [0.05 + 0.05 * i for i in range(19)] + [0.975, 0.99]
        if self.hub == 'USA-scenario':
            for category in ['case', 'death', 'hospitalization']:
                dict[category] = [0.01, 0.025] + [0.05 + 0.05 * i for i in range(19)] + [0.975, 0.99]
        return dict

    def def_period_dict(self):
        # defines the default periods for each category
        dict = {}
        if self.hub == 'USA':
            dict['case'] = 'weekly'
            dict['death'] = 'weekly'
            dict['hospitalization'] = 'daily'
        if self.hub == 'Germany':
            dict['case'] = 'weekly'
            dict['death'] = 'weekly'
            dict['hospitalization'] = 'weekly'
        if self.hub == 'USA-scenario':
            dict['case'] = 'weekly'
            dict['death'] = 'weekly'
            dict['hospitalization'] = 'weekly'
        return dict

    def append_user_dict(self, category, model):
        if 'forecast_hub' not in model.user_dict:
            model.user_dict['forecast_hub'] = {}
        if self.forecast_date not in model.user_dict['forecast_hub']:
            model.user_dict['forecast_hub'][self.forecast_date] = {}
        model.user_dict['forecast_hub'][self.forecast_date][category] = {}
        model.user_dict['forecast_hub'][self.forecast_date][category]['point_estimates'] = self.point_estimates[category]
        model.user_dict['forecast_hub'][self.forecast_date][category]['quantiles'] = self.quantiles[category]
        model.user_dict['forecast_hub'][self.forecast_date][category]['inc_periods'] = self.inc_periods[category]
        model.user_dict['forecast_hub'][self.forecast_date][category]['sim_alphas'] = self.sim_alphas

    def get_quantiles(self, categories, n_periods_dict, model, n_rep=10, scale_std_alpha=1., back_up=7, rescale=False,
                      t0 = datetime.date(2020, 3, 1), fall_back=False):
        # categories: an array containing a combination of 'case', 'death', 'hospitalization'
        # n_periods_for_categories: a dictionary containing number of weeks for case/death and number of days for hospitalizations
        # n_rep: number of repetitions to produce quantiles
        # scale_std_alpha: if necessary the alpha parameter variation can be adjusted by this scaling factor
        # fall_back: allow the norm_day (start of data production) begin before the last transition date
        
        # The following are dictionaries indexed by categories
        self.point_estimates = {}
        self.quantiles = {}
        self.inc_periods = {}
        
        self.sim_alphas = []

        if scale_std_alpha != 1.:
            print('Warning: scaling the uncertainty for alpha:', str(scale_std_alpha))

        # epi-week period calculations:
        day_of_week = self.forecast_date.weekday()
        days_after_t0 = (self.forecast_date - t0).days
        first_sunday = days_after_t0
        if day_of_week == 0:
            first_sunday -= 1
        elif day_of_week < 6:
            first_sunday += 6 - day_of_week

        # find maximum number of days of simulation required
        n_days = 0
        for category in categories:
            if self.period_dict[category] == 'weekly':
                n_days_c = first_sunday + n_periods_dict[category] * 7
                n_days = max(n_days, n_days_c)
            elif self.period_dict[category] == 'daily':
                n_days_c = days_after_t0 + n_periods_dict[category]
                n_days = max(n_days, n_days_c)
            self.point_estimates[category] = {}
            self.quantiles[category] = {}
            self.inc_periods[category] = [[] for i in range(n_periods_dict[category])]

        # norm_day is when the expectation propagation ends and the data generation begins
        # back this up by a week or more, so that reporting noise issues do not generate immediate
        # negative bias (for daily reports) and to produce variation in deaths
        norm_day = max(1,days_after_t0 - back_up)

        # quantile estimates
        # Do many repititions: each time use a new alpha_last and renormalize expectation
        # then produce data and use difference from point estimate and data at last day to force last
        # day of data to have same starting point for forecast

        # find last alpha transition, and min and max of all alphas (prior to forecast date)
        last_time = 0
        alpha_par = None
        alpha_min = model.parameters['alpha_0'].get_value()
        alpha_max = alpha_min
        for trans_name in model.transitions:
            if 'rate' in trans_name:
                trans = model.transitions[trans_name]
                if trans.enabled:
                    time = trans.transition_time.get_value()
                    if time < days_after_t0:
                        alpha_trans = trans.parameter_after.get_value()
                        alpha_min = min(alpha_min, alpha_trans)
                        alpha_max = max(alpha_max, alpha_trans)
                        time = trans.transition_time.get_value()
                        if time >= last_time:
                            alpha_par = trans.parameter_after
                            last_time = time
        alpha_min = min(0.11, alpha_min)
        alpha_name = alpha_par.name
        alpha = alpha_par.get_value()
        alpha_err = 0.001
        if alpha_par.std_estimator is not None:
            alpha_err = alpha_par.std_estimator
        else:
            print('Warning: no std_estimator found for', alpha_par.name)

        if not fall_back and last_time >= norm_day:
            print('Warning: backup parameter too large: last_time (',last_time,') >= norm_day (',norm_day,')')
            norm_day = last_time + 1
            print(' -> norm_day changed to',norm_day)

        # adjust alpha uncertainty
        alpha_err *= scale_std_alpha

        # find expectations at start of forecast - to correct for variance produced in last section of generated data
        model.reset()
        model.evolve_expectations(days_after_t0-1)
        case_history = model.populations['reported'].history
        expected_cases = case_history[-1] - case_history[0]

        # define time for which copies of simulation models are made
        sim_ref_time = last_time
        if fall_back:
            sim_ref_time = min(norm_day,last_time)

        # run model with expectations to step before last alpha transition (prior to forecast date)
        model.reset()
        model.evolve_expectations(sim_ref_time - 1)
        sim_model_ref = copy.deepcopy(model)

        for i_rep in range(n_rep):
            sim_alpha = stats.norm.rvs(loc=alpha, scale=alpha_err)
            sim_trial = sim_alpha

            # protect against alphas far outside of observed range
            #if category == 'case':
            sim_alpha = max(sim_alpha, alpha_min * 0.9)
            sim_alpha = min(sim_alpha, alpha_max * 1.1)
            if sim_trial != sim_alpha:
                print('simulated alpha truncated: from', sim_trial, 'to', sim_alpha)

            # start with expectations up until last transition day
            sim_model = copy.deepcopy(sim_model_ref)
            sim_model.parameters[alpha_name].set_value(sim_alpha)

            # include independent variance from other sources
            if 'interval_maker' in model.user_dict:
                if 'smearing parameters' in model.user_dict['interval_maker']:
                    for par_name in model.user_dict['interval_maker']['smearing parameters']:
                        mean = model.parameters[par_name].get_value()
                        sigma = model.parameters[par_name].std_estimator
                        new_value = stats.norm.rvs(loc=mean, scale=sigma)
                        # if this is a parameter that is limited to between 0 and 1, assume a beta distribution
                        p_min = model.parameters[par_name].get_min()
                        p_max = model.parameters[par_name].get_max()
                        if p_min == 0. and p_max == 1.:
                            term = (mean * (1.-mean))/sigma**2 - 1.
                            a = term * mean
                            b = term * (1.-mean)
                            new_value = stats.beta.rvs(a,b)
                        # if this is an integer parameter: convert
                        if model.parameters[par_name].parameter_type == int:
                            new_value = int(new_value)
                        sim_model.parameters[par_name].set_value(new_value)

            # evolve_expectations up until the renormalization day (just before forecast period)
            sim_model.evolve_expectations(norm_day - sim_ref_time + 1, from_step=sim_ref_time - 1)

            # no renormalization performed... just convert to integers...

            for key in sim_model.populations:
                pop = sim_model.populations[key]
                nu = pop.history[norm_day]
                pop.history[norm_day] = int(round(nu))
                pop.scale_future(1., expectations=False)

            # now generate data starting from norm_day
            sim_model.generate_data(n_days - 1 - norm_day, from_step=norm_day)

            # derive a scaling factor to correct for variance produced since norm_day
            scale_factor = 1.
            if rescale:
                sim_case_history = sim_model.populations['reported'].history
                sim_cases = sim_case_history[days_after_t0] - sim_case_history[0]
                scale_factor = expected_cases / sim_cases

            for category in categories:
                population_name = self.population_name_dict[category]

                sim_population_history = sim_model.populations[population_name].history
                if self.period_dict[category] == 'weekly':
                    for week in range(n_periods_dict[category]):
                        end_epiweek = first_sunday - 1 + 7 * (week + 1)
                        inc_week = (sim_population_history[end_epiweek] -
                                    sim_population_history[end_epiweek - 7])
                        self.inc_periods[category][week].append(inc_week*scale_factor)
                elif self.period_dict[category] == 'daily':
                    for day in range(n_periods_dict[category]):
                        end_day = days_after_t0 + day
                        inc_day = (sim_population_history[end_day] -
                                   sim_population_history[end_day - 1])
                        self.inc_periods[category][day].append(inc_day*scale_factor)

            self.sim_alphas.append(sim_alpha)

        # point estimates
        model.reset()
        model.evolve_expectations(n_days - 1)

        #        quant = None
        #        if category in self.quantile_dict:
        #            quant = self.quantile_dict[category]

        for category in categories:
            population_name = self.population_name_dict[category]
            quant = self.quantile_dict[category]

            if self.period_dict[category] == 'weekly':
                for week in range(n_periods_dict[category]):
                    population_history = model.populations[population_name].history
                    end_epiweek = first_sunday - 1 + 7 * (week + 1)
                    value = population_history[end_epiweek] - population_history[end_epiweek - 7]
                    self.point_estimates[category][str(week + 1)] = value
    
                    quant_dict = {}
                    for quantile in quant:
                        value = np.percentile(self.inc_periods[category][week], quantile * 100.)
                        if value < 0.:
                            value = 0.
                        quantile_text = '{0:0.3f}'.format(quantile)
                        quant_dict[quantile_text] = value
                    self.quantiles[category][str(week + 1)] = quant_dict
    
            elif self.period_dict[category] == 'daily':
                for day in range(n_periods_dict[category]):
                    population_history = model.populations[population_name].history
                    end_day = days_after_t0 + day
                    value = population_history[end_day] - population_history[end_day - 1]
                    self.point_estimates[category][str(day + 1)] = value
    
                    quant_dict = {}
                    for quantile in quant:
                        value = np.percentile(self.inc_periods[category][day], quantile * 100.)
                        if value < 0.:
                            value = 0.
                        quantile_text = '{0:0.3f}'.format(quantile)
                        quant_dict[quantile_text] = value
                    self.quantiles[category][str(day + 1)] = quant_dict

