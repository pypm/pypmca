"""
Interval_maker:
 - produce simulated samples, with a distribution of transmission rates, to calculate quantiles

 - this analysis is required for inclusion into the forecast hub
 - adds necessary data into model file (user_dict['forecast_hub'])

"""

import datetime
import numpy as np
from scipy import stats
import copy
from pypmca import Model


class IntervalMaker:
    """ Run model many times to find intervals
    """

    def __init__(self, hub, forecast_date):
        # hub = "USA" or "Germany" for now
        # forecast_date: intended date for forecast (typically upcoming Sunday): datetime.date object
        self.hub = hub
        self.forecast_date = forecast_date
        self.quant = []
        self.point_estimates = {}
        self.quantiles = {}
        self.inc_periods = None
        self.sim_alphas = []

    def set_quantiles(self, category):
        self.quant = None
        if self.hub == 'USA':
            if category in ['case']:
                self.quant = [0.025, 0.1, 0.25, 0.5, 0.75, 0.9, 0.975]
            elif category in ['death', 'hospitalization']:
                self.quant = [0.01, 0.025] + [0.05 + 0.05 * i for i in range(19)] + [0.975, 0.99]
        if self.hub == 'Germany':
            self.quant = [0.01, 0.025] + [0.05 + 0.05 * i for i in range(19)] + [0.975, 0.99]

    def append_user_dict(self, category, model):
        if 'forecast_hub' not in model.user_dict:
            model.user_dict['forecast_hub'] = {}
        if self.forecast_date not in model.user_dict['forecast_hub']:
            model.user_dict['forecast_hub'][self.forecast_date] = {}
        model.user_dict['forecast_hub'][self.forecast_date][category] = {}
        model.user_dict['forecast_hub'][self.forecast_date][category]['point_estimates'] = self.point_estimates
        model.user_dict['forecast_hub'][self.forecast_date][category]['quantiles'] = self.quantiles
        model.user_dict['forecast_hub'][self.forecast_date][category]['inc_periods'] = self.inc_periods
        model.user_dict['forecast_hub'][self.forecast_date][category]['sim_alphas'] = self.sim_alphas

    def get_quantiles(self, category, n_periods, model, n_rep=10, scale_std_alpha=1.):
        # category: 'case', 'death', 'hospitalization'
        # n_periods: number of weeks for case/death and number of days for hospitalizations
        # n_rep: number of repetitions to produce quantiles
        self.point_estimates = {}
        self.quantiles = {}
        self.inc_periods = None
        self.sim_alphas = []

        if scale_std_alpha != 1.:
            print('Warning: scaling the uncertainty for alpha:', str(scale_std_alpha))

        t0 = datetime.date(2020, 3, 1)

        # epi-week period calculations:
        day_of_week = self.forecast_date.weekday()
        days_after_t0 = (self.forecast_date - t0).days
        first_sunday = days_after_t0
        if day_of_week == 0:
            first_sunday -= 1
        elif day_of_week < 6:
            first_sunday += 6 - day_of_week

        n_days = 0
        if category in ['case', 'death']:
            n_days = first_sunday + n_periods * 7
        elif category in ['hospitalization']:
            n_days = days_after_t0 + n_periods

        self.set_quantiles(category)

        population_name = None
        if category == 'case':
            population_name = 'reported'
        elif category == 'death':
            population_name = 'deaths'
        elif category == 'hospitalization':
            population_name = 'hospitalized'

        # norm_day is when the expectation propagation ends and the data generation begins
        # back this up by a few days, so that reporting noise issues do not generate immediate
        # negative bias (for daily reports)
        norm_day = days_after_t0 - 4

        # quantile estimates
        # Do many repititions: each time use a new alpha_last and renormalize expectation
        # then produce data and use difference from point estimate and data at last day to force last
        # day of data to have same starting point for forecast

        # find last alpha transition, and min and max of all alphas
        last_time = 0
        alpha_par = None
        alpha_min = model.parameters['alpha_0'].get_value()
        alpha_max = alpha_min
        for trans_name in model.transitions:
            if 'rate' in trans_name:
                trans = model.transitions[trans_name]
                if trans.enabled:
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

        # adjust alpha uncertainty
        alpha_err *= scale_std_alpha

        # compare two methods
        # #2 simulate expectations until just before last alpha transition - save that model
        # Then, each repetition only simulate data from that point onwards (much less simulation needed)
        # This is like a boot that brings state to the last alpha transition

        # run model with expectations to step before last alpha transition
        model.reset()
        model.evolve_expectations(last_time - 1)
        sim_model_ref = copy.deepcopy(model)

        # run model with expectations to normalization day
        model.evolve_expectations(norm_day - last_time + 1, from_step=last_time - 1)

        self.inc_periods = [[] for i in range(n_periods)]

        for i_rep in range(n_rep):
            sim_alpha = stats.norm.rvs(loc=alpha, scale=alpha_err)
            sim_trial = sim_alpha

            # protect against alphas outside of observed range (for case only)
            if category == 'case':
                sim_alpha = max(sim_alpha, alpha_min * 0.9)
                sim_alpha = min(sim_alpha, alpha_max * 1.1)
            if sim_trial != sim_alpha:
                print('simulated alpha truncated: from', sim_trial, 'to', sim_alpha)

            # start with expectations up until last transition day
            sim_model = copy.deepcopy(sim_model_ref)
            sim_model.parameters[alpha_name].set_value(sim_alpha)
            # evolve_expectations up until the renormalization day (just before forecast period)
            sim_model.evolve_expectations(norm_day - last_time + 1, from_step=last_time - 1)

            # no renormalization performed... just convert to integers...

            for key in sim_model.populations:
                pop = sim_model.populations[key]
                nu = pop.history[norm_day]
                pop.history[norm_day] = int(round(nu))
                pop.scale_future(1., expectations=False)

            # now generate data starting from norm_day
            sim_model.generate_data(n_days - 1 - norm_day, from_step=norm_day)

            sim_population_history = sim_model.populations[population_name].history
            if category in ['case', 'death']:
                for week in range(n_periods):
                    end_epiweek = first_sunday - 1 + 7 * (week + 1)
                    inc_week = (sim_population_history[end_epiweek] -
                                sim_population_history[end_epiweek - 7])
                    self.inc_periods[week].append(inc_week)
            elif category in ['hospitalization']:
                for day in range(n_periods):
                    end_day = days_after_t0 + day
                    inc_day = (sim_population_history[end_day] -
                               sim_population_history[end_day - 1])
                    self.inc_periods[day].append(inc_day)
            self.sim_alphas.append(sim_alpha)

        # point estimates
        model.reset()
        model.evolve_expectations(n_days - 1)

        if category in ['case', 'death']:
            for week in range(n_periods):
                population_history = model.populations[population_name].history
                end_epiweek = first_sunday - 1 + 7 * (week + 1)
                value = population_history[end_epiweek] - population_history[end_epiweek - 7]
                self.point_estimates[str(week + 1)] = value

                quant_dict = {}
                for quantile in self.quant:
                    value = np.percentile(self.inc_periods[week], quantile * 100.)
                    if value < 0.:
                        value = 0.
                    quantile_text = '{0:0.3f}'.format(quantile)
                    quant_dict[quantile_text] = value
                self.quantiles[str(week + 1)] = quant_dict

        elif category in ['hospitalization']:
            for day in range(n_periods):
                population_history = model.populations[population_name].history
                end_day = days_after_t0 + day
                value = population_history[end_day] - population_history[end_day - 1]
                self.point_estimates[str(day + 1)] = value

                quant_dict = {}
                for quantile in self.quant:
                    value = np.percentile(self.inc_periods[day], quantile * 100.)
                    if value < 0.:
                        value = 0.
                    quantile_text = '{0:0.3f}'.format(quantile)
                    quant_dict[quantile_text] = value
                self.quantiles[str(day + 1)] = quant_dict

