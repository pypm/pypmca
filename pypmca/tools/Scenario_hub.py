# -*- coding: utf-8 -*-
"""
Produce .csv file suitable for submission to covid19 scenario-hub (https://github.com/midas-network/covid19-scenario-modeling-hub)

@author: karlen
"""
from pathlib import Path
import numpy as np
import pandas as pd
from scipy import stats
import requests
import datetime
from datetime import timedelta
import copy

from pypmca import Model, Population, Delay, Parameter, Multiplier, Propagator, \
    Splitter, Adder, Subtractor, Chain, Modifier, Injector, Ensemble


class Scenario_hub:
    """ Scenario_hub: produce .csv file for covid19 scenario-hub
    """

    def __init__(self, folder, model_names, scenario_names, scenario_abbrevs, scenario_ids):

        self.model_dir = Path(folder).resolve()
        self.model_names = model_names
        self.scenario_names = scenario_names
        self.scenario_abbrevs = scenario_abbrevs
        self.scenario_ids = scenario_ids
        self.pd_dict = self.get_data()
        self.buff = [['model_projection_date', 'target', 'target_end_date', 'location', 'type', 'quantile', 'scenario_name', 'scenario_id', 'value']]

        self.fips_code = {
            'AL': '01',
            'AK': '02',
            'AZ': '04',
            'AR': '05',
            'CA': '06',
            'CO': '08',
            'CT': '09',
            'DE': '10',
            'DC': '11',
            'FL': '12',
            'GA': '13',
            'HI': '15',
            'ID': '16',
            'IL': '17',
            'IN': '18',
            'IA': '19',
            'KS': '20',
            'KY': '21',
            'LA': '22',
            'ME': '23',
            'MD': '24',
            'MA': '25',
            'MI': '26',
            'MN': '27',
            'MS': '28',
            'MO': '29',
            'MT': '30',
            'NE': '31',
            'NV': '32',
            'NH': '33',
            'NJ': '34',
            'NM': '35',
            'NY': '36',
            'NC': '37',
            'ND': '38',
            'OH': '39',
            'OK': '40',
            'OR': '41',
            'PA': '42',
            'PR': '72',
            'RI': '44',
            'SC': '45',
            'SD': '46',
            'TN': '47',
            'TX': '48',
            'UT': '49',
            'VT': '50',
            'VA': '51',
            'WA': '53',
            'WV': '54',
            'WI': '55',
            'WY': '56'
        }

    def get_data(self):
        data_folder = 'data/covid19/USA'

        success = True
        pd_dict = {}
        try:
            data_desc_resp = requests.get('http://data.ipypm.ca/get_data_desc/' + data_folder)
        except requests.exceptions.RequestException as error:
            print('Error retrieving data description over network:')
            success = False
        if success:
            data_description = data_desc_resp.json()
            data_description['folder'] = data_folder
            print(data_description['description'])

            # load the data into a panda dictionary
            for filename in data_description['files']:
                path = data_folder + '/' + filename
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

    def get_csv(self,forecast_date,us_deaths,us_cases,cor_scale=1.):
        t0 = datetime.date(2020, 3, 1)
        forecast_date_text = forecast_date.isoformat()
        day_of_week = forecast_date.weekday()
        days_after_t0 = (forecast_date - t0).days
        first_sunday = days_after_t0
        if day_of_week == 0:
            first_sunday -= 1
        elif day_of_week < 6:
            first_sunday += 6 - day_of_week

        names = ['case', 'death', 'hosp']
        dict_names = ['case', 'death', 'hospitalization']
        periods = ['wk', 'wk', 'wk']
        data_types = ['cum','cum','inc']
        inc_types_list = [['inc','cum inc'], ['inc', 'cum'], ['inc','cum inc']]

        for i_s, scenario_name in enumerate(self.scenario_names):
            scenario_abbrev = self.scenario_abbrevs[i_s]
            scenario_id = self.scenario_ids[i_s]

            # collect information for full US
            us_inc_point_est_dict = {}
            us_inc_periods_dict = {}
            us_cum_point_est_dict = {}
            us_cum_periods_dict = {}
            us_total = {}
            state_sums = {}
            for i,dict_name in enumerate(dict_names):
                if 'inc' in inc_types_list[i]:
                    us_inc_point_est_dict[dict_name] = {}
                    us_inc_periods_dict[dict_name] = {}
                if 'cum' in inc_types_list[i] or 'cum inc' in inc_types_list[i]:
                    us_cum_point_est_dict[dict_name] = {}
                    us_cum_periods_dict[dict_name] = {}
                if data_types[i] == 'cum':
                    state_sums[dict_name] = 0

            us_total[dict_names[0]] = us_cases
            us_total[dict_names[1]] = us_deaths

            for state in self.fips_code:
            #for state in ['TX','SC','FL']:
                abbrev = state.lower()
                location = self.fips_code[state]

                state_data = {}
                state_data[dict_names[0]] = self.pd_dict['usa-jhu-pypm.csv'][abbrev.upper() + '-pt'].fillna(0).values
                state_data[dict_names[1]] = self.pd_dict['usa-jhu-pypm.csv'][abbrev.upper() + '-dt'].fillna(0).values

                success = False
                for model_name in self.model_names:
                    try:
                        filename = abbrev + model_name + scenario_abbrev + '.pypm'
                        path_model = self.model_dir / scenario_name / filename
                        case_model = Model.open_file(path_model)
                        success = True
                        break
                    except:
                        pass
                for model_name in self.model_names:
                    for model_suffix in ['_d','_h','']:
                        try:
                            filename = abbrev + model_name + model_suffix + scenario_abbrev + '.pypm'
                            path_model = self.model_dir / scenario_name / filename
                            death_model = Model.open_file(path_model)
                            break
                        except:
                            pass
                for model_name in self.model_names:
                    for model_suffix in ['_h','']:
                        try:
                            filename = abbrev + model_name + model_suffix + scenario_abbrev + '.pypm'
                            path_model = self.model_dir / scenario_name / filename
                            hosp_model = Model.open_file(path_model)
                            break
                        except:
                            pass
                if not success:
                    raise RuntimeError('No model for: ', abbrev)

                # These correspond to the names, dict_names above
                models = [case_model, death_model, hosp_model]
                print(models[0].name, models[1].name, models[2].name)

                for i in range(len(models)):
                    model = models[i]
                    inc_types = inc_types_list[i]

                    hub_dict = model.user_dict['forecast_hub'][forecast_date][dict_names[i]]
                    point_est_dict = hub_dict['point_estimates']
                    quantile_dict = hub_dict['quantiles']
                    inc_periods = hub_dict['inc_periods']
                    cum_periods = []
                    if data_types[i] == 'cum':
                        cum_periods = [state_data[dict_names[i]][first_sunday-1]] * len(inc_periods[0])
                    else:
                        cum_periods = [0] * len(inc_periods[0])
                    for inc_type in inc_types:
                        sum = 0
                        if data_types[i] == 'cum':
                            sum = state_data[dict_names[i]][first_sunday-1]
                        if inc_type == 'cum' or inc_type == 'cum inc':
                            if data_types[i] == 'cum':
                                state_sums[dict_names[i]] += sum
                        for i_period in point_est_dict:
                            index = int(i_period) - 1
                            target = i_period+' '+periods[i]+' ahead '+inc_type+' '+names[i]
                            days = int(i_period)
                            if periods[i] == 'wk':
                                days *= 7
                            else:
                                days += 1
                            target_end_date = (t0 + timedelta(days=(first_sunday + days - 1))).isoformat()

                            # point estimates
                            value = point_est_dict[i_period]
                            if inc_type == 'cum' or inc_type == 'cum inc':
                                sum += point_est_dict[i_period]
                                value = sum

                                if i_period in us_cum_point_est_dict[dict_names[i]]:
                                    us_cum_point_est_dict[dict_names[i]][i_period] += sum
                                else:
                                    us_cum_point_est_dict[dict_names[i]][i_period] = sum

                            self.add_record(forecast_date_text,target,target_end_date,location,'point','NA',value,scenario_name, scenario_id)

                            # quantiles
                            if inc_type == 'inc':
                                if i_period in us_inc_point_est_dict[dict_names[i]]:
                                    us_inc_point_est_dict[dict_names[i]][i_period] += value
                                else:
                                    us_inc_point_est_dict[dict_names[i]][i_period] = value

                                for quant in quantile_dict[i_period]:
                                    value = quantile_dict[i_period][quant]
                                    self.add_record(forecast_date_text, target, target_end_date, location, 'quantile', quant, value, scenario_name, scenario_id)

                                if i_period in us_inc_periods_dict[dict_names[i]]:
                                    ip_len = len(inc_periods[index])
                                    for i_rep in range(len(us_inc_periods_dict[dict_names[i]][i_period])):
                                        # in case there are fewer repetitions, loop over others:
                                        ip_index = i_rep % ip_len
                                        us_inc_periods_dict[dict_names[i]][i_period][i_rep] += inc_periods[index][ip_index]
                                else:
                                    us_inc_periods_dict[dict_names[i]][i_period] = copy.copy(inc_periods[index])
                            elif inc_type == 'cum' or inc_type == 'cum inc':
                                for i_rep in range(len(inc_periods[index])):
                                    cum_periods[i_rep] += inc_periods[index][i_rep]
                                for quant in quantile_dict[i_period]:
                                    value = np.percentile(cum_periods, float(quant)*100.)
                                    self.add_record(forecast_date_text, target, target_end_date, location, 'quantile', quant, value, scenario_name, scenario_id)

                                if i_period in us_cum_periods_dict[dict_names[i]]:
                                    ip_len = len(cum_periods)
                                    for i_rep in range(len(us_cum_periods_dict[dict_names[i]][i_period])):
                                        # in case there are fewer repetitions, loop over others:
                                        ip_index = i_rep % ip_len
                                        us_cum_periods_dict[dict_names[i]][i_period][i_rep] += cum_periods[ip_index]
                                else:
                                    us_cum_periods_dict[dict_names[i]][i_period] = [cum_periods[i_rep] for i_rep in range(len(cum_periods))]

            # return 'US' summary:
            # there are additional regions used to define total US deaths and deaths Need to correct by adding the additional amounts here
            additional = {}
            for i in range(len(dict_names)):
                if data_types[i] == 'cum':
                    additional[dict_names[i]] = us_total[dict_names[i]] - state_sums[dict_names[i]]
                    print('US total: additional',dict_names[i],'s included:', additional[dict_names[i]])

            location = 'US'
            for i in range(3):
                inc_types = inc_types_list[i]
                quants = []
                if dict_names[i] in ['case', 'death', 'hospitalization']:
                    quants = [0.01, 0.025] + [0.05 + 0.05 * i for i in range(19)] + [0.975, 0.99]
                for inc_type in inc_types:
                    for i_period in us_inc_point_est_dict[dict_names[i]]:
                        target = i_period + ' ' + periods[i] + ' ahead ' + inc_type + ' ' + names[i]
                        days = int(i_period)
                        if periods[i] == 'wk':
                            days *= 7
                        else:
                            days += 1
                        target_end_date = (t0 + timedelta(days=(first_sunday + days - 1))).isoformat()

                        # point estimates
                        value = us_inc_point_est_dict[dict_names[i]][i_period]
                        if inc_type == 'cum' or inc_type == 'cum inc':
                            value = us_cum_point_est_dict[dict_names[i]][i_period]
                            if data_types[i] == 'cum':
                                value += additional[dict_names[i]]
                        self.add_record(forecast_date_text, target, target_end_date, location, 'point', 'NA', value, scenario_name, scenario_id)

                        # quantiles
                        if inc_type == 'inc':
                            median = np.percentile(us_inc_periods_dict[dict_names[i]][i_period], 50.)
                            for quant in quants:
                                value = np.percentile(us_inc_periods_dict[dict_names[i]][i_period], float(quant)*100.)
                                # increase the width of the interval for the US as a whole, to account for correlations between states
                                scaled_value = cor_scale*(value-median) + median
                                quant_text = '{0:0.3f}'.format(quant)
                                self.add_record(forecast_date_text, target, target_end_date, location, 'quantile', quant_text,
                                                scaled_value, scenario_name, scenario_id)

                        elif inc_type == 'cum' or inc_type == 'cum inc':
                            median = np.percentile(us_cum_periods_dict[dict_names[i]][i_period], 50.)
                            if data_types[i] == 'cum':
                                median += additional[dict_names[i]]
                            for quant in quants:
                                value = np.percentile(us_cum_periods_dict[dict_names[i]][i_period], float(quant) * 100.)
                                if data_types[i] == 'cum':
                                    value += additional[dict_names[i]]
                                scaled_value = cor_scale * (value - median) + median
                                quant_text = '{0:0.3f}'.format(quant)
                                self.add_record(forecast_date_text, target, target_end_date, location, 'quantile', quant_text,
                                                scaled_value, scenario_name, scenario_id)

        return self.buff


    def add_record(self, forecast_date, target, target_end_date, location, record_type, quantile_text, value, scenario_name, scenario_id):
        record = [forecast_date, target, target_end_date, location, record_type, quantile_text, scenario_name, scenario_id]
        record.append('{0:0.1f}'.format(value))
        self.buff.append(record)

scenario_names = ['optimistic','moderate','fatigue','counterfactual']
scenario_abbrevs = ['_opt','_mod','_fat','_cfs']
scenario_ids = ['A-2020-12-22','B-2020-12-22','C-2020-12-22','D-2020-12-22']

my_scenario = Scenario_hub('/Users/karlen/pypm-temp/usa-scenario', ['_2_7_1229'], scenario_names, scenario_abbrevs, scenario_ids)
# Indicate the total US deaths (up to and including Saturday) here:
us_deaths = 331909
us_cases = 18982634


my_csv = my_scenario.get_csv(datetime.date(2020, 12, 27), us_deaths, us_cases, cor_scale=1.5)
pass
with open('/Users/karlen/pypm-temp/test-scenario.csv','w') as out:
    for line in my_csv:
        record = ','.join(line)
        out.write(record + '\n')

i=1