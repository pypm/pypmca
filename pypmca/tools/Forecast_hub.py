# -*- coding: utf-8 -*-
"""
Produce .csv file suitable for submission to covid19-forecast-hub (https://github.com/reichlab/covid19-forecast-hub)



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


class Forecast_hub:
    """ Forecast_hub: produce .csv file for covid19-forecast-hub
    """

    def __init__(self, folder, model_names):

        self.model_dir = Path(folder).resolve()
        self.model_names = model_names
        self.pd_dict = self.get_data()
        self.buff = [['forecast_date', 'target', 'target_end_date', 'location', 'type', 'quantile', 'value']]

        # number of daysfor case_normalization
        self.case_norm = 14
        # number of days for death_normalization
        self.death_norm = 14

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

    def test_Model_copy_values_from(self):
        path_model_2 = self.model_dir / 'ref_model_2.pypm'
        ref_model_2 = Model.open_file(path_model_2)

    def get_csv(self,forecast_date,us_deaths,cor_scale=1.):
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
        periods = ['wk', 'wk', 'day']
        inc_types_list = [['inc'], ['inc', 'cum'], ['inc']]

        # collect information for full US
        us_inc_point_est_dict = {}
        us_inc_periods_dict = {}
        us_cum_point_est_dict = {}
        us_cum_periods_dict = {}
        for dict_name in dict_names:
            us_inc_point_est_dict[dict_name] = {}
            us_inc_periods_dict[dict_name] = {}
        us_cum_point_est_dict['death'] = {}
        us_cum_periods_dict['death'] = {}
        state_deaths = 0

        for state in self.fips_code:
        #for state in ['TX','SC','FL']:
            abbrev = state.lower()
            location = self.fips_code[state]
            deaths = self.pd_dict['usa-jhu-pypm.csv'][abbrev.upper() + '-dt'].fillna(0).values

            success = False
            for model_name in self.model_names:
                try:
                    filename = abbrev + model_name + '.pypm'
                    path_model = self.model_dir / filename
                    case_model = Model.open_file(path_model)
                    success = True
                    break
                except:
                    pass
            for model_name in self.model_names:
                for model_suffix in ['_d','_h','']:
                    try:
                        filename = abbrev + model_name + model_suffix + '.pypm'
                        path_model = self.model_dir / filename
                        death_model = Model.open_file(path_model)
                        break
                    except:
                        pass
            for model_name in self.model_names:
                for model_suffix in ['_h','']:
                    try:
                        filename = abbrev + model_name + model_suffix + '.pypm'
                        path_model = self.model_dir / filename
                        hosp_model = Model.open_file(path_model)
                        break
                    except:
                        pass
            if not success:
                raise RuntimeError('No model for: ', abbrev)

            models = [case_model, death_model, hosp_model]
            print(models[0].name, models[1].name, models[2].name)

            for i in range(3):
                model = models[i]
                inc_types = inc_types_list[i]

                hub_dict = model.user_dict['forecast_hub'][forecast_date][dict_names[i]]
                point_est_dict = hub_dict['point_estimates']
                quantile_dict = hub_dict['quantiles']
                inc_periods = hub_dict['inc_periods']
                cum_periods = [deaths[first_sunday-1]] * len(inc_periods[0])
                for inc_type in inc_types:
                    # for cum record (so far, only deaths needs this)
                    sum = deaths[first_sunday-1]
                    if inc_type == 'cum':
                        state_deaths += sum
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
                        if inc_type == 'cum':
                            sum += point_est_dict[i_period]
                            value = sum

                            if i_period in us_cum_point_est_dict['death']:
                                us_cum_point_est_dict['death'][i_period] += sum
                            else:
                                us_cum_point_est_dict['death'][i_period] = sum

                        self.add_record(forecast_date_text,target,target_end_date,location,'point','NA',value)

                        # quantiles
                        if inc_type == 'inc':
                            if i_period in us_inc_point_est_dict[dict_names[i]]:
                                us_inc_point_est_dict[dict_names[i]][i_period] += value
                            else:
                                us_inc_point_est_dict[dict_names[i]][i_period] = value

                            for quant in quantile_dict[i_period]:
                                value = quantile_dict[i_period][quant]
                                self.add_record(forecast_date_text, target, target_end_date, location, 'quantile', quant, value)

                            if i_period in us_inc_periods_dict[dict_names[i]]:
                                ip_len = len(inc_periods[index])
                                for i_rep in range(len(us_inc_periods_dict[dict_names[i]][i_period])):
                                    # in case there are fewer repetitions, loop over others:
                                    ip_index = i_rep % ip_len
                                    us_inc_periods_dict[dict_names[i]][i_period][i_rep] += inc_periods[index][ip_index]
                            else:
                                us_inc_periods_dict[dict_names[i]][i_period] = copy.copy(inc_periods[index])
                        elif inc_type == 'cum':
                            for i_rep in range(len(inc_periods[index])):
                                cum_periods[i_rep] += inc_periods[index][i_rep]
                            for quant in quantile_dict[i_period]:
                                value = np.percentile(cum_periods, float(quant)*100.)
                                self.add_record(forecast_date_text, target, target_end_date, location, 'quantile', quant, value)

                            if i_period in us_cum_periods_dict[dict_names[i]]:
                                ip_len = len(cum_periods)
                                for i_rep in range(len(us_cum_periods_dict[dict_names[i]][i_period])):
                                    # in case there are fewer repetitions, loop over others:
                                    ip_index = i_rep % ip_len
                                    us_cum_periods_dict[dict_names[i]][i_period][i_rep] += cum_periods[ip_index]
                            else:
                                us_cum_periods_dict[dict_names[i]][i_period] = [cum_periods[i_rep] for i_rep in range(len(cum_periods))]

        # return 'US' summary:
        # there are additional regions used to define total US deaths Need to correct by adding the additional deaths here
        additional_deaths = us_deaths - state_deaths
        print('US total: additional deaths included:', additional_deaths)

        location = 'US'
        for i in range(3):
            inc_types = inc_types_list[i]
            if dict_names[i] in ['case']:
                quants = [0.025, 0.1, 0.25, 0.5, 0.75, 0.9, 0.975]
            elif dict_names[i] in ['death', 'hospitalization']:
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
                    if inc_type == 'cum':
                        value = us_cum_point_est_dict[dict_names[i]][i_period] + additional_deaths
                    self.add_record(forecast_date_text, target, target_end_date, location, 'point', 'NA', value)

                    # quantiles
                    if inc_type == 'inc':
                        median = np.percentile(us_inc_periods_dict[dict_names[i]][i_period], 50.)
                        for quant in quants:
                            value = np.percentile(us_inc_periods_dict[dict_names[i]][i_period], float(quant)*100.)
                            # increase the width of the interval for the US as a whole, to account for correlations between states
                            scaled_value = max(0., cor_scale*(value-median) + median)
                            quant_text = '{0:0.3f}'.format(quant)
                            self.add_record(forecast_date_text, target, target_end_date, location, 'quantile', quant_text,
                                            scaled_value)

                    elif inc_type == 'cum':
                        median = np.percentile(us_cum_periods_dict[dict_names[i]][i_period], 50.) + additional_deaths
                        for quant in quants:
                            value = np.percentile(us_cum_periods_dict[dict_names[i]][i_period], float(quant) * 100.) + additional_deaths
                            scaled_value = max(0., cor_scale * (value - median) + median)
                            quant_text = '{0:0.3f}'.format(quant)
                            self.add_record(forecast_date_text, target, target_end_date, location, 'quantile', quant_text,
                                            scaled_value)

        return self.buff



    def add_record(self, forecast_date, target, target_end_date, location, record_type, quantile_text, value):
        record = [forecast_date, target, target_end_date, location, record_type, quantile_text]
        record.append('{0:0.1f}'.format(value))
        self.buff.append(record)

my_forecast = Forecast_hub('/Users/karlen/pypm-temp/usa', ['_2_9_0801'])
# Indicate the total US deaths (up to and including Saturday) here:
us_deaths = 613157

# changed to 3 for Feb 28 - all have variant - large correlated uncertainty! - Apr 25: return to 1.5 (single strain)
my_csv = my_forecast.get_csv(datetime.date(2021, 8, 1), us_deaths, cor_scale=1.5)
pass
with open('/Users/karlen/pypm-temp/test-forecast.csv','w') as out:
    for line in my_csv:
        record = ','.join(line)
        out.write(record + '\n')
