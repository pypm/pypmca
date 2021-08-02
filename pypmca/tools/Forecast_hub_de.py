# -*- coding: utf-8 -*-
"""
Produce .csv file suitable for submission to the Germany covid19-forecast-hub
(https://github.com/KITmetricslab/covid19-forecast-hub-de)

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

        self.regional_abbreviations = {
            'Baden-Wurttemberg': 'bw',
            'Bavaria': 'by',
            'Berlin': 'be',
            'Brandenburg': 'bb',
            'Bremen': 'hb',
            'Hamburg': 'hh',
            'Hesse': 'he',
            'Lower Saxony': 'ni',
            'Mecklenburg-Vorpommern': 'mv',
            'North Rhine-Westphalia': 'nw',
            'Rhineland-Palatinate': 'rp',
            'Saarland': 'sl',
            'Saxony': 'sn',
            'Saxony-Anhalt': 'st',
            'Schleswig-Holstein': 'sh',
            'Thuringia': 'th',
        }

        self.regional_locations = {
            'Baden-Wurttemberg': 'GM01',
            'Bavaria': 'GM02',
            'Bremen': 'GM03',
            'Hamburg': 'GM04',
            'Hesse': 'GM05',
            'Lower Saxony': 'GM06',
            'North Rhine-Westphalia': 'GM07',
            'Rhineland-Palatinate': 'GM08',
            'Saarland': 'GM09',
            'Schleswig-Holstein': 'GM10',
            'Brandenburg': 'GM11',
            'Mecklenburg-Vorpommern': 'GM12',
            'Saxony': 'GM13',
            'Saxony-Anhalt': 'GM14',
            'Thuringia': 'GM15',
            'Berlin': 'GM16',
            'Germany': 'GM'
        }

    def get_data(self):
        data_folder = 'data/covid19/Germany'

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

    def get_csv(self, forecast_date, category,cor_scale=1.):
        t0 = datetime.date(2020, 3, 1)
        forecast_date_text = forecast_date.isoformat()
        day_of_week = forecast_date.weekday()
        days_after_t0 = (forecast_date - t0).days
        first_sunday = days_after_t0
        if day_of_week == 0:
            first_sunday -= 1
        elif day_of_week < 6:
            first_sunday += 6 - day_of_week

        #names = ['case', 'death']
        #dict_names = ['case', 'death']
        #periods = ['wk', 'wk']
        #inc_types_list = [['inc', 'cum'], ['inc', 'cum']]

        period = 'wk'
        inc_types = ['inc', 'cum']

        # collect information for all Germany
        de_inc_point_est_dict = {}
        de_inc_periods_dict = {}
        de_cum_point_est_dict = {}
        de_cum_periods_dict = {}
        de_inc_point_est_dict[category] = {}
        de_inc_periods_dict[category] = {}
        de_cum_point_est_dict[category] = {}
        de_cum_periods_dict[category] = {}
        state_sums = 0

        for state in self.regional_abbreviations:
            abbrev = self.regional_abbreviations[state]
            location = self.regional_locations[state]
            state_data = None
            if category == 'case':
                state_data = self.pd_dict['germany-rki-pypm.csv'][abbrev + '-pt'].fillna(0).values
            elif category == 'death':
                state_data = self.pd_dict['germany-rki-pypm.csv'][abbrev + '-dt'].fillna(0).values

            success = False
            model = None

            if category == 'case':
                for model_name in self.model_names:
                    try:
                        filename = abbrev + model_name + '.pypm'
                        path_model = self.model_dir / filename
                        model = Model.open_file(path_model)
                        success = True
                        break
                    except:
                        pass
            elif category == 'death':
                for model_name in self.model_names:
                    for model_suffix in ['_d','_h','']:
                        try:
                            filename = abbrev + model_name + model_suffix + '.pypm'
                            path_model = self.model_dir / filename
                            model = Model.open_file(path_model)
                            success = True
                            break
                        except:
                            pass

            if not success:
                raise RuntimeError('No model for: ', abbrev)

            print(model.name)

            hub_dict = model.user_dict['forecast_hub'][forecast_date][category]
            point_est_dict = hub_dict['point_estimates']
            quantile_dict = hub_dict['quantiles']
            inc_periods = hub_dict['inc_periods']
            cum_periods = [state_data[first_sunday-1]] * len(inc_periods[0])
            for inc_type in inc_types:
                # for cum record
                sum = state_data[first_sunday-1]
                if inc_type == 'cum':
                    state_sums += sum
                for i_period in point_est_dict:
                    index = int(i_period) - 1
                    target = i_period+' '+period+' ahead '+inc_type+' '+category
                    days = int(i_period)
                    if period == 'wk':
                        days *= 7
                    else:
                        days += 1
                    target_end_date = (t0 + timedelta(days=(first_sunday + days - 1))).isoformat()

                    # point estimates
                    value = point_est_dict[i_period]
                    if inc_type == 'cum':
                        sum += point_est_dict[i_period]
                        value = sum

                        if i_period in de_cum_point_est_dict[category]:
                            de_cum_point_est_dict[category][i_period] += sum
                        else:
                            de_cum_point_est_dict[category][i_period] = sum

                    self.add_record(forecast_date_text,target,target_end_date,location,'point','NA',value)

                    # quantiles
                    if inc_type == 'inc':
                        if i_period in de_inc_point_est_dict[category]:
                            de_inc_point_est_dict[category][i_period] += value
                        else:
                            de_inc_point_est_dict[category][i_period] = value

                        for quant in quantile_dict[i_period]:
                            value = quantile_dict[i_period][quant]
                            self.add_record(forecast_date_text, target, target_end_date, location, 'quantile', quant, value)

                        if i_period in de_inc_periods_dict[category]:
                            ip_len = len(inc_periods[index])
                            for i_rep in range(len(de_inc_periods_dict[category][i_period])):
                                # in case there are fewer repetitions, loop over others:
                                ip_index = i_rep % ip_len
                                de_inc_periods_dict[category][i_period][i_rep] += inc_periods[index][ip_index]
                        else:
                            de_inc_periods_dict[category][i_period] = copy.copy(inc_periods[index])
                    elif inc_type == 'cum':
                        for i_rep in range(len(inc_periods[index])):
                            cum_periods[i_rep] += inc_periods[index][i_rep]
                        for quant in quantile_dict[i_period]:
                            value = np.percentile(cum_periods, float(quant)*100.)
                            self.add_record(forecast_date_text, target, target_end_date, location, 'quantile', quant, value)

                        if i_period in de_cum_periods_dict[category]:
                            ip_len = len(cum_periods)
                            for i_rep in range(len(de_cum_periods_dict[category][i_period])):
                                # in case there are fewer repetitions, loop over others:
                                ip_index = i_rep % ip_len
                                de_cum_periods_dict[category][i_period][i_rep] += cum_periods[ip_index]
                        else:
                            de_cum_periods_dict[category][i_period] = [cum_periods[i_rep] for i_rep in range(len(cum_periods))]

        # return 'Germany' summary:
        # there is a separate file with Germany wide deaths. Need to correct by adding the additional deaths here
        #additional_deaths = de_deaths - state_deaths
        additional_counts = 0
        if category == 'death':
            print('Germany total deaths included:',state_sums)
        elif category == 'case':
            print('Germany total cases included:', state_sums)
        #print('Germany total: additional deaths included:', additional_deaths)

        location = self.regional_locations['Germany']
        quants = [0.01, 0.025] + [0.05 + 0.05 * i for i in range(19)] + [0.975, 0.99]
        #for i in range(2):
        for inc_type in inc_types:
            for i_period in de_inc_point_est_dict[category]:
                target = i_period + ' ' + period + ' ahead ' + inc_type + ' ' + category
                days = int(i_period)
                if period == 'wk':
                    days *= 7
                else:
                    days += 1
                target_end_date = (t0 + timedelta(days=(first_sunday + days - 1))).isoformat()

                # point estimates
                value = de_inc_point_est_dict[category][i_period]
                if inc_type == 'cum':
                    value = de_cum_point_est_dict[category][i_period] + additional_counts
                self.add_record(forecast_date_text, target, target_end_date, location, 'point', 'NA', value)

                # quantiles
                if inc_type == 'inc':
                    median = np.percentile(de_inc_periods_dict[category][i_period], 50.)
                    for quant in quants:
                        value = np.percentile(de_inc_periods_dict[category][i_period], float(quant)*100.)
                        # increase the width of the interval for the US as a whole, to account for correlations between states
                        scaled_value = max(0., cor_scale * (value - median) + median)
                        quant_text = '{0:0.3f}'.format(quant)
                        self.add_record(forecast_date_text, target, target_end_date, location, 'quantile', quant_text,
                                        scaled_value)

                elif inc_type == 'cum':
                    median = np.percentile(de_cum_periods_dict[category][i_period], 50.) + additional_counts
                    for quant in quants:
                        value = np.percentile(de_cum_periods_dict[category][i_period], float(quant) * 100.) + additional_counts
                        scaled_value = max(0., cor_scale * (value - median) + median)
                        quant_text = '{0:0.3f}'.format(quant)
                        self.add_record(forecast_date_text, target, target_end_date, location, 'quantile', quant_text,
                                        scaled_value)

        return self.buff



    def add_record(self, forecast_date, target, target_end_date, location, record_type, quantile_text, value):
        record = [forecast_date, target, target_end_date, location, record_type, quantile_text]
        record.append('{0:0.1f}'.format(value))
        self.buff.append(record)

my_forecast = Forecast_hub('/Users/karlen/pypm-temp/germany', ['_2_9_0704'])
# Indicate the total US deaths (up to and including Saturday) here:
#de_deaths = xxx

# 1.5 -> 3.0 on introduction of variant -> 1.5 (back to 1 strain)
my_csv = my_forecast.get_csv(datetime.date(2021, 7, 4), 'case',cor_scale=1.5)
pass
with open('/Users/karlen/pypm-temp/test-germany-forecast-case.csv','w') as out:
    for line in my_csv:
        record = ','.join(line)
        out.write(record + '\n')

my_forecast = Forecast_hub('/Users/karlen/pypm-temp/germany', ['_2_9_0704'])

my_csv = my_forecast.get_csv(datetime.date(2021, 7, 4), 'death',cor_scale=1.5)
pass
with open('/Users/karlen/pypm-temp/test-germany-forecast.csv','w') as out:
    for line in my_csv:
        record = ','.join(line)
        out.write(record + '\n')
