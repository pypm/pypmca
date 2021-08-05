# -*- coding: utf-8 -*-
"""
Produce .csv file suitable for submission to the EU covid19-forecast-hub
(https://github.com/epiforecasts/covid19-forecast-hub-europe)

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
        'Austria': 'AT',
        'Belgium': 'BE',
        'Bulgaria': 'BG',
        'Croatia': 'HR',
        'Cyprus': 'CY',
        'Czechia': 'CZ',
        'Denmark': 'DK',
        'Estonia': 'EE',
        'Finland': 'FI',
        'France': 'FR',
        'Germany': 'DE',
        'Greece': 'GR',
        'Hungary': 'HU',
        'Iceland': 'IS',
        'Ireland': 'IE',
        'Italy': 'IT',
        'Latvia': 'LV',
        'Liechtenstein': 'LI',
        'Lithuania': 'LT',
        'Luxembourg': 'LU',
        'Malta': 'MT',
        'Netherlands': 'NL',
        'Norway': 'NO',
        'Poland': 'PL',
        'Portugal': 'PT',
        'Romania': 'RO',
        'Slovakia': 'SK',
        'Slovenia': 'SI',
        'Spain': 'ES',
        'Sweden': 'SE',
        'Switzerland': 'CH',
        'United Kingdom': 'GB'
        }

        self.hospitalization_states = ['BE','HR','DK','EE','FR','NL','NO','SI','ES','GB','IE']

    def get_data(self):
        data_folder = 'data/covid19/EU'

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

    def get_csv(self, forecast_date):
        t0 = datetime.date(2020, 3, 1)
        forecast_date_text = forecast_date.isoformat()
        day_of_week = forecast_date.weekday()
        days_after_t0 = (forecast_date - t0).days
        first_sunday = days_after_t0
        if day_of_week == 0:
            first_sunday -= 1
        elif day_of_week < 6:
            first_sunday += 6 - day_of_week

        period = 'wk'
        inc_type = 'inc'

        for state in self.regional_abbreviations:
            abbrev = self.regional_abbreviations[state].lower()
            location = self.regional_abbreviations[state]

            categories = ['case','death']
            if location in self.hospitalization_states:
                categories = ['case', 'death', 'hospitalization']

            for category in categories:

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
                elif category == 'hospitalization':
                    for model_name in self.model_names:
                        for model_suffix in ['_h','']:
                            try:
                                filename = abbrev + model_name + model_suffix + '.pypm'
                                path_model = self.model_dir / filename
                                model = Model.open_file(path_model)
                                success = True
                                break
                            except:
                                pass

                if not success:
                    print('No model for: ', abbrev)
                else:
                    print(category,model.name)

                    hub_dict = model.user_dict['forecast_hub'][forecast_date][category]
                    point_est_dict = hub_dict['point_estimates']
                    quantile_dict = hub_dict['quantiles']

                    for i_period in point_est_dict:
                        # eu has different label for hospitalization...
                        label = category
                        if category == 'hospitalization':
                            label = 'hosp'
                        target = i_period+' '+period+' ahead '+inc_type+' '+label
                        days = int(i_period)
                        if period == 'wk':
                            days *= 7
                        else:
                            days += 1
                        target_end_date = (t0 + timedelta(days=(first_sunday + days - 1))).isoformat()

                        # point estimates
                        value = point_est_dict[i_period]
                        self.add_record(forecast_date_text,target,target_end_date,location,'point','NA',value)

                        # quantiles
                        for quant in quantile_dict[i_period]:
                            value = quantile_dict[i_period][quant]
                            self.add_record(forecast_date_text, target, target_end_date, location, 'quantile', quant, value)

        return self.buff

    def add_record(self, forecast_date, target, target_end_date, location, record_type, quantile_text, value):
        record = [forecast_date, target, target_end_date, location, record_type, quantile_text]
        record.append('{0:0.0f}'.format(value))
        self.buff.append(record)

my_forecast = Forecast_hub('/Users/karlen/pypm-temp/eu', ['_2_9_0801'])

my_csv = my_forecast.get_csv(datetime.date(2021, 8, 1))
pass
with open('/Users/karlen/pypm-temp/test-eu-forecast.csv','w') as out:
    for line in my_csv:
        record = ','.join(line)
        out.write(record + '\n')
