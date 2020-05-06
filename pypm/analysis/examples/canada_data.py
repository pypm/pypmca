# -*- coding: utf-8 -*-
"""
Example of data description file: data.py

@author: karlen
"""

def get_data_description():
    """ Define data provided by Virihealth
    """
    data = {}
    data['nation'] = 'Canada'
    data['description'] = 'Canada-wide data by Province and Territory'
    data['source'] = 'Virihealth'
    data['source_url'] = 'www.Virihealth.com'

    # common regional abbreviations used in the data files
    regional_abbreviations = {
        'BC':'BC',
        'Alberta':'AB',
        'Saskatchewan':'SK',
        'Manitoba':'MB',
        'Ontario':'ON',
        'Quebec':'QC',
        'New Brunswick':'NB',
        'Newfoundland':'NL',
        'Nova Scotia':'NS',
        'PEI':'PE',
        'Yukon':'YT',
        'NWT':'NT',
        'Nunavut':'NU'
        }

    files_data = {}
    filenames = ['CV - TR.csv', 'CV - H.csv']
    # since the Canada file formats are in common
    # we provide the same information for each file, using a for loop
    for filename in filenames:
        file_data = {}
        file_data['source'] = 'More detailed information if folder has multiple sources'
        file_data['date header'] = 'Date'
        file_data['date start'] = [2020, 3, 1]

        files_data[filename] = file_data
    data['files'] = files_data

    #overall country combined data

    populations_data = {}
    # First for the CV - TR.csv file:
    tr_populations = ['reported', 'deaths']
    filename = filenames[0]
    for population in tr_populations:
        pop_data_daily = {}
        pop_data_daily['filename'] = filename
        header = ''
        if population == 'reported':
            header = 'Cases'
        if population == 'deaths':
            header = 'Deaths'
        pop_data_daily['header'] = header

        pop_data_total = {}
        pop_data_total['filename'] = filename
        header = ''
        if population == 'reported':
            header = 'Cumulative Cases'
        if population == 'deaths':
            header = 'Cumulative Deaths'
        pop_data_total['header'] = header

        population_data = {'daily': pop_data_daily, 'total': pop_data_total}
        populations_data[population] = population_data

    # Now treat the CV - H.csv file
    h_populations = ['in_hospital', 'in_icu']
    filename = filenames[1]
    for population in h_populations:
        pop_data_total = {}
        pop_data_total['filename'] = filename
        header = ''
        if population == 'in_hospital':
            header = 'TotH'
        if population == 'in_icu':
            header = 'TotI'
        pop_data_total['header'] = header

        population_data = {'total': pop_data_total}
        populations_data[population] = population_data

    data['national_data'] = populations_data

    #for each region, define location of each population data (file/column header)
    #either daily or total (or both) data can be provided

    regions_data = {}
    for region in regional_abbreviations:

        populations_data = {}
        # First for the CV - TR.csv file:
        tr_populations = ['reported', 'deaths']
        filename = filenames[0]
        for population in tr_populations:
            pop_data_daily = {}
            pop_data_daily['filename'] = filename
            header = ''
            if population == 'reported':
                header = regional_abbreviations[region]
            if population == 'deaths':
                header = regional_abbreviations[region]+'.1'
            pop_data_daily['header'] = header

            pop_data_total = {}
            pop_data_total['filename'] = filename
            header = ''
            if population == 'reported':
                header = regional_abbreviations[region]+' C'
            if population == 'deaths':
                header = regional_abbreviations[region]+'C.1'
            pop_data_total['header'] = header

            population_data = {'daily': pop_data_daily, 'total': pop_data_total}
            populations_data[population] = population_data

        # Now treat the CV - H.csv file
        h_populations = ['in_hospital', 'in_icu']
        filename = filenames[1]
        for population in h_populations:
            pop_data_total = {}
            pop_data_total['filename'] = filename
            header = ''
            if population == 'in_hospital':
                header = regional_abbreviations[region]+'H'
            if population == 'in_icu':
                header = regional_abbreviations[region]+'I'
            pop_data_total['header'] = header

            population_data = {'total': pop_data_total}
            populations_data[population] = population_data

        regions_data[region] = populations_data

    data['regional_data'] = regions_data

    return data
