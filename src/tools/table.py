# -*- coding: utf-8 -*-
"""
Show the details of a pypm model by a set of tables for:
    - connectors
    - populations
    - parameters (fixed and variable)
    - variable parameters - not yet done
    - transitions: modifiers and injectors separately

@author: karlen
"""
from texttable import Texttable

from Model import Model
from Parameter import Parameter

def print_all(model):
    print(connector_table(model))
    print('\n')
    print(population_table(model))
    print('\n')
    print(parameter_table(model))
    print('\n')
    print(modifier_table(model))
    print('\n')
    print(injector_table(model))    

def table_setup(model, width):
    if not isinstance(model, Model):
        raise TypeError('Error in model_table. '+
                            ': model argument must be a Model object')

# show table of modifiers

    table = Texttable()
    table.set_deco(Texttable.HEADER | Texttable.HLINES)
    table.set_max_width(width)
    return table

def injector_table(model, width=120):    
    table = table_setup(model, width)
    
 #   transition_name, time_spec, transition_time,
 #                to_population, injection

    header = ['Injector', 'Time parameter', 'Time', 'spec', 'to (Population)', 
              'Injection parameter', 'value']
    table.set_cols_dtype(['t', 't', 'a', 't', 't', 't', 'a'])
  
    table.set_header_align(['l', 'l', 'l', 'l', 'l', 'l', 'l'])
    
    rows = [header]
    for key in model.transitions:
        trans = model.transitions[key]
        if type(trans).__name__ == 'Injector':
            row = [str(trans), str(trans.transition_time),
                   trans.transition_time.get_value(), trans.time_spec,
                   str(trans.to_population),
                   str(trans.injection), trans.injection.get_value()]
            rows.append(row)
        
    table.add_rows(rows)
    return table.draw()

def modifier_table(model, width=120):
    table = table_setup(model, width)

    header = ['Modifier', 'Time parameter', 'Time', 'spec', 'Parameter', 
              'before value', 'after parameter', 'after value']
    table.set_cols_dtype(['t', 't', 'a', 't', 't', 'a', 't', 'a'])
  
    table.set_header_align(['l', 'l', 'l', 'l', 'l', 'l', 'l', 'l'])
    
    rows = [header]
    for key in model.transitions:
        trans = model.transitions[key]
        if type(trans).__name__ == 'Modifier':
            row = [str(trans), str(trans.transition_time),
                   trans.transition_time.get_value(), trans.time_spec,
                   str(trans.parameter_before), trans.parameter_before.get_value(),
                   str(trans.parameter_after), trans.parameter_after.get_value()]
            rows.append(row)
        
    table.add_rows(rows)
    return table.draw()

def get_connector_row(connector):
    con = connector
    row = [str(con), type(con).__name__]
    for pop in [con.from_population, con.to_population]:
        buff = []
        if isinstance(pop, list):
            for item in pop:
                buff.append(str(item))
        else:
            buff.append(str(pop))
        row.append('\n'.join(buff))
    names = []
    types = []
    values = []
    status = []
    for par_key in con.parameters:
        par = con.parameters[par_key]
        names.append(str(par))
        types.append(par.parameter_type)
        if par.parameter_type in ['int', 'bool']:
            values.append(str(par.get_value()))
        else:
            values.append("{:10.3f}".format(par.get_value()))
        status.append(par.get_status())
    row.append('\n'.join(names))
    row.append('\n'.join(values))
    return row
    

def connector_table(model, width=120):
    table = table_setup(model, width)

    header = ['Connector','Type','From\n(populations)','To\n(populations)',
              'Parameters','Values']
    table.set_cols_dtype(['t', 't', 't', 't',
                          't', 't'])
    table.set_header_align(['l', 'l', 'l', 'l',
                            'l', 'r'])
    
    rows = [header]
    for key in model.connectors:
        con = model.connectors[key]
        rows.append(get_connector_row(con))
        if type(con).__name__ == 'Chain':
            for chain_con in con.chain:
                row = get_connector_row(chain_con)
                row[0] = str(con)+': '+row[0]
                row[1] = '->'+row[1]
                rows.append(row)
    
    table.add_rows(rows)
    return table.draw()

def population_table(model, width=120):
    table = table_setup(model, width)

    header = ['Population', 'Description', 'Parameter', 'Start value']
    table.set_cols_dtype(['t', 't', 't', 'a'])
    table.set_header_align(['l', 'l', 'l', 'l'])
    
    rows = [header]
    
    for key in model.populations:
        pop = model.populations[key]
        row = [str(pop), pop.description]
        init = pop.initial_value
        if isinstance(init, Parameter):
            row.append(str(init))
            row.append(init.get_value())
        else:
            row.append(' ')
            row.append(init)
        
        rows.append(row)
        
    table.add_rows(rows)
        
    return table.draw()

def get_par_row(par):
    buff = []
    buff.append(str(par))
    buff.append(par.description)
    buff.append(par.get_min())
    buff.append(par.get_max())
    buff.append(par.initial_value)
    buff.append(par.get_value())
    return buff


def parameter_table(model, width=120):
    table = table_setup(model, width)

    header = ['Parameter', 'Description', 'min', 'max', 'Init', 'Value']
    table.set_cols_dtype(['t', 't', 'a', 'a', 'a', 'a'])
    table.set_header_align(['l', 'l', 'l', 'l', 'l', 'l'])

    par_dict = {}
    
    for key in model.connectors:
        con = model.connectors[key]
        for par_key in con.parameters:
            par = con.parameters[par_key]
            row = get_par_row(par)
            if len(row) > 0:
                par_dict[row[0]] = row
        if type(con).__name__ == 'Chain':
            for chain_con in con.chain:
                for par_key in chain_con.parameters:
                    par = chain_con.parameters[par_key]
                    row = get_par_row(par)
                    if len(row) > 0:
                        par_dict[row[0]] = row

    for key in model.populations:
        pop = model.populations[key]
        init = pop.initial_value
        if isinstance(init, Parameter):
            row = get_par_row(init)
            par_dict[row[0]] = row
        
    for key in model.transitions:
        trans = model.transitions[key]
        for par_key in trans.parameters:
            par = trans.parameters[par_key]
            row = get_par_row(par)
            if len(row) > 0:
                par_dict[row[0]] = row

    par_list = list(par_dict.keys())
    par_list.sort()
    
    rows = [header]
    for item in par_list:
        rows.append(par_dict[item])
        
    table.add_rows(rows)
        
    return table.draw()


#my_model = Model.open_file('../test/model2v1.pypm')
#print(parameter_table(my_model))
#print_all(my_model)