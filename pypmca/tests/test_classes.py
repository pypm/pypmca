#!/usr/bin/env python

"""Tests for `pypmca` package."""

import pytest

from pypmca import Model, Population, Delay, Parameter, Multiplier, Propagator, \
    Splitter, Adder, Subtractor, Chain, Modifier, Injector
import numpy as np
from pathlib import Path

example_dir = Path('../../examples/').resolve()
path_model_2_8 = example_dir / 'ref_model_2_8.pypm'


def test_class_Model():
    """tests to ensure the behaviour class Model"""
    test_model = Model('test')


def test_class_Population():
    """tests to ensure the behaviour class Population"""
    init_value = 100
    test_pop = Population('test population', init_value, description='For testing populations',
                          hidden=True, color='black', show_sim=False, report_noise=False,
                          report_noise_par=None)
    assert test_pop.history[0] == init_value

    model = Model.open_file(path_model_2_8)
    test_pop.set_model(model)

    incoming = 10
    test_pop.update_future_fast(incoming)
    assert test_pop.future[0] == incoming

    test_pop.do_time_step()
    assert test_pop.history[1] == init_value + incoming
    assert len(test_pop.future) == 0 or test_pop.future == 0

    scale_factor = 0.5
    test_pop.scale_history(scale_factor)
    assert test_pop.history[0] == init_value * scale_factor
    assert test_pop.history[1] == (init_value + incoming) * scale_factor

    # check that noise in reporting makes sense
    # Expectations are not effected - data should be
    # but the total reported should be the same after all reporting complete
    noise_factor = Parameter('noise_factor', 0.3)
    backlog_factor = Parameter('backlog_factor', 0.8)
    # restart - back to initial value
    for expectations in [True, False]:
        test_pop.reset()
        future = [0, 10, 100, 100, 100, 10, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
        future_sum = np.sum(np.array(future))
        test_pop.set_report_noise(True, noise_factor, backlog_factor, None)
        test_pop.future = future
        for i in range(len(future) + 5):
            test_pop.do_time_step(expectations=expectations)
        assert test_pop.history[-1] == init_value + future_sum

    report_days = Parameter('report_days',127,parameter_min=-7,parameter_max=127,parameter_type='int')
    for report_day_value in [127,63,-1,-2,-5,-7]:
        report_days.set_value(report_day_value)
        test_pop.reset()
        future = [0, 10, 100, 100, 100, 10, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
        future_sum = np.sum(np.array(future))
        test_pop.set_report_noise(True, noise_factor, backlog_factor, report_days)
        test_pop.future = future
        for i in range(len(future) + 5):
            test_pop.do_time_step(expectations=expectations)
        assert test_pop.history[-1] == init_value + future_sum


def test_class_Parameter():
    """tests to ensure the behaviour class Parameter"""
    test_parameter = Parameter('test_par', 10., parameter_min=5., parameter_max=100.)
    for value in [50]:
        error_caught = False
        try:
            test_parameter.set_value(value)
        except ValueError:
            error_caught = True
        except TypeError:
            error_caught = True
        assert error_caught


def test_class_Delay():
    """tests to ensure the behaviour class Delay"""
    test_model = Model('test_model')
    EPS = 0.1
    mean = 10.
    std_dev = 4.
    half_width = float(std_dev * np.sqrt(12.) / 2.)
    k_vals = [1, 2, 3]
    for delay_type in ['norm', 'uniform', 'erlang', 'gamma']:
        k_s = [1]
        if delay_type == 'erlang':
            k_s = k_vals
        for k in k_s:
            delay_pars = {
                'mean': Parameter('mean', mean, parameter_min=-100., parameter_max=100.),
                'sigma': Parameter('sigma', std_dev, parameter_min=-100., parameter_max=100.),
                'half_width': Parameter('hw', half_width, parameter_min=-100., parameter_max=100.),
                'k': Parameter('k', k, parameter_type='int', parameter_min=1, parameter_max=100)
            }
            for time_step in [1., 1. / 4.]:
                test_model.set_time_step(time_step)
                # The delay is created after the model: since it is not associated with a connector, the model
                # does not know to call its update method. (The model does not have stand alone delays.)
                test_delay = Delay('test_delay', delay_type, delay_parameters=delay_pars, model=test_model)
                distribution = test_delay.future_expectations
                sum_p = 0.
                sum_tp = 0.
                sum_ttp = 0.
                for i in range(len(distribution)):
                    sum_p += distribution[i]
                    sum_tp += i * distribution[i] * time_step
                    sum_ttp += i * i * distribution[i] * time_step ** 2
                assert np.abs(sum_p - 1.) < EPS
                est_mean = sum_tp
                assert np.abs(est_mean - mean) < EPS
                est_sigma = np.sqrt(sum_ttp - sum_tp * sum_tp)
                if delay_type != 'erlang':
                    assert np.abs(est_sigma - std_dev) < EPS
                else:
                    assert np.abs(est_sigma - mean / np.sqrt(1. * k)) < EPS


def test_class_Adder():
    """tests to ensure the behaviour class Adder"""
    for itest in range(5):
        f_over_t = 2
        to_init = 10
        from_init = f_over_t * to_init
        from_pop = Population('from_pop', from_init)
        from_next = 40
        from_pop.future = [from_next]
        to_pop = Population('to_pop', to_init)
        if itest == 0:
            test_adder = Adder('test_add', from_pop, to_pop)
            test_adder.update_expectation()
            assert to_pop.future[0] == from_next
        elif itest == 1:
            sf = 2.
            scale_factor = Parameter('scale_factor', sf, 0., 10.)
            test_adder = Adder('test_add', from_pop, to_pop, scale_factor=scale_factor)
            test_adder.update_expectation()
            assert to_pop.future[0] == from_next*sf
        elif itest == 2:
            sf = 3.
            ratio_pops = [from_pop, to_pop]
            scale_factor = Parameter('scale_factor', sf, 0., 10.)
            test_adder = Adder('test_add', from_pop, to_pop, scale_factor=scale_factor, ratio_populations=ratio_pops)
            test_adder.update_expectation()
            assert to_pop.future[0] == from_next*sf*f_over_t
        elif itest == 3:
            sf = 3
            ratio_pops = [from_pop, to_pop]
            scale_factor = Parameter('scale_factor', sf, 0, 10,' ','int')
            test_adder = Adder('test_add', from_pop, to_pop, scale_factor=scale_factor, ratio_populations=ratio_pops)
            test_adder.update_data()
            assert to_pop.future[0] == from_next*sf*f_over_t
        elif itest == 4:
            sf = 3.1
            ratio_pops = [from_pop, to_pop]
            scale_factor = Parameter('scale_factor', sf, 0., 10.)
            test_adder = Adder('test_add', from_pop, to_pop, scale_factor=scale_factor, ratio_populations=ratio_pops)
            test_adder.update_data()
            assert to_pop.future[0] >= from_next*int(sf)*f_over_t

def test_class_Injector():
    """tests to ensure the behaviour class Injector"""
    test_model = Model('test_model')
    number = 50.
    inject = Parameter('inject', number, parameter_min=0., parameter_max=1000.)
    time = 5
    trans_time = Parameter('time', time, parameter_type='int', parameter_min=0, parameter_max=1000)
    to_pop = Population('to_pop', 0)
    test_injector = Injector('injector', 'rel_days', trans_time, to_pop, inject, model=test_model)
    test_model.add_transition(test_injector)
    for time_step in [1., 1. / 4.]:
        test_model.set_time_step(time_step)
        to_pop.reset()
        test_injector.take_action()
        assert to_pop.future[0] == number
        assert np.abs(test_injector.trigger_step - time / time_step) < 0.1


def test_class_Modifier():
    """tests to ensure the behaviour class Modifier"""
    test_model = Model('test_model')
    mod_time = Parameter('time', 5, parameter_type='int', parameter_min=0, parameter_max=1000)
    par_val = 0.3
    par_0_val = 0.5
    par_1_val = 0.7
    parameter = Parameter('par', par_val)
    parameter_0 = Parameter('par_0', par_0_val)
    parameter_1 = Parameter('par_1', par_1_val)
    test_modifier = Modifier('test_modifier', 'rel_days', mod_time, parameter, parameter_0, parameter_1,
                             model=test_model)
    test_model.add_transition(test_modifier)
    for time_step in [1., 1. / 4.]:
        test_model.set_time_step(time_step)
        parameter.reset()
        assert parameter.get_value() == par_val
        test_modifier.take_action()
        assert parameter.get_value() == par_1_val
        test_modifier.reset()
        assert parameter.get_value() == par_0_val
        assert np.abs(test_modifier.trigger_step - 5 / time_step) < 0.1


def test_class_Multiplier():
    """tests to ensure the behaviour class Multiplier"""
    test_model = Model('test_model')
    EPS = 1.
    n1 = 50.
    n2 = 20.
    n3 = 2.
    scale = 0.1
    f_pops = [Population('f1_pop', n1), Population('f2_pop', n2), Population('f3_pop', n3)]
    to_pop = Population('to_pop', 0.)
    scale_par = Parameter('alpha', scale)
    delay = Delay('fast', 'fast')
    test_multiplier = Multiplier('test_multiplier', f_pops, to_pop, scale_par, delay, model=test_model)
    test_model.add_connector(test_multiplier)
    for time_step in [1., 1. / 4.]:
        test_model.set_time_step(time_step)
        # expectation:
        expected = n1 * n2 / n3 * scale * time_step
        to_pop.reset()
        test_multiplier.set_distribution('poisson', None)
        test_multiplier.update_expectation()
        assert to_pop.future[0] == expected

        # Poisson
        n_rep = 1000
        n_list = []
        for i in range(n_rep):
            to_pop.reset()
            test_multiplier.update_data()
            n_list.append(to_pop.future[0])
        assert np.abs(np.mean(n_list) - expected) < EPS
        assert np.abs(np.std(n_list) - np.sqrt(expected)) < EPS

        # Negative binomial
        p_nb = 0.2
        nbinom_par = Parameter('nb', p_nb)
        test_multiplier.set_distribution('nbinom', nbinom_par)

        n_rep = 1000
        n_list = []
        for i in range(n_rep):
            to_pop.reset()
            test_multiplier.update_data()
            n_list.append(to_pop.future[0])
        assert np.abs(np.mean(n_list) - expected) < EPS
        assert np.abs(np.std(n_list) - np.sqrt(expected / p_nb)) < EPS


def test_class_Propagator():
    """tests to ensure the behaviour class Propagator"""
    test_model = Model('test_model')
    start = 100000
    from_pop = Population('from_pop', 0)
    from_pop.future = [start]
    to_pop = Population('to_pop', 0.)
    frac = 0.4
    fraction = Parameter('frac', frac)
    mean = 10.
    std_dev = 4.
    delay_pars = {
        'mean': Parameter('mean', mean, parameter_min=-100., parameter_max=100.),
        'sigma': Parameter('sigma', std_dev, parameter_min=-100., parameter_max=100.)
    }
    test_delay = Delay('test_delay', 'norm', delay_parameters=delay_pars, model=test_model)

    test_prop = Propagator('test_prop', from_pop, to_pop, fraction, test_delay)
    test_model.add_connector(test_prop)
    for time_step in [1., 1. / 4.]:
        test_model.set_time_step(time_step)
        for func in [test_prop.update_expectation, test_prop.update_data]:
            EPS = 0.01
            if func == test_prop.update_data:
                EPS = 0.1
            to_pop.reset()
            func()
            distribution = to_pop.future
            sum_p = 0.
            sum_tp = 0.
            sum_ttp = 0.
            for i in range(len(distribution)):
                sum_p += distribution[i]
                sum_tp += i * distribution[i] * time_step
                sum_ttp += i * i * distribution[i] * time_step ** 2
            assert np.abs(sum_p - start * frac) < EPS * start * frac
            est_mean = sum_tp / (start * frac)
            assert np.abs(est_mean - mean) < EPS * mean
            est_sigma = np.sqrt(sum_ttp / (start * frac) - est_mean ** 2)
            assert np.abs(est_sigma - std_dev) < EPS * std_dev


def test_class_Splitter():
    """tests to ensure the behaviour class Splitter"""
    test_model = Model('test_model')
    start = 100000
    from_pop = Population('from_pop', 0)
    from_pop.future = [start]
    to_pops = [Population('to_pop1', 0.), Population('to_pop2', 0.)]
    fracs = [0.4, 0.6]
    fraction = Parameter('frac', fracs[0])
    mean = 10.
    std_dev = 4.
    delay_pars = {
        'mean': Parameter('mean', mean, parameter_min=-100., parameter_max=100.),
        'sigma': Parameter('sigma', std_dev, parameter_min=-100., parameter_max=100.)
    }
    test_delay = Delay('test_delay', 'norm', delay_parameters=delay_pars, model=test_model)

    test_split = Splitter('test_prop', from_pop, to_pops, [fraction], test_delay)
    test_model.add_connector(test_split)
    for time_step in [1., 1. / 4.]:
        test_model.set_time_step(time_step)
        for func in [test_split.update_expectation, test_split.update_data]:
            EPS = 0.01
            if func == test_split.update_data:
                EPS = 0.1
            to_pops[0].reset()
            to_pops[1].reset()
            func()
            total = 0
            for i in range(2):
                distribution = to_pops[i].future
                frac = fracs[i]
                sum_p = 0.
                sum_tp = 0.
                sum_ttp = 0.
                for i in range(len(distribution)):
                    sum_p += distribution[i]
                    sum_tp += i * distribution[i] * time_step
                    sum_ttp += i * i * distribution[i] * time_step ** 2
                total += sum_p
                assert np.abs(sum_p - start * frac) < EPS * start * frac
                est_mean = sum_tp / (start * frac)
                assert np.abs(est_mean - mean) < EPS * mean
                est_sigma = np.sqrt(sum_ttp / (start * frac) - est_mean ** 2)
                assert np.abs(est_sigma - std_dev) < EPS * std_dev
            assert np.abs(total - start) < 0.1


def test_class_Chain():
    """tests to ensure the behaviour class Chain"""
    test_model = Model('test_model')
    from_pop = Population('from_pop', 0.)
    to_pop = Population('to_pop', 0.)
    pop_a = Population('pop_a', 0.)
    pop_b = Population('pop_b', 0.)
    start = 100000
    from_pop.future = [start]

    mean = 10.
    sigma = 2.
    delay_pars = {'mean': Parameter('mean', mean, parameter_min=0.1, parameter_max=100.),
                  'sigma': Parameter('sigma', sigma, parameter_min=0.1, parameter_max=1000.)}
    delay = Delay('delay', 'norm', delay_pars, test_model)
    frac = 0.8
    fraction = Parameter('frac', frac)
    chain = []
    chain.append(Propagator('prop_0', from_pop, pop_a, fraction, delay))
    chain.append(Propagator('prop_1', pop_a, pop_b, fraction, delay))

    test_chain = Chain('test_chain', from_pop, to_pop, chain, fraction, delay, test_model)
    test_model.add_connector(test_chain)

    for time_step in [1., 1. / 4.]:
        test_model.set_time_step(time_step)
        for func in [test_chain.update_expectation, test_chain.update_data]:
            EPS = 0.02
            if func == test_chain.update_data:
                EPS = 0.2
            to_pop.reset()
            pop_a.reset()
            pop_b.reset()
            func()
            for pop in [to_pop, pop_b]:
                distribution = pop.future
                total = start * (1. - frac ** 2) * frac
                ave = mean
                std_dev = sigma
                if pop == pop_b:
                    total = start * frac ** 2
                    ave = 2. * mean
                    std_dev = np.sqrt(2.) * sigma
                sum_p = 0.
                sum_tp = 0.
                sum_ttp = 0.
                for i in range(len(distribution)):
                    sum_p += distribution[i]
                    sum_tp += i * distribution[i] * time_step
                    sum_ttp += i * i * distribution[i] * time_step ** 2
                assert np.abs(sum_p - total) < EPS * total
                est_mean = sum_tp / total
                assert np.abs(est_mean - ave) < EPS * ave
                est_sigma = np.sqrt(sum_ttp / total - est_mean ** 2)
                assert np.abs(est_sigma - std_dev) < EPS * std_dev


def test_class_Subtractor():
    """tests to ensure the behaviour class Subtractor"""
    from_init = 20
    to_init = 10
    from_pop = Population('from_pop', from_init)
    to_pop = Population('to_pop', to_init)
    to_next = 40
    to_pop.future = [to_next]
    test_subtractor = Subtractor('test_sub', from_pop, to_pop)
    test_subtractor.update_expectation()
    assert from_pop.future[0] == -1 * to_next
