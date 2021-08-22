#!/usr/bin/env python

"""Tests for `pypmca` package."""

import pytest

from pypmca import Model, Population, Delay, Parameter, Multiplier, Propagator, \
    Splitter, Adder, Subtractor, Chain, Modifier, Injector, Ensemble
from pypmca.analysis.Optimizer import Optimizer
from pypmca.analysis.Trajectory import Trajectory
from pypmca.tools.IntervalMaker import IntervalMaker
import numpy as np
import copy
from pathlib import Path
import datetime

example_dir = Path('../../examples/').resolve()
path_model_2_9 = example_dir / 'ref_model_2_9.pypm'
path_model_2_8 = example_dir / 'ref_model_2_8.pypm'
path_model_2_7 = example_dir / 'ref_model_2_7.pypm'
path_model_2_6 = example_dir / 'ref_model_2_6.pypm'
path_model_2_5 = example_dir / 'ref_model_2_5.pypm'
path_model_2_4 = example_dir / 'ref_model_2_4.pypm'
path_model_2_3 = example_dir / 'ref_model_2_3.pypm'
path_model_2_2 = example_dir / 'ref_model_2_2.pypm'
path_model_2 = example_dir / 'ref_model_2.pypm'
path_model_1 = example_dir / 'ref_model_1.pypm'

def test_Model_copy_values_from():
    ref_model_2 = Model.open_file(path_model_2)
    alberta_model_2 = Model.open_file('ab_2_0514.pypm')

    ref_model_2.copy_values_from(alberta_model_2)

    ref_model_2.reset()
    ref_model_2.evolve_expectations(80)
    alberta_model_2.reset()
    alberta_model_2.evolve_expectations(80)

    for pop_name in ref_model_2.populations:
        pop = ref_model_2.populations[pop_name]
        # Bug fixed in ref_model_2, not in alberta_model_2 accounts for small differences in icu
        if 'vent' not in pop.name and 'in_' not in pop.name:
            ref_hist = pop.history
            alberta_hist = alberta_model_2.populations[pop_name].history
            for i in range(len(ref_hist)):
                if np.abs(ref_hist[i] - alberta_hist[i]) > 0.0001:
                    print(pop.name, i, ref_hist[i], alberta_hist[i], np.abs(ref_hist[i] - alberta_hist[i]))
                assert np.abs(ref_hist[i] - alberta_hist[i]) < 0.001

def test_Model_properties():
    """tests to ensure the properties of Model"""
    ref_model_2 = Model.open_file(path_model_2)
    ref_model_1 = Model.open_file(path_model_1)

    for test_model in [ref_model_1, ref_model_2]:
        # check simple scaling of initial contagious population
        # Will not be exact due to bootstrap
        EPS = 0.01
        cont_0 = test_model.parameters['cont_0'].get_value()
        test_model.reset()
        test_model.evolve_expectations(200)
        n_reported1 = test_model.populations['reported'].history[-1]
        assert test_model.populations['contagious'].history[-1] < 10.

        test_model.parameters['cont_0'].set_value(2 * cont_0)
        test_model.reset()
        test_model.evolve_expectations(200)
        n_reported2 = test_model.populations['reported'].history[-1]
        assert np.abs(n_reported2 - n_reported1 * 2) < EPS * n_reported1 * 2

    # check that the mean of many data runs is near expectation

    ref_model = Model.open_file(path_model_2)
    ref_model.reset()
    n_days = 60
    ref_model.evolve_expectations(n_days)
    sim_models = []
    n_rep = 100
    for i in range(n_rep):
        sim_model = copy.deepcopy(ref_model)
        sim_model.reset()
        sim_model.generate_data(n_days)
        sim_models.append(sim_model)
    for pop_name in ref_model.populations:
        pop = ref_model.populations[pop_name]
        if pop.show_sim:
            results = []
            for sim_model in sim_models:
                results.append(sim_model.populations[pop_name].history[-1])
            mean = np.mean(np.array(results))
            std = np.std(np.array(results))
            error = std / np.sqrt(1. * n_rep)
            expect = ref_model.populations[pop_name].history[-1]
            assert np.abs(expect - mean) < 8. * error


def test_Ensemble_properties_identical():
    """tests to ensure the properties of Ensemble"""
    # Test that the ensemble of two identical models
    # behalves like twice a single model.
    # independent of the contact matrix
    test_a = Model.open_file(path_model_2)
    test_a.name = 'test_a'
    test_b = Model.open_file(path_model_2)
    test_b.name = 'test_b'
    reference = Model.open_file(path_model_2)
    reference.name = 'reference'
    single = Model.open_file(path_model_2)
    single.name = 'single'

    test_ensemble = Ensemble('test_ensemble', reference)
    test_ensemble.upload_models([test_a, test_b])

    contacts = [
        [[1., 0.0], [0.0, 1.]],
        [[1., 1.0], [1.0, 1.]],
        [[1., 0.5], [0.7, 1.]],
    ]

    for contact in contacts:
        test_ensemble.define_cross_transmission('infection cycle', 'infected',
                                                'susceptible', 'total',
                                                'contagious', 'alpha',
                                                contact_type='fixed',
                                                contact=contact)

        n_days = 100
        test_ensemble.reset()
        test_ensemble.evolve_expectations(n_days)

        single.reset()
        single.evolve_expectations(n_days)

        for pop_name in test_ensemble.populations:
            pop = test_ensemble.populations[pop_name]
            if pop.show_sim:
                ens_hist = test_ensemble.populations[pop_name].history
                single_hist = single.populations[pop_name].history
                for i in range(len(ens_hist)):
                    ratio = ens_hist[i] / single_hist[i]
                    assert np.abs(ratio - 2.) < 0.001

def test_Ensemble_properties_equality():
    """tests to ensure the properties of Ensemble"""
    # Test that the ensemble of two identical models, with different scales
    # and equality contact matrix behaves like the same ensemble
    # with independent contact matrix
    test_a1 = Model.open_file(path_model_2_3)
    test_a1.name = 'test_a1'
    test_b1 = Model.open_file(path_model_2_3)
    test_b1.name = 'test_b1'
    reference1 = Model.open_file(path_model_2_3)
    reference1.name = 'reference1'
    cont_0 = test_b1.parameters['cont_0'].get_value()
    test_b1.parameters['cont_0'].set_value(cont_0*2.)
    N_0 = test_b1.parameters['N_0'].get_value()
    test_b1.parameters['N_0'].set_value(N_0*2)

    test_ensemble1 = Ensemble('test_ensemble1', reference1)
    test_ensemble1.upload_models([test_a1, test_b1])

    test_ensemble1.define_cross_transmission('infection cycle', 'infected',
                                            'susceptible', 'total',
                                            'contagious', 'alpha',
                                            contact_type='diagonal')

    test_a2 = Model.open_file(path_model_2_3)
    test_a2.name = 'test_a2'
    test_b2 = Model.open_file(path_model_2_3)
    test_b2.name = 'test_b2'
    reference2 = Model.open_file(path_model_2_3)
    reference2.name = 'reference2'
    cont_0 = test_b2.parameters['cont_0'].get_value()
    test_b2.parameters['cont_0'].set_value(cont_0*2.)
    N_0 = test_b2.parameters['N_0'].get_value()
    test_b2.parameters['N_0'].set_value(N_0*2)

    test_ensemble2 = Ensemble('test_ensemble2', reference2)
    test_ensemble2.upload_models([test_a2, test_b2])

    test_ensemble2.define_cross_transmission('infection cycle', 'infected',
                                            'susceptible', 'total',
                                            'contagious', 'alpha',
                                            contact_type='equality')

    n_days = 100

    test_ensemble1.reset()
    test_ensemble1.evolve_expectations(n_days)

    test_ensemble2.reset()
    test_ensemble2.evolve_expectations(n_days)

    for pop_name in test_ensemble1.populations:
        pop1 = test_ensemble1.populations[pop_name]
        if pop1.show_sim:
            pop1h = test_ensemble1.populations[pop_name].history
            pop2h = test_ensemble2.populations[pop_name].history
            for i in range(len(pop1h)):
                if pop2h[i] > 0.:
                    ratio = pop1h[i]/pop2h[i]
                    assert np.abs(ratio - 1.) < 0.03

def test_Ensemble_properties_different():
    """tests to ensure the properties of Ensemble with different sub models"""
    # Test that the ensemble of two identical models
    # behalves like twice a single model.
    # Only for the independent (ie. diagonal) case
    # Note: this test would fail if a fixed contact_matrix
    # is passed which happens to be the identity matrix
    # The difference is that when specified as independent
    # each model is booted independently
    # If a contact matrix is specified, then the boot goal
    # is the combined total of all models. So there will be a different starting point.
    test_a = Model.open_file(path_model_2)
    test_a.name = 'test_a'
    test_b = Model.open_file(path_model_2)
    test_b.name = 'test_b'
    test_b.parameters['alpha_0'].set_value(0.7)
    test_c = Model.open_file(path_model_2)
    test_c.name = 'test_c'
    test_d = Model.open_file(path_model_2)
    test_d.name = 'test_d'
    test_d.parameters['alpha_0'].set_value(0.7)
    reference = Model.open_file(path_model_2)
    reference.name = 'reference'

    test_ensemble = Ensemble('test_ensemble', reference)
    test_ensemble.upload_models([test_a, test_b])
    test_ensemble.define_cross_transmission('infection cycle', 'infected',
                                            'susceptible', 'total',
                                            'contagious', 'alpha',
                                            contact_type='diagonal')

    n_days = 100
    test_c.reset()
    test_c.evolve_expectations(n_days)
    test_d.reset()
    test_d.evolve_expectations(n_days)
    test_ensemble.reset()
    test_ensemble.evolve_expectations(n_days)
    for pop_name in test_ensemble.populations:
        pop = test_ensemble.populations[pop_name]
        if pop.show_sim:
            ens_hist = test_ensemble.populations[pop_name].history
            tc_hist = test_c.populations[pop_name].history
            td_hist = test_d.populations[pop_name].history
            for i in range(len(ens_hist)):
                ratio = ens_hist[i] / (tc_hist[i] + td_hist[i])
                assert np.abs(ratio - 1.) < 0.01

def test_Ensemble_data():
    """tests to check data from an ensemble"""
    # Test that the ensemble of two identical models
    # behalves like twice a single model.
    # Only for the independent (ie. diagonal) case
    # Note: this test would fail if a fixed contact_matrix
    # is passed which happens to be the identity matrix
    # The difference is that when specified as independent
    # each model is booted independently
    # If a contact matrix is specified, then the boot goal
    # is the combined total of all models. So there will be a different starting point.
    test_a = Model.open_file(path_model_2)
    test_a.name = 'test_a'
    test_b = Model.open_file(path_model_2)
    test_b.name = 'test_b'
    test_b.parameters['alpha_0'].set_value(0.7)

    reference = Model.open_file(path_model_2)
    reference.name = 'reference'

    off_diag = Parameter('off_diagonal', 0.1,
                         parameter_min=0., parameter_max=1., description='off diagonal element of contact matrix')
    off_diags = [off_diag]

    test_ensemble = Ensemble('test_ensemble', reference)
    test_ensemble.upload_models([test_a, test_b])
    test_ensemble.define_cross_transmission('infection cycle', 'infected',
                                            'susceptible', 'total',
                                            'contagious', 'alpha',
                                            contact_type='simple',contact=off_diags)

    n_days = 100
    norm_day = 50
    test_ensemble.reset()
    test_ensemble.evolve_expectations(norm_day)
    for key in test_ensemble.populations:
        pop = test_ensemble.populations[key]
        nu = pop.history[norm_day]
        pop.history[norm_day] = int(round(nu))
        pop.scale_future(1., expectations=False)
    for model_name in test_ensemble.models:
        model = test_ensemble.models[model_name]
        for pop_name in model.populations:
            pop = model.populations[pop_name]
            nu = pop.history[norm_day]
            pop.history[norm_day] = int(round(nu))
            pop.scale_future(1., expectations=False)

    test_ensemble.generate_data(n_days,norm_day)
    for pop_name in test_ensemble.populations:
        pop = test_ensemble.populations[pop_name]
        if pop.show_sim:
            ens_hist = test_ensemble.populations[pop_name].history

def test_point_estimate():
    start_day = 12
    end_day = 60
    ref_2 = Model.open_file(path_model_2_2)
    sim_2 = Model.open_file(path_model_2_2)

    # do fit of alpha_0, alpha_1, cont_0, trans_rate_1_time
    for par_name in ['alpha_0', 'alpha_1', 'cont_0']:
        par = ref_2.parameters[par_name]
        par.set_variable(None, None)

    par = ref_2.parameters['trans_rate_1_time']
    par.set_variable(None, None)
    par.set_min(13)
    par.set_max(19)

    sim_2.reset()
    sim_2.generate_data(end_day)
    sim_2.populations['reported'].history[47] = np.inf
    optimizer = Optimizer(ref_2, 'total reported', sim_2.populations['reported'].history, [start_day, end_day],skip_data='42,45:48')
    optimizer.reset_variables()

    scan_dict = optimizer.i_fit()
    assert ref_2.parameters['trans_rate_1_time'].get_value() in [15,16,17]

    par = ref_2.parameters['trans_rate_1_time']
    par.set_fixed()

    popt, pcov = optimizer.fit()
    assert np.abs(ref_2.parameters['alpha_0'].get_value()-ref_2.parameters['alpha_0'].initial_value) < 0.06
    assert np.abs(ref_2.parameters['alpha_1'].get_value() - ref_2.parameters['alpha_1'].initial_value) < 0.02
    assert np.abs(ref_2.parameters['cont_0'].get_value() - ref_2.parameters['cont_0'].initial_value) < 20.

def test_point_estimate_local():
    start_day = 12
    end_day = 60
    ref_2 = Model.open_file(path_model_2_2)
    sim_2 = Model.open_file(path_model_2_2)

    # do fit of alpha_1, trans_rate_1_time
    for par_name in ['alpha_1']:
        par = ref_2.parameters[par_name]
        par.set_variable(None, None)

    par = ref_2.parameters['trans_rate_1_time']
    par.set_variable(None, None)
    par.set_min(13)
    par.set_max(19)

    sim_2.reset()
    sim_2.generate_data(end_day)
    sim_2.populations['reported'].history[47] = np.inf
    optimizer = Optimizer(ref_2, 'total reported', sim_2.populations['reported'].history, [start_day, end_day],cumul_reset=True,skip_data='42,45:48')
    optimizer.reset_variables()

    scan_dict = optimizer.i_fit()
    assert ref_2.parameters['trans_rate_1_time'].get_value() in [15,16,17]

    par = ref_2.parameters['trans_rate_1_time']
    par.set_fixed()

    popt, pcov = optimizer.fit()
    assert np.abs(ref_2.parameters['alpha_1'].get_value() - ref_2.parameters['alpha_1'].initial_value) < 0.02

def test_point_estimate_daily():

    def delta(cumul):
        diff = []
        for i in range(1, len(cumul)):
            diff.append(cumul[i] - cumul[i - 1])
        # first daily value is repeated since val(t0-1) is unknown
        diff.insert(0,diff[0])
        return diff

    start_day = 12
    end_day = 60
    ref_2 = Model.open_file(path_model_2_2)
    sim_2 = Model.open_file(path_model_2_2)

    # do fit of alpha_0, alpha_1, cont_0, trans_rate_1_time
    for par_name in ['alpha_0', 'alpha_1', 'cont_0']:
        par = ref_2.parameters[par_name]
        par.set_variable(None, None)

    par = ref_2.parameters['trans_rate_1_time']
    par.set_variable(None, None)
    par.set_min(13)
    par.set_max(19)

    sim_2.reset()
    sim_2.generate_data(end_day)
    daily_data = delta(sim_2.populations['reported'].history)
    daily_data[47] = np.inf
    optimizer = Optimizer(ref_2, 'daily reported', daily_data, [start_day, end_day],skip_data='42,45:48')
    optimizer.reset_variables()

    scan_dict = optimizer.i_fit()
    assert ref_2.parameters['trans_rate_1_time'].get_value() in [15,16,17]

    par = ref_2.parameters['trans_rate_1_time']
    par.set_fixed()

    popt, pcov = optimizer.fit()
    assert np.abs(ref_2.parameters['alpha_0'].get_value()-ref_2.parameters['alpha_0'].initial_value) < 0.06
    assert np.abs(ref_2.parameters['alpha_1'].get_value() - ref_2.parameters['alpha_1'].initial_value) < 0.02
    assert np.abs(ref_2.parameters['cont_0'].get_value() - ref_2.parameters['cont_0'].initial_value) < 20.

def test_point_estimates_repeated():
    start_day = 12
    end_day = 60
    ref_2 = Model.open_file(path_model_2_2)
    sim_2 = Model.open_file(path_model_2_2)

    # do fit of alpha_0, alpha_1, cont_0
    par_names = ['alpha_0', 'alpha_1', 'cont_0']
    sums = {}
    sum2s = {}
    for par_name in par_names:
        par = ref_2.parameters[par_name]
        par.set_variable(None, None)
        sums[par_name] = 0.
        sum2s[par_name] = 0.

    n_rep = 10
    fit_stat_list = []
    for i in range(n_rep):
        sim_2.reset()
        sim_2.generate_data(end_day)
        optimizer = Optimizer(ref_2, 'total reported', sim_2.populations['reported'].history, [start_day, end_day])
        optimizer.reset_variables()
        popt, pcov = optimizer.fit()
        fit_stat_list.append(optimizer.fit_statistics)
        for par_name in par_names:
            value = ref_2.parameters[par_name].get_value()
            sums[par_name] += value
            sum2s[par_name] += value**2

    ass_std = {}
    ass_std['alpha_0'] = 0.03
    ass_std['alpha_1'] = 0.01
    ass_std['cont_0'] = 10.

    means = {}
    std = {}
    for par_name in par_names:
        means[par_name] = sums[par_name]/n_rep
        std[par_name] = np.sqrt(sum2s[par_name]/n_rep - means[par_name]**2)
        assert std[par_name] < ass_std[par_name]
        truth = ref_2.parameters[par_name].initial_value
        assert np.abs((means[par_name]-truth)/std[par_name]/np.sqrt(1.*n_rep)) < 3.

    ndof = fit_stat_list[0]['ndof']
    chi2_list = [fit_stat_list[i]['chi2'] for i in range(n_rep)]
    chi2_mean = np.mean(chi2_list)
    assert np.abs(chi2_mean - ndof) < 8.
    acor_list = [fit_stat_list[i]['acor'] for i in range(n_rep)]
    acor_mean = np.mean(acor_list)
    assert np.abs(acor_mean) < 0.2

def test_sim_gof():
    start_day = 12
    end_day = 60
    ref_2 = Model.open_file(path_model_2_2)
    sim_2 = Model.open_file(path_model_2_2)

    # do fit of alpha_0, alpha_1, cont_0
    par_names = ['alpha_0', 'alpha_1', 'cont_0']
    for par_name in par_names:
        par = ref_2.parameters[par_name]
        par.set_variable(None, None)

    sim_2.reset()
    sim_2.generate_data(end_day)
    sim_2.populations['reported'].history[47] = np.inf
    optimizer = Optimizer(ref_2, 'total reported', sim_2.populations['reported'].history, [start_day, end_day],skip_data='42,45:48')
    optimizer.reset_variables()
    popt, pcov = optimizer.fit()
    fit_statistics = optimizer.fit_statistics

    optimizer.calc_chi2s = False
    optimizer.calc_chi2f = True
    n_rep = 10
    optimizer.calc_sim_gof(n_rep)

    fit_stat_list = optimizer.fit_stat_list
    ndof = fit_stat_list[0]['ndof']
    chi2_list = [fit_stat_list[i]['chi2'] for i in range(n_rep)]
    chi2_mean = np.mean(chi2_list)
    assert np.abs(chi2_mean - ndof) < 8.
    acor_list = [fit_stat_list[i]['acor'] for i in range(n_rep)]
    acor_mean = np.mean(acor_list)
    assert np.abs(acor_mean) < 0.2

def test_sim_gof_local():
    start_day = 12
    end_day = 60
    ref_2 = Model.open_file(path_model_2_2)
    sim_2 = Model.open_file(path_model_2_2)

    # do fit of alpha_0, alpha_1, cont_0
    par_names = ['alpha_1']
    for par_name in par_names:
        par = ref_2.parameters[par_name]
        par.set_variable(None, None)

    sim_2.reset()
    sim_2.generate_data(end_day)
    sim_2.populations['reported'].history[47] = np.inf
    optimizer = Optimizer(ref_2, 'total reported', sim_2.populations['reported'].history, [start_day, end_day],cumul_reset=True,skip_data='42,45:48')
    optimizer.reset_variables()
    popt, pcov = optimizer.fit()
    fit_statistics = optimizer.fit_statistics

    optimizer.calc_chi2s = False
    optimizer.calc_chi2f = True
    n_rep = 10
    optimizer.calc_sim_gof(n_rep)

    fit_stat_list = optimizer.fit_stat_list
    ndof = fit_stat_list[0]['ndof']
    chi2_list = [fit_stat_list[i]['chi2'] for i in range(n_rep)]
    chi2_mean = np.mean(chi2_list)
    assert np.abs(chi2_mean - ndof) < 8.E6
    acor_list = [fit_stat_list[i]['acor'] for i in range(n_rep)]
    acor_mean = np.mean(acor_list)
    assert np.abs(acor_mean) < 0.2

def test_report_noise():
    for report_noise_weekly in [False, True]:
        start_day = 12
        end_day = 80
        ref_2 = Model.open_file(path_model_2_2)
        ref_2.parameters['report_noise'].set_value(0.1)
        sim_2 = copy.deepcopy(ref_2)
        sim_2.populations['reported'].report_noise_weekly = report_noise_weekly

        # do fit of alpha_0, alpha_1, cont_0
        par_names = ['alpha_0', 'alpha_1', 'cont_0']
        sums = {}
        sum2s = {}
        for par_name in par_names:
            par = ref_2.parameters[par_name]
            par.set_variable(None, None)
            sums[par_name] = 0.
            sum2s[par_name] = 0.

        n_rep = 10
        fit_stat_list = []
        for i in range(n_rep):
            sim_2.reset()
            sim_2.generate_data(end_day)
            optimizer = Optimizer(ref_2, 'total reported', sim_2.populations['reported'].history, [start_day, end_day])
            optimizer.reset_variables()
            popt, pcov = optimizer.fit()
            fit_stat_list.append(optimizer.fit_statistics)
            for par_name in par_names:
                value = ref_2.parameters[par_name].get_value()
                sums[par_name] += value
                sum2s[par_name] += value**2

        ass_std = {}
        if not report_noise_weekly:
            ass_std['alpha_0'] = 0.03
            ass_std['alpha_1'] = 0.01
            ass_std['cont_0'] = 10.
        else:
            ass_std['alpha_0'] = 0.06
            ass_std['alpha_1'] = 0.01
            ass_std['cont_0'] = 10.

        means = {}
        std = {}
        for par_name in par_names:
            means[par_name] = sums[par_name]/n_rep
            std[par_name] = np.sqrt(sum2s[par_name]/n_rep - means[par_name]**2)
        for par_name in par_names:
            assert std[par_name] < ass_std[par_name]
            truth = ref_2.parameters[par_name].initial_value
            assert np.abs((means[par_name]-truth)/std[par_name]/np.sqrt(1.*n_rep)) < 3.

        ndof = fit_stat_list[0]['ndof']
        chi2_list = [fit_stat_list[i]['chi2'] for i in range(n_rep)]
        chi2_mean = np.mean(chi2_list)
        assert chi2_mean - ndof > 100.
        if not report_noise_weekly:
            acor_list = [fit_stat_list[i]['acor'] for i in range(n_rep)]
            acor_mean = np.mean(acor_list)
            assert acor_mean < -0.3

def test_report_noise_days():
    for report_noise_weekly in [False, True]:
        start_day = 12
        end_day = 80
        ref_2 = Model.open_file(path_model_2_2)
        ref_2.parameters['report_noise'].set_value(0.1)
        # BC: no reporting on Sundays
        ref_2.parameters['report_days'].set_value(63)
        sim_2 = copy.deepcopy(ref_2)
        sim_2.populations['reported'].report_noise_weekly = report_noise_weekly

        # do fit of alpha_0, alpha_1, cont_0
        par_names = ['alpha_0', 'alpha_1', 'cont_0']
        sums = {}
        sum2s = {}
        for par_name in par_names:
            par = ref_2.parameters[par_name]
            par.set_variable(None, None)
            sums[par_name] = 0.
            sum2s[par_name] = 0.

        n_rep = 10
        fit_stat_list = []
        for i in range(n_rep):
            sim_2.reset()
            sim_2.generate_data(end_day)
            optimizer = Optimizer(ref_2, 'total reported', sim_2.populations['reported'].history, [start_day, end_day])
            optimizer.reset_variables()
            popt, pcov = optimizer.fit()
            fit_stat_list.append(optimizer.fit_statistics)
            for par_name in par_names:
                value = ref_2.parameters[par_name].get_value()
                sums[par_name] += value
                sum2s[par_name] += value ** 2

        ass_std = {}
        ass_std['alpha_0'] = 0.05
        ass_std['alpha_1'] = 0.01
        ass_std['cont_0'] = 10.

        means = {}
        std = {}
        for par_name in par_names:
            means[par_name] = sums[par_name] / n_rep
            std[par_name] = np.sqrt(sum2s[par_name] / n_rep - means[par_name] ** 2)
        for par_name in par_names:
            assert std[par_name] < ass_std[par_name]
            truth = ref_2.parameters[par_name].initial_value
            assert np.abs((means[par_name] - truth) / std[par_name] / np.sqrt(1. * n_rep)) < 3.

        ndof = fit_stat_list[0]['ndof']
        chi2_list = [fit_stat_list[i]['chi2'] for i in range(n_rep)]
        chi2_mean = np.mean(chi2_list)
        assert chi2_mean/ndof > 2.5
        if not report_noise_weekly:
            acor_list = [fit_stat_list[i]['acor'] for i in range(n_rep)]
            acor_mean = np.mean(acor_list)
            assert acor_mean < -0.2

def test_trajectory():
    ref_2 = Model.open_file(path_model_2_2)
    trajectory = Trajectory(ref_2, 'contagious', 'trans_rate_1', [0.03, 0.75])
    alpha_c = trajectory.get_alpha(0.)
    delta_1 = trajectory.get_delta(0.1)
    delta_2 = trajectory.get_delta(0.2)
    assert np.abs(alpha_c - 0.152) < 0.001
    assert np.abs(trajectory.get_delta(alpha_c)) < 0.00001
    assert np.abs(delta_1 + 0.0513) < 0.0001
    assert np.abs(delta_2-0.0389) < 0.0001

def test_linear_modifier():
    ref_2_4 = Model.open_file(path_model_2_4)
    trajectory = Trajectory(ref_2_4, 'contagious', 'trans_rate_1', [0.02, 2.0])
    alpha_c = trajectory.get_alpha(0.)
    ref_2_4.parameters['alpha_1'].set_value(alpha_c)
    ref_2_4.parameters['to_icu_delay_mean'].set_value(0.5)
    ref_2_4.parameters['to_icu_delay_sigma'].set_value(0.5)
    ref_2_4.transitions['mod_icu_frac'].enabled = True
    ref_2_4.parameters['icu_frac_time'].set_value(60)
    ref_2_4.parameters['icu_frac_0'].set_value(0.1)
    # the meaning of this parameter changed to the end_value
    # ref_2_4.parameters['icu_frac_slope'].set_value(0.01)
    ref_2_4.parameters['icu_frac_slope'].set_value(0.2)
    ref_2_4.parameters['icu_frac_nstep'].set_value(10)
    ref_2_4.reset()
    ref_2_4.evolve_expectations(100)
    pop_hist = ref_2_4.populations['icu admissions'].history
    assert np.abs(pop_hist[55]-pop_hist[54]-pop_hist[50]+pop_hist[49]) < 0.01
    assert np.abs(pop_hist[75] - pop_hist[74] - pop_hist[70] + pop_hist[69]) < 0.1
    assert np.abs((pop_hist[75] - pop_hist[74])/(pop_hist[55]-pop_hist[54]) - 2.) < 0.1

def test_model_2_5():
    ref_2_5 = Model.open_file(path_model_2_5)
    ref_2_3 = Model.open_file(path_model_2_3)
    #ref_2_5.transitions['vaccination'].enabled = True
    #ref_2_5.parameters['vaccination_number'].set_value(10000.)
    ref_2_5.reset()
    ref_2_5.evolve_expectations(200)
    ref_2_3.reset()
    ref_2_3.evolve_expectations(200)

    for pop_name in ref_2_3.populations:
        pop = ref_2_3.populations[pop_name]
        if pop.show_sim:
            print(pop_name)
            hist_2_3 = ref_2_3.populations[pop_name].history
            hist_2_5 = ref_2_5.populations[pop_name].history
            for i in range(len(hist_2_3)):
                if hist_2_3[i] > 0:
                    ratio = hist_2_5[i] / hist_2_3[i]
                    assert np.abs(ratio - 1.) < 0.001

def test_max_mult():
    ref_2_5 = Model.open_file(path_model_2_5)
    ref_2_5.transitions['vaccination'].enabled = True
    ref_2_5.parameters['vaccination_number'].set_value(800000.)
    ref_2_5.parameters['vaccination_time'].set_value(5)
    ref_2_5.reset()
    ref_2_5.boot(expectations=True)
    ref_2_5.evolve_expectations(200)
    i=1

def test_interval_maker():

    hub_date = datetime.date(2020, 4, 1)
    my_IntervalMaker = IntervalMaker("USA", hub_date)
    categories = ['case','death','hospitalization']
    n_period_dict = {'case':5, 'death':5, 'hospitalization':30}
    n_rep = 10
    scale_std_alpha = 2.
    model = Model.open_file(path_model_2_6)
    if 'interval_maker' not in model.user_dict:
        model.user_dict['interval_maker'] = {}
        model.user_dict['interval_maker']['smearing parameters'] = ['non_icu_hosp_frac','recover_frac']
    model.parameters['non_icu_hosp_frac'].std_estimator = 0.002
    model.parameters['recover_frac'].set_value(0.99)
    model.parameters['recover_frac'].std_estimator = 0.01

    my_IntervalMaker.get_quantiles(categories, n_period_dict, model, n_rep=n_rep,
                                               scale_std_alpha=scale_std_alpha, back_up=21, fall_back=True, rescale=True)
    for category in categories:
        my_IntervalMaker.append_user_dict(category, model)
    i=1

def test_model_2_7():
    ref_2_7 = Model.open_file(path_model_2_7)
    ref_2_7.transitions['vaccination_1'].enabled = True
    ref_2_7.transitions['vaccination_2'].enabled = True
    ref_2_7.parameters['vacc_time_2'].set_value(80)
    ref_2_7.parameters['vacc_number_2'].set_value(-5.)
#    ref_2_7.reset()
#    ref_2_7.evolve_expectations(200)
    ref_2_7.reset()
    ref_2_7.generate_data(200)

    i=1

def test_model_2_8():
    ref_2_8 = Model.open_file(path_model_2_8)
    ref_2_8.transitions['outbreak_v'].enabled = True
    ref_2_8.parameters['outbreak_v_time'].set_value(30)
    ref_2_8.parameters['outbreak_v_number'].set_value(2.)
    ref_2_8.parameters['alpha_0'].set_value(0.35)
#    ref_2_8.reset()
#    ref_2_8.evolve_expectations(200)
    ref_2_8.reset()
    ref_2_8.generate_data(200)

    i=1

def test_model_2_9():
    ref_2_9 = Model.open_file(path_model_2_9)
    ref_2_9.parameters['alpha_0'].set_value(0.35)

    ref_2_9.transitions['outbreak_v'].enabled = True
    ref_2_9.parameters['outbreak_v_time'].set_value(30)
    ref_2_9.parameters['outbreak_v_number'].set_value(2.)
    ref_2_9.parameters['alpha_0_v'].set_value(0.35)

    ref_2_9.transitions['outbreak_w'].enabled = True
    ref_2_9.parameters['outbreak_w_time'].set_value(60)
    ref_2_9.parameters['outbreak_w_number'].set_value(4.)
    ref_2_9.parameters['alpha_0_w'].set_value(0.45)

    ref_2_9.transitions['vaccination_1'].enabled = True
    ref_2_9.parameters['vacc_time_1'].set_value(80)
    ref_2_9.parameters['vacc_number_1'].set_value(1000.)

#    ref_2_9.reset()
#    ref_2_9.evolve_expectations(200)

    ref_2_9.reset()
    ref_2_9.generate_data(200)

    i=1

def test_model_az_2_9():
    az_2_8 = Model.open_file('az_2_8_0530.pypm')
    az_2_9 = Model.open_file(path_model_2_9)

    az_2_9.copy_values_from(az_2_8)
    az_2_9.name = 'az_2_9_0530'

    #v_hesitant = 0.14
    #N0 = az_2_8.parameters['N_0'].get_value()
    #az_2_8.populations['vacc cand'].set_initial_value(int(N0 * (1. - v_hesitant)))
    #az_2_8.populations['sus vacc cand'].set_initial_value(int(N0 * (1. - v_hesitant)))

    #az_2_8.parameters['vacc_time_1'].set_value(290)
    #az_2_9.parameters['vacc_time_1'].set_value(290)

    az_2_8.reset()
    az_2_8.evolve_expectations(289)
    az_2_8.evolve_expectations(500-289, from_step=289)

    az_2_9.parameters['vac_waned_delay_mean'].set_value(10000.)
    az_2_9.parameters['nat_waned_delay_mean'].set_value(10000.)
    az_2_9.parameters['vac_waned_frac'].set_value(0.)

    #v_hesitant = 0.14
    #N0 = az_2_9.parameters['N_0'].get_value()
    #az_2_9.populations['vacc cand'].set_initial_value(int(N0 * (1. - v_hesitant)))
    #az_2_9.populations['sus vacc cand'].set_initial_value(int(N0 * (1. - v_hesitant)))

    az_2_9.reset()
    #az_2_9.evolve_expectations(500)
    az_2_9.evolve_expectations(289)
    az_2_9.evolve_expectations(500-289, from_step=289)

    reps = []
    rep_vs = []
    wanned = []

    pops = ['susceptible','immunized','usefully vaccinated','daily vaccinated','vacc cand','sus vacc cand']
    compares = {}
    for pop in pops:
        compares[pop] = []

    for day in [288,289,290,291,292,293,500]:
        reps_m = []
        rep_vs_m = []

        for model in [az_2_8, az_2_9]:
            rep = model.populations['reported'].history[day] - model.populations['reported'].history[day - 1]
            reps_m.append(rep)
            rep_v = model.populations['reported_v'].history[day]-model.populations['reported_v'].history[day-1]
            rep_vs_m.append(rep_v)

        for pop in pops:
            compares[pop].append([az_2_8.populations[pop].history[day],az_2_9.populations[pop].history[day]])

        wanned.append(az_2_9.populations['waned_immunity'].history[day])
        reps.append(reps_m)
        rep_vs.append(rep_vs_m)

    i=1

def test_mixing_evolve_data():
    ref_2_9 = Model.open_file(path_model_2_9)
    ref_2_9.reset()
    ref_2_9.generate_data(100,from_step=0,data_start=50)
    reported = ref_2_9.populations['reported'].history
    i=1