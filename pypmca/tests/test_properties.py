#!/usr/bin/env python

"""Tests for `pypmca` package."""

import pytest

from pypmca import Model, Population, Delay, Parameter, Multiplier, Propagator, \
    Splitter, Adder, Subtractor, Chain, Modifier, Injector, Ensemble
import numpy as np
import copy
from pathlib import Path

example_dir = Path('../../examples/').resolve()
path_model_2 = example_dir / 'ref_model_2.pypm'
path_model_1 = example_dir / 'ref_model_1.pypm'

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
