#!/usr/bin/env python

"""Tests for `pypm` package."""

import pytest

from pypm import Model, Population, Delay, Parameter, Multiplier, Propagator, \
    Splitter, Adder, Subtractor, Chain, Modifier, Injector
import numpy as np

def test_Model_properties():
    """tests to ensure the properties of Model"""
    test_model44 = Model.open_file('../examples/model_v4_4.pypm')
    test_model42 = Model.open_file('../examples/model_v4_2.pypm')

    for test_model in [test_model44, test_model42]:
        # check simple scaling of initial contagious population
        # Will not be exact due to bootstrap
        EPS = 0.01
        cont_0 = test_model.parameters['cont_0'].get_value()
        test_model.reset()
        test_model.evolve_expectations(200)
        n_reported1 = test_model.populations['reported'].history[-1]
        assert test_model.populations['contagious'].history[-1] < 10.

        test_model.parameters['cont_0'].set_value(2*cont_0)
        test_model.reset()
        test_model.evolve_expectations(200)
        n_reported2 = test_model.populations['reported'].history[-1]
        assert np.abs(n_reported2 - n_reported1*2) < EPS*n_reported1*2

    # check that the mean of many data runs is near expectation


def test_Population_properties():
    """tests to ensure the properties of Population"""


