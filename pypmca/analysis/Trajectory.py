# -*- coding: utf-8 -*-
"""
Trajectory: Given a model: characterize the relation between the model transmission rate parameter (alpha) and the
infection trajectory growth rate (delta: fractional growth per day).

This is evaluated by running the model for different values of the growth rate parameter.

After doing so, the following methods can be called:
- get_delta(alpha)
- get_alpha(delta)

It is assumed that the model has a consistent relation between alpha and delta throughout (time delays not changing)

@author: karlen
"""

import copy
from scipy.interpolate import interp1d
import numpy as np

class Trajectory:
    """ Trajectory: relationship between model transmission rate parameter and infection trajectory growth parameter
    """

    def __init__(self, model, contagious_pop_name: str, rate_modifier_name: str, rate_range: list):
        self.model = model
        self.contagious_pop_name = contagious_pop_name
        self.rate_modifier_name = rate_modifier_name
        self.rate_range = rate_range
        self.alphas = None
        self.deltas = None
        self.delta_inter = None
        self.alpha_inter = None

        self.calc_deltas()

    def calc_deltas(self):
        # use 100 points for the interpolation
        n_point = 100
        length = self.rate_range[1] - self.rate_range[0]
        self.alphas = np.arange(self.rate_range[0], self.rate_range[1], length/n_point)

        # make a copy of model, since it needs to be modified
        tr_model = copy.deepcopy(self.model)
        rate_modifier = tr_model.transitions[self.rate_modifier_name]
        rate_step = rate_modifier.trigger_step
        # turn off all transitions on or after the rate_modifier
        for trans_name in tr_model.transitions:
            if trans_name != self.rate_modifier_name:
                transition = tr_model.transitions[trans_name]
                if transition.trigger_step >= rate_step:
                    transition.enabled = False

        rate_par = rate_modifier.parameter_after

        # run the simulation 20 steps after the transition - use the last points to work out growth rate
        n_step = rate_step + 20

        self.deltas = []
        for alpha in self.alphas:
            rate_par.set_value(alpha)
            tr_model.reset()
            tr_model.evolve_expectations(n_step)
            ratio = tr_model.populations[self.contagious_pop_name].history[n_step] / tr_model.populations[self.contagious_pop_name].history[n_step-1]
            delta = ratio - 1
            self.deltas.append(delta)

        self.alpha_inter = interp1d(self.deltas, self.alphas, kind='cubic')
        self.delta_inter = interp1d(self.alphas, self.deltas, kind='cubic')

    def get_delta(self, alpha):
        return self.delta_inter(alpha)

    def get_alpha(self, delta):
        return self.alpha_inter(delta)