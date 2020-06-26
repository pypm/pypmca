# -*- coding: utf-8 -*-
"""
Trajectory: Given a model: characterize the relation between the model transmission rate parameter (alpha) and the
infection trajectory growth rate (delta: fractional growth per day).

This is evaluated by running the model for different values of the growth rate parameter.
The interpolation function is done as delta(ln_a), but that is hidden from user

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
        if rate_range is None:
            raise TypeError('Error in constructing trajectory'+
                            ': rate_range cannot be None')
        if not isinstance(rate_range, list):
            raise TypeError('Error in constructing trajectory'+
                            ': rate_range must be a list')
        if len(rate_range) != 2:
            raise TypeError('Error in constructing trajectory'+
                            ': rate_range must be a list of length 2')
        if rate_range[0] <= 0. or rate_range[1] <= 0. or rate_range[1] <= rate_range[0]:
            raise ValueError('Error in constructing trajectory'+
                            ': rate_range list elements must be positive and second element larger than first')
        self.rate_range = rate_range
        self.log_alphas = None
        self.deltas = None
        self.delta_inter = None
        self.log_alpha_inter = None

        self.calc_deltas()

    def calc_deltas(self):
        # use 100 points for the interpolation (using log alpha)
        n_point = 100
        length = np.log(self.rate_range[1]) - np.log(self.rate_range[0])
        self.log_alphas = np.arange(np.log(self.rate_range[0]), np.log(self.rate_range[1]), length/n_point)

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
        for log_alpha in self.log_alphas:
            alpha = np.exp(log_alpha)
            rate_par.set_value(alpha)
            tr_model.reset()
            tr_model.evolve_expectations(n_step)
            ratio = tr_model.populations[self.contagious_pop_name].history[n_step] / tr_model.populations[self.contagious_pop_name].history[n_step-1]
            delta = ratio - 1
            self.deltas.append(delta)

        self.log_alpha_inter = interp1d(self.deltas, self.log_alphas, kind='cubic')
        self.delta_inter = interp1d(self.log_alphas, self.deltas, kind='cubic')

    def get_delta(self, alpha):
        return self.delta_inter(np.log(alpha)).item()

    def get_alpha(self, delta):
        return np.exp(self.log_alpha_inter(delta).item())