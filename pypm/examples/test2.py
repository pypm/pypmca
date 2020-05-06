# -*- coding: utf-8 -*-
"""
Created on Wed Apr  8 10:57:09 2020

@author: karlen
"""

import pickle
from pypm import Model

my_model=None

with open('model.pickle', 'rb') as f:
    # The protocol version used is detected automatically, so we do not
    # have to specify it.
    my_model = pickle.load(f)
    
i=1