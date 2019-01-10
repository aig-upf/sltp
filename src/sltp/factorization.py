#!/usr/bin/env python

#  Copyright (C) 2018-<date> Blai Bonet
#
#  Permission is hereby granted to distribute this software for
#  non-commercial research purposes, provided that this copyright
#  notice is included with any such distribution.
#
#  THIS SOFTWARE IS PROVIDED "AS IS" WITHOUT WARRANTY OF ANY KIND,
#  EITHER EXPRESSED OR IMPLIED, INCLUDING, BUT NOT LIMITED TO, THE
#  IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
#  PURPOSE.  THE ENTIRE RISK AS TO THE QUALITY AND PERFORMANCE OF THE
#  SOFTWARE IS WITH YOU.  SHOULD THE PROGRAM PROVE DEFECTIVE, YOU
#  ASSUME THE COST OF ALL NECESSARY SERVICING, REPAIR OR CORRECTION.
#
#  Blai Bonet, bonet@ldc.usb.ve, bonetblai@gmail.com
#  Guillem Frances, guillem.frances@unibas.ch
import logging
import os

from .errors import CriticalPipelineError


BASEDIR = os.path.dirname(os.path.realpath(__file__))


def compute_factorization(transitions, num_features, max_effects_per_action):
    # Create the variables
    pass

    # Create the constraints
    pass




def learn_factorization(config, data, rng):
    optimization = config.optimization if hasattr(config, "optimization") else OptimizationPolicy.NUM_FEATURES

    state_ids = data.state_ids
    transitions = data.transitions
    max_effects_per_action = config.factorization_max_effects_per_action  # M
    num_features = config.factorization_num_features

    logging.info("Searching for state factorization. Num. features: {}, max. effects per action: {}"
                 .format(num_features, max_effects_per_action))

    return dict(factorization=compute_factorization(transitions, num_features, max_effects_per_action))


