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
import itertools
import logging
import os
from collections import defaultdict
from enum import Enum
from signal import signal, SIGPIPE, SIG_DFL
import time

import numpy as np

from errors import CriticalPipelineError
from tarski.dl import FeatureValueChange, MinDistanceFeature
from util.console import print_header, print_lines
from util.command import count_file_lines, remove_duplicate_lines, read_file
from solvers import solve
from util.performance import print_memory_usage
from util.serialization import serialize_to_string


def generate_encoding(config, data, rng):
    config.qnp_abstraction_filename
    config.qnp_prefix


    states, actions, features = data.cnf_translator.compute_action_model(data.cnf_solution.assignment, data.features, config)
    data.cnf_translator.compute_qnp(states, actions, features, config, data)
    # return dict(abstract_states=states, abstract_actions=actions)
    return dict()
