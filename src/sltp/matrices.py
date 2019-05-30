""" Formatting and printing of the output of the concept-based feature generation process,
along with some other related output necessary for subsequent steps in the pipeline """
import itertools
import logging
import math
import sys

import numpy as np

from tarski.dl import EmpiricalBinaryConcept, NullaryAtomFeature, ConceptCardinalityFeature, MinDistanceFeature


PRUNE_DUPLICATE_FEATURES = True
NP_FEAT_VALUE_TYPE = np.int8  # Keep it allowing negative values, so that we can subtract without overflow!


def next_power_of_two(x):
    """ Return the smallest power of two Z such that Z >= x """
    if x == 0:
        return 0
    return 2 ** (math.ceil(math.log2(x)))


def cast_feature_value_to_numpy_value(value):
    """ Cast a given feature value into a suitable numpy value, if possible, or raise error if not """
    assert value >= 0
    max_ = np.iinfo(NP_FEAT_VALUE_TYPE).max
    if value == sys.maxsize or value == 2147483647:  # std::numeric_limits<int>::max(). Yes, this is not portable :-)
        return max_

    if value >= max_:  # Max value reserved to denote infty.
        raise RuntimeError("Cannot cast feature value {} into numpy value".format(value))

    return value


def print_feature_info(config, features):
    filename = config.feature_info_filename
    logging.info("Printing feature info for {} features to '{}'".format(len(features), filename))

    with open(filename, 'w') as f:
        for feat in features:
            print("{} {}".format(feat, feat.complexity()), file=f)


def print_sat_feature_matrix(config, state_ids, goals, features, models):
    filename = config.sat_feature_matrix_filename
    logging.info("Printing SAT feature matrix of {} features x {} states to '{}'".
                 format(len(features), len(state_ids), filename))

    # Separate features into binary and numeric
    binary, numeric = [], []
    for i, feat in enumerate(features, 0):
        if isinstance(feat, (NullaryAtomFeature, EmpiricalBinaryConcept)):
            binary.append(i)
        elif isinstance(feat, (ConceptCardinalityFeature, MinDistanceFeature)):
            numeric.append(i)
        else:
            assert 0, "Unknown feature type"

    num_numeric = len(numeric)
    all_feature_idxs = numeric + binary
    nfeatures = len(features)
    ngoals = len(goals)

    with open(filename, 'w') as f:
        # Header row: <#states> <#features> <#goals> <1 + index-last-numerical-feature> <index-first-boolean-feature>
        print("{} {} {} {} {}".format(len(state_ids), nfeatures, ngoals, num_numeric, num_numeric), file=f)

        # Line #2:: <#features> <list of feature names>
        print("{} {}".format(nfeatures, " ".join(str(features[f]).replace(" ", "") for f in all_feature_idxs)), file=f)

        # Line #3: <#features> <list of feature costs>
        print("{} {}".format(nfeatures, " ".join(str(features[f].complexity()) for f in all_feature_idxs)), file=f)

        # Line #4: <# goal states> <list of goal state IDs>
        print("{} {}".format(ngoals, " ".join(map(str, goals))), file=f)

        # next lines: one per each state with format: <state-index> <#features-in-state> <list-features>
        # each feature has format: <feature-index>:<value>
        for s in state_ids:
            feature_values = ((i, models[s].denotation(features[f])) for i, f in enumerate(all_feature_idxs, 0))
            nonzero_features = [(i, v) for i, v in feature_values if v != 0]
            print("{} {} {}".format(s, len(nonzero_features), " ".join(
                "{}:{}".format(i, int(v)) for i, v in nonzero_features
            )), file=f)

    return all_feature_idxs  # A mapping between the new and the old feature indexes


def log_feature_denotations(state_ids, features, models, feature_denotation_filename, selected=None):
    selected = selected or features
    selected = ((str(f), f) for f in selected)
    selected = sorted(selected, key=lambda x: x[0])  # Sort features by name

    with open(feature_denotation_filename, 'w') as file:
        for s, (fname, f) in itertools.product(state_ids, selected):
            val = models[s].denotation(f)
            print("s_{}[{}] = {}".format(s, fname, val), file=file)
    logging.info("Logged feature denotations at '{}'".format(feature_denotation_filename))


def printer(feature, value):
    return "1" if feature.bool_value(value) else "0"


def int_printer(value):
    return str(int(value))
