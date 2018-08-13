""" Formatting and printing of the output of the concept-based feature generation process,
along with some other related output necessary for subsequent steps in the pipeline """
import itertools
import logging
import math
import numpy as np

from tarski.dl import EmpiricalBinaryConcept, NullaryAtomFeature, ConceptCardinalityFeature, MinDistanceFeature

from features import Model

PRUNE_DUPLICATE_FEATURES = True


def next_power_of_two(x):
    """ Return the smallest power of two Z such that Z >= x """
    if x == 0:
        return 0
    return 2 ** (math.ceil(math.log2(x)))


def print_transition_matrix(state_ids, transitions, transitions_filename):
    num_transitions = sum(len(targets) for targets in transitions.values())
    logging.info("Printing transition matrix with {} states and {} transitions to '{}'".
                 format(len(state_ids), num_transitions, transitions_filename))
    with open(transitions_filename, 'w') as f:
        for s, succ in transitions.items():
            print("{} {}".format(s, " ".join("{}".format(sprime) for sprime in succ)), file=f)


def print_sat_transition_matrix(state_ids, transitions, transitions_filename):
    num_transitions = sum(len(targets) for targets in transitions.values())
    num_states = len(state_ids)
    num_expanded_states = len(transitions.keys())
    logging.info("Printing SAT transition matrix with {} states, {} expanded states and {} transitions to '{}'".
                 format(len(state_ids), num_expanded_states, num_transitions, transitions_filename))
    with open(transitions_filename, 'w') as f:
        # first line: <#states> <#transitions>
        print("{} {}".format(num_states, num_transitions), file=f)

        # second line: <#expanded states>
        print("{}".format(num_expanded_states), file=f)

        for s, succ in transitions.items():
            print("{} {} {}".format(s, len(succ), " ".join("{}".format(sprime) for sprime in sorted(succ))), file=f)


def print_feature_matrix(config, state_ids, features, model):
    filename = config.feature_matrix_filename
    logging.info("Printing feature matrix of {} features x {} states to '{}'".
                 format(len(features), len(state_ids), filename))

    with open(filename, 'w') as f:
        for s in state_ids:
            line = "{} {}".format(s, " ".join(int_printer(model.compute_feature_value(feat, s)) for feat in features))
            print(line, file=f)

    # Convert features to a numpy array with n rows and m columns, where n=num states, m=num features
    matrix = np.empty(shape=(len(state_ids), len(features)), dtype=np.int8)
    bin_matrix = np.empty(shape=(len(state_ids), len(features)), dtype=np.bool_)
    for i, s in enumerate(state_ids, 0):
        assert i == s
        matrix[i] = [model.compute_feature_value(feat, s) for feat in features]
        bin_matrix[i] = [feat.bool_value(model.compute_feature_value(feat, s)) for feat in features]

    np.save(config.feature_matrix_filename, matrix)
    np.save(config.bin_feature_matrix_filename, bin_matrix)

    np.save(config.feature_complexity_filename,
            np.array([f.complexity() for f in features], dtype=np.uint8))

    np.save(config.feature_names_filename,
            np.array([str(f) for f in features], dtype=np.unicode_))

    # np.savez(config.feature_matrix_filename, matrix)


def print_feature_info(config, features):
    filename = config.feature_info_filename
    logging.info("Printing feature info for {} features to '{}'".format(len(features), filename))

    with open(filename, 'w') as f:
        for feat in features:
            print("{} {}".format(feat, feat.complexity()), file=f)


def print_sat_matrix(config, state_ids, features, model):
    filename = config.sat_feature_matrix_filename
    logging.info("Printing SAT feature matrix of {} features x {} states to '{}'".
                 format(len(features), len(state_ids), filename))

    # Separate features into binary and numeric
    binary, numeric = [], []
    for feat in features:
        if isinstance(feat, (NullaryAtomFeature, EmpiricalBinaryConcept)):
            binary.append(feat)
        elif isinstance(feat, (ConceptCardinalityFeature, MinDistanceFeature)):
            numeric.append(feat)
        else:
            assert 0, "Unknown feature type"

    num_numeric = len(numeric)
    npt = next_power_of_two(num_numeric)
    num_dummy = npt-num_numeric
    all_with_dummy = numeric + ["dummy()"]*num_dummy + binary

    with open(filename, 'w') as f:
        # Header row: <#states> <#features> <1 + index-last-numerical-feature> <index-first-boolean-feature>
        print("{} {} {} {}".format(len(state_ids), len(all_with_dummy), num_numeric, npt), file=f)

        # second line: <#features> <list of feature names>
        print("{} {}".format(len(all_with_dummy), " ".join(str(f).replace(" ", "") for f in all_with_dummy)), file=f)

        # next lines: one per each state with format: <state-index> <#features-in-state> <list-features>
        # each feature has format: <feature-index>:<value>
        for s in state_ids:
            feature_values = ((i, model.compute_feature_value(f, s)) for i, f in enumerate(all_with_dummy, 0)
                              if not isinstance(f, str))
            nonzero_features = [(i, v) for i, v in feature_values if v != 0]
            print("{} {} {}".format(s, len(nonzero_features), " ".join(
                "{}:{}".format(i, int(v)) for i, v in nonzero_features
            )), file=f)


def print_goal_states(goal_states, filename):
    with open(filename, 'w') as f:
        print(" ".join(str(s) for s in goal_states), file=f)


def generate_features(config, data):
    state_ids, features, model = compute_features(config, data)
    print_feature_info(config, features)
    print_feature_matrix(config, state_ids, features, model)
    print_transition_matrix(state_ids, data.transitions, config.transitions_filename)
    print_sat_matrix(config, state_ids, features, model)
    print_sat_transition_matrix(state_ids, data.transitions, config.sat_transitions_filename)
    print_goal_states(data.goal_states, config.goal_states_filename)
    log_features(features, config.feature_filename)
    return dict(state_ids=state_ids, features=features, feature_model=model)


def compute_features(config, data):
    state_ids = sorted(list(data.states.keys()))

    # First keep only those features which are able to distinguish at least some pair of states
    features, model = compute_feature_extensions(state_ids, data.features, data.extensions)

    # DEBUGGING: FORCE A CERTAIN SET OF FEATURES:
    # selected = [translator.features[i] for i in [8, 9, 25, 26, 1232]]
    # selected = [translator.features[i] for i in [8, 9, 25, 26, 2390]]
    # selected = [translator.features[i] for i in [1232]]
    # translator.features = selected
    selected = None
    log_feature_denotations(state_ids, features, model, config.feature_denotation_filename, selected)

    return state_ids, features, model


def log_features(features, feature_filename):
    logging.info("Feature set saved in file {}".format(feature_filename))
    with open(feature_filename, 'w') as f:
        for i, feat in enumerate(features, 0):
            f.write("{}:\t{}\n".format(i, feat))
            # serialized = serialize_to_string(feat)
            # f.write("{}:\t{}\n".format(feat, serialized))
            # assert deserialize_from_string(serialized) == feat


def compute_feature_extensions(states, features, concept_extensions):
    """ Cache all feature denotations and prune those which have constant denotation at the same time """
    model = Model(concept_extensions)
    ns, nf = len(states), len(features)
    logging.info("Computing feature denotations over a total of "
                 "{} states and {} features ({:0.1f}K matrix entries)".format(ns, nf, nf*ns/1000))
    accepted = []
    traces = dict()
    for f in features:
        all_equal, all_0_or_1 = True, True
        previous = None
        all_denotations = []
        for s in states:
            denotation = model.compute_feature_value(f, s)
            if previous is not None and previous != denotation:
                all_equal = False
            if denotation not in (0, 1):
                all_0_or_1 = False
            previous = denotation
            all_denotations.append(denotation)

        # If the denotation of the feature is exactly the same in all states, we remove it
        if all_equal:
            logging.debug("Feature \"{}\" has constant denotation ({}) over all states and will be ignored"
                          .format(f, previous))
            continue

        if PRUNE_DUPLICATE_FEATURES:
            trace = tuple(all_denotations)
            duplicated = traces.get(trace, None)
            if duplicated is None:
                traces[trace] = f
            else:
                logging.debug('Feature "{}" has same denotation trace as feature "{}" and will be ignored'
                              .format(f, duplicated))
                # logging.debug(" ".join(map(str, trace)))
                continue

        if all_0_or_1 and not isinstance(f, (NullaryAtomFeature, MinDistanceFeature)):
            accepted.append(EmpiricalBinaryConcept(f))
        else:
            accepted.append(f)

    logging.info("{}/{} features have constant or duplicated denotations and have been pruned"
                 .format(len(features)-len(accepted), len(features)))
    return accepted, model


def log_feature_denotations(state_ids, features, model, feature_denotation_filename, selected=None):
    selected = selected or features
    selected = ((str(f), f) for f in selected)
    selected = sorted(selected, key=lambda x: x[0])  # Sort features by name

    with open(feature_denotation_filename, 'w') as file:
        for s, (fname, f) in itertools.product(state_ids, selected):
            val = model.compute_feature_value(f, s)
            print("s_{}[{}] = {}".format(s, fname, val), file=file)


def printer(feature, value):
    return "1" if feature.bool_value(value) else "0"


def int_printer(value):
    return str(int(value))

