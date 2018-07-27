""" Formatting and printing of the output of the concept-based feature generation process,
along with some other related output necessary for subsequent steps in the pipeline """
import logging
from collections import defaultdict

from tarski.dl import EmpiricalBinaryConcept, NullaryAtomFeature

from features import Model

PRUNE_DUPLICATE_FEATURES = True


def print_transition_matrix(state_ids, transitions, transitions_filename):
    logging.info("Logging transition matrix with {} states to '{}'".
                 format(len(state_ids), transitions_filename))
    with open(transitions_filename, 'w') as f:
        for s, succ in transitions.items():
            for sprime in succ:
                print("{} {}".format(s, sprime), file=f)


def print_feature_matrix(config, state_ids, features, model):
    logging.info("Printing feature matrix of {} features x {} states".format(len(features), len(state_ids)))
    log_features(features, config.feature_filename)
    log_feature_matrix(features, state_ids, model, config.feature_matrix_filename)


def generate_features(config, data):
    state_ids, features, model = compute_features(config, data)
    print_feature_matrix(config, state_ids, features, model)
    print_transition_matrix(state_ids, data.transitions, config.transitions_filename)
    return dict(state_ids=state_ids, features=features, feature_model=model)


def compute_features(config, data):
    state_ids = sorted(list(data.states.keys()))

    # First keep only those features which are able to distinguish at least some pair of states
    features, model = compute_feature_extensions(state_ids, data.features, data.extensions)

    # DEBUGGING: FORCE A CERTAIN SET OF FEATURES:
    # selected = None
    # selected = [translator.features[i] for i in [8, 9, 25, 26, 1232]]
    # selected = [translator.features[i] for i in [8, 9, 25, 26, 2390]]
    # selected = [translator.features[i] for i in [1232]]
    # translator.features = selected
    # translator.log_feature_denotations(config.feature_denotation_filename, selected)

    return state_ids, features, model


def log_features(features, feature_filename):
    logging.info("Feature set saved in file {}".format(feature_filename))
    with open(feature_filename, 'w') as f:
        for i, feat in enumerate(features, 0):
            f.write("{}:\t{}\n".format(i, feat))
            # serialized = serialize_to_string(feat)
            # f.write("{}:\t{}\n".format(feat, serialized))
            # assert deserialize_from_string(serialized) == feat


def log_feature_matrix(features, state_ids, model, feature_matrix_filename):
    logging.info("Logging feature matrix with {} features and {} states to '{}'".
                 format(len(features), len(state_ids), feature_matrix_filename))

    def feature_printer(feature, value):
        return "1" if feature.bool_value(value) else "0"

    def int_feature_printer(value):
        return str(int(value))

    with open(feature_matrix_filename + ".int", 'w') as f_int:
        with open(feature_matrix_filename, 'w') as f:
            print(" ".join(map(str, state_ids)), file=f)   # Header row with state IDs
            for feat in features:
                print(" ".join(feature_printer(feat, model.compute_feature_value(feat, s)) for s in state_ids), file=f)
                print(" ".join(int_feature_printer(model.compute_feature_value(feat, s)) for s in state_ids), file=f_int)


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
                # logging.debug('Feature "{}" has same denotation trace as feature "{}" and will be ignored'
                #               .format(f, duplicated))
                # print(" ".join(trace))
                continue

        if all_0_or_1 and not isinstance(f, NullaryAtomFeature):
            accepted.append(EmpiricalBinaryConcept(f))
        else:
            accepted.append(f)

    logging.info("{}/{} features have constant or duplicated denotations and have been pruned"
                 .format(len(features)-len(accepted), len(features)))
    return accepted, model


def log_feature_denotations(self, feature_denotation_filename, selected=None):
    selected = selected or self.features
    d = defaultdict(list)
    for (f, s), v in self.model.cache.feature_values.items():
        if f in selected:
            d[s].append((f, v))

    states = sorted(d.keys())
    with open(feature_denotation_filename, 'w') as file:
        for s in states:
            for (f, v) in sorted(d[s], key=lambda x: str(x[0])):
                print("{} on state {}: {}".format(f, s, v), file=file)
