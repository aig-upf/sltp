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
import logging
import itertools

import os
import stat
from collections import defaultdict

from tarski.dl import PrimitiveConcept, UniversalConcept, NullaryAtom, NominalConcept, GoalConcept, GoalRole, \
    EmptyConcept, GoalNullaryAtom

from . import BASE_DIR
from .matrices import NP_FEAT_VALUE_TYPE, cast_feature_value_to_numpy_value, log_feature_denotations
from .models import DLModel
from .features import parse_all_instances, compute_models, InstanceInformation
from .util.command import execute
from .returncodes import ExitCode


def run(config, data, rng):
    return extract_features(config, data.sample)


def convert_tuple(universe, name, values):
    return (name, ) + tuple(universe.value(i) for i in (values if isinstance(values, tuple) else (values, )))


def goal_predicate_name(name):
    return "{}_g".format(name)


def process_predicate_atoms(universe, p, extension, atoms):
    if (isinstance(p, PrimitiveConcept) and p.name == 'object') or \
            isinstance(p, (UniversalConcept, EmptyConcept)) or isinstance(p, NominalConcept):
        return

    name = goal_predicate_name(p.name) \
        if isinstance(p, (GoalConcept, GoalRole, GoalNullaryAtom)) else p.name  # HACK HACK HACK

    if isinstance(p, NullaryAtom):
        if extension is True:
            atoms.append((name, ))
    else:
        for point in extension:
            atoms.append(convert_tuple(universe, name, point))


def serialize_dl_model(model: DLModel, info: InstanceInformation):
    serialized = []
    universe = info.universe

    # Add fluent information
    for k, v in model.primitive_denotations.items():
        process_predicate_atoms(universe, k, v, serialized)

    # Add static information
    for k, v in model.statics.items():
        process_predicate_atoms(universe, k, v, serialized)

    return serialized


def serialize_static_info(model: DLModel, info: InstanceInformation):
    serialized = []
    universe = info.universe

    # Add fluent information
    for k, v in model.primitive_denotations.items():
        for point in v:
            serialized.append(convert_tuple(universe, k.name, point))

    # Add static information
    for k, v in model.statics.items():
        if (isinstance(k, PrimitiveConcept) and k.name == 'object') or isinstance(k, UniversalConcept):
            continue

        for point in v:
            serialized.append(convert_tuple(universe, k.name, point))

    return serialized


def print_sample_info(sample, infos, model_cache, all_predicates, all_functions, nominals,
                      all_objects, goal_predicate_info, workspace):

    sample_fn = os.path.join(workspace, "sample.io")
    state_info = []
    atoms_per_instance = defaultdict(set)
    # Iterate over all states and collect the necessary information
    for expected_id, (id_, state) in enumerate(sample.states.items(), 0):
        # Check all ids are consecutive, as expected
        assert expected_id == id_
        instance_id = sample.instance[id_]  # The instance to which the state belongs
        full_state = serialize_dl_model(model_cache.models[id_], infos[instance_id])
        atoms_per_instance[instance_id].update(full_state)
        state_info.append([str(instance_id)] + [",".join(atom) for atom in full_state])

    logging.info("Printing sample information to {}".format(sample_fn))
    with open(sample_fn, "w") as f:
        # First line: sample name
        print("dummy-sample-name", file=f)

        # Second line: all predicate names
        print(" ".join("{}/{}".format(name, arity) for name, arity in sorted(all_predicates)), file=f)

        # Third line: all function names
        print(" ".join("{}/{}".format(name, arity) for name, arity in sorted(all_functions)), file=f)

        # Next: List of all predicates and functions mentioned in the goal
        goal_predicate_info = {}
        print(" ".join("{}".format(name) for name, arity in sorted(goal_predicate_info)), file=f)

        # Next: per-instance information.
        # Number of instances:
        print(len(infos), file=f)

        assert len(infos) == len(all_objects) == len(atoms_per_instance)
        for instance_id, objects in enumerate(all_objects, start=0):
            # all object names in instance i
            print(" ".join(sorted(objects)), file=f)

            # all possible atoms in instance i
            print("\t".join(",".join(atom) for atom in sorted(atoms_per_instance[instance_id])), file=f)

        # Next: all states. One state per line. first column is instance_id of state, rest are all atoms in that state,
        # including static atoms and type-predicate atoms
        for stinfo in state_info:
            # print one line per state with all state atoms, e.g. at,bob,shed   at,spanner,location;
            print("\t".join(stinfo), file=f)

    nominals_fn = os.path.join(workspace, "nominals.io")
    if nominals:
        logging.info("Printing information on nominal concepts to {}".format(nominals_fn))
        with open(nominals_fn, "w") as f:
            # Print off the desired nominals
            print(" ".join("{}".format(name) for name in sorted(x.symbol for x in nominals)), file=f)
    else:
        open(nominals_fn, 'w').close()  # Just write an empty file


def transform_generator_output(config, sample, matrix_filename, info_filename):
    logging.info("Transforming generator output to numpy format...")
    import numpy as np

    logging.info("Reading feature information from {}".format(info_filename))
    with open(info_filename, "r") as f:
        # Line  # 1: feature names
        names = f.readline().rstrip().split("\t")

        # Line #2: feature complexities
        complexities = [int(x) for x in f.readline().rstrip().split("\t")]

        # Line #3: feature types (0: boolean; 1: numeric)
        types = [int(x) for x in f.readline().rstrip().split("\t")]

        # Line #4: whether feature is goal feature (0: No; 1: yes)
        in_goal_features = [bool(int(x)) for x in f.readline().rstrip().split("\t")]
        in_goal_features = set(i for i, v in enumerate(in_goal_features, 0) if v)  # Transform into a set of indexes

    num_features = len(names)
    assert num_features == len(complexities) == len(types)
    num_numeric = sum(types)
    num_binary = num_features - num_numeric
    logging.info("Read {} features: {} numeric + {} binary".format(num_features, num_numeric, num_binary))

    logging.info("Reading denotation matrix from {}".format(matrix_filename))
    # Convert features to a numpy array with n rows and m columns, where n=num states, m=num features

    with open(matrix_filename, "r") as f:
        num_states = sum(1 for _ in f)

    matrix = np.empty(shape=(num_states, num_features), dtype=NP_FEAT_VALUE_TYPE)

    with open(matrix_filename, "r") as f:
        # One line per state with the numeric denotation of all features
        for i, line in enumerate(f, start=0):
            data = line.rstrip().split(' ')
            matrix[i] = [cast_feature_value_to_numpy_value(int(x)) for x in data]

    logging.info("Read denotation matrix with dimensions {}".format(matrix.shape))
    sat_feature_mapping, in_goal_features = print_actual_output(sample, config, in_goal_features, names, complexities, matrix, types)
    return sat_feature_mapping, in_goal_features, len(names)


def generate_output_from_handcrafted_features(sample, config, features, model_cache):
    import numpy as np
    num_features = len(features)
    num_states = len(sample.states)

    names = [str(f) for f in features]
    types = [1]*num_features  # Let us harmlessly assume that all features are numeric
    complexities = [f.complexity() for f in features]
    matrix = np.empty(shape=(num_states, num_features), dtype=NP_FEAT_VALUE_TYPE)

    for i, (sid, atoms) in enumerate(sample.states.items(), start=0):
        assert i == sid
        model = model_cache.get_feature_model(sid)
        denotations = [model.denotation(f) for f in features]
        matrix[i] = [cast_feature_value_to_numpy_value(int(x)) for x in denotations]

    # These next 3 lines just to print the denotation of all features
    # state_ids = sample.get_sorted_state_ids()
    # models = {sid: model_cache.get_feature_model(sid) for sid in sample.states}
    # log_feature_denotations(state_ids, features, models, config.feature_denotation_filename, None)

    sat_feature_mapping, in_goal_features = print_actual_output(sample, config, [], names, complexities, matrix, types)
    return sat_feature_mapping, in_goal_features, len(names)


def print_actual_output(sample, config, in_goal_features, names, complexities, matrix, types):
    import numpy as np
    bin_matrix = np.array(matrix, dtype=np.bool)
    np.save(config.feature_matrix_filename, matrix)
    np.save(config.bin_feature_matrix_filename, bin_matrix)
    np.save(config.feature_complexity_filename, np.array(complexities, dtype=NP_FEAT_VALUE_TYPE))
    np.save(config.feature_names_filename, np.array(names, dtype=np.unicode_))

    state_ids = sample.get_sorted_state_ids()
    sat_feature_mapping = print_blai_sat_feature_matrix(config.sat_feature_matrix_filename,
                                                        matrix, state_ids, sample.goals,
                                                        names, complexities, types)

    print_feature_matrix(config.feature_matrix_filename, matrix, state_ids, sample.goals, names, complexities)
    return sat_feature_mapping, in_goal_features


def generate_debug_scripts(target_dir, exe, arguments):
    # If generating a debug build, create some debug script helpers
    shebang = "#!/usr/bin/env bash"
    args = ' '.join(arguments)
    debug_script = "{}\n\n cgdb -ex=run --args {} {}".format(shebang, exe, args)
    memleaks = "{}\n\n valgrind --leak-check=full --show-leak-kinds=all --num-callers=50 --track-origins=yes " \
               "--log-file=\"valgrind-output.$(date '+%H%M%S').txt\" {} {}"\
        .format(shebang, exe, args)

    memprofile = "{}\n\n valgrind --tool=massif {} {}".format(shebang, exe, args)

    make_script(os.path.join(target_dir, 'debug.sh'), debug_script)
    make_script(os.path.join(target_dir, 'memleaks.sh'), memleaks)
    make_script(os.path.join(target_dir, 'memprofile.sh'), memprofile)


def make_script(filename, code):
    with open(filename, 'w') as f:
        print(code, file=f)
    st = os.stat(filename)
    os.chmod(filename, st.st_mode | stat.S_IEXEC)


def extract_features(config, sample):
    logging.info("Generating non-redundant concepts from sample set: {}".format(sample.info()))

    parsed_problems = parse_all_instances(config.domain, config.instances)  # Parse all problem instances

    language, nominals, model_cache, infos = compute_models(
        config.domain, sample, parsed_problems, config.parameter_generator)

    # If user provides handcrafted features, no need to go further than here
    if config.feature_generator is not None:
        logging.info('Skipping automatic feature generation: User provided set of handcrafted features')
        features = config.feature_generator(language)
        generate_output_from_handcrafted_features(sample, config, features, model_cache)
        return ExitCode.Success, dict(enforced_feature_idxs=[], in_goal_features=[], sat_feature_mapping={},
                                      model_cache=model_cache,)

    all_goal_predicates = set(itertools.chain.from_iterable(info.goal_predicates for info in infos))
    if any(all_goal_predicates != info.goal_predicates for info in infos):
        logging.warning("Not all instances in the training set use the same goal predicate")

    all_objects = []
    all_predicates, all_functions = set(), set()
    goal_predicate_info = set()
    for problem, lang, _ in parsed_problems:
        all_objects.append(set(c.symbol for c in lang.constants()))
        all_predicates.update((p.name, p.arity) for p in lang.predicates if not p.builtin)
        all_functions.update((p.name, p.arity) for p in lang.functions if not p.builtin)

        # Add goal predicates
        goal_predicate_info = set((goal_predicate_name(p.name), p.arity)
                                  for p in lang.predicates if not p.builtin and p.name in all_goal_predicates)
        all_predicates.update(goal_predicate_info)

        # Add type predicates
        all_predicates.update((p.name, 1) for p in lang.sorts if not p.builtin and p != lang.Object)

    # Write sample information
    print_sample_info(sample, infos, model_cache, all_predicates, all_functions, nominals,
                      all_objects, goal_predicate_info, config.experiment_dir)

    # Invoke C++ feature generation module
    logging.info('Invoking C++ feature generation module'.format())
    featuregen_location = os.path.join(BASE_DIR, "..", "features")
    cmd = os.path.realpath(os.path.join(featuregen_location, "featuregen"))
    args = [str(config.max_concept_size), config.experiment_dir]
    generate_debug_scripts(config.experiment_dir, cmd, args)
    retcode = execute([cmd] + args)

    if retcode != 0:
        return ExitCode.FeatureGenerationUnknownError, dict()

    # Read off the output of the module and transform it into the numpy matrices to be consumed
    # by the next pipeline step
    sat_feature_mapping, in_goal_features, nfeatures = transform_generator_output(
        config, sample,
        os.path.join(config.experiment_dir, "feature-matrix.io"),
        os.path.join(config.experiment_dir, "feature-info.io"),)

    return ExitCode.Success, dict(enforced_feature_idxs=[],
                                  in_goal_features=in_goal_features,
                                  num_features=nfeatures,
                                  model_cache=model_cache,
                                  sat_feature_mapping=sat_feature_mapping)


def print_blai_sat_feature_matrix(filename, matrix, state_ids, goals, names, complexities, types):
    ngoals = len(goals)
    nfeatures = len(names)
    num_numeric = sum(types)

    logging.info("Printing SAT matrix of {} features x {} states to '{}'".format(nfeatures, len(state_ids), filename))

    # Separate features into binary and numeric
    binary, numeric = [], []
    for i, t in enumerate(types, 0):
        if t == 0:
            binary.append(i)
        elif t == 1:
            numeric.append(i)
        else:
            assert 0, "Unknown feature type"

    all_feature_idxs = numeric + binary  # Shuffle the idxs so that all numeric features come first

    with open(filename, 'w') as f:
        # Header row: <#states> <#features> <#goals> <1 + index-last-numerical-feature> <index-first-boolean-feature>
        print("{} {} {} {} {}".format(len(state_ids), nfeatures, ngoals, num_numeric, num_numeric), file=f)

        # Line #2:: <#features> <list of feature names>
        print("{} {}".format(nfeatures, " ".join(names[f].replace(" ", "") for f in all_feature_idxs)), file=f)

        # Line #3: <#features> <list of feature costs>
        print("{} {}".format(nfeatures, " ".join(str(complexities[f]) for f in all_feature_idxs)), file=f)

        # Line #4: <# goal states> <list of goal state IDs>
        print("{} {}".format(ngoals, " ".join(map(str, goals))), file=f)

        # next lines: one per each state with format: <state-index> <#features-in-state> <list-features>
        # each feature has format: <feature-index>:<value>
        for s in state_ids:
            feature_values = ((i, matrix[s][f]) for i, f in enumerate(all_feature_idxs, 0))
            nonzero_features = [(i, v) for i, v in feature_values if v != 0]
            print("{} {} {}".format(s, len(nonzero_features), " ".join(
                "{}:{}".format(i, int(v)) for i, v in nonzero_features)), file=f)

    return all_feature_idxs  # A mapping between the new and the old feature indexes


def print_feature_matrix(filename, matrix, state_ids, goals, names, complexities):
    ngoals, nfeatures = len(goals), len(names)
    assert nfeatures == len(complexities) == matrix.shape[1]
    logging.info("Printing matrix of {} features x {} states to '{}'".format(nfeatures, len(state_ids), filename))

    with open(filename, 'w') as f:
        # Header row: <#states> <#features> <#goals>
        print("{} {} {}".format(len(state_ids), nfeatures, ngoals), file=f)

        # Line #2:: <list of feature names>
        print("{}".format(" ".join(name.replace(" ", "") for name in names)), file=f)

        # Line #3: <list of feature costs>
        print("{}".format(" ".join(str(c) for c in complexities)), file=f)

        # Line #4: <list of goal state IDs>
        print("{}".format(" ".join(map(str, goals))), file=f)

        # next lines: one per each state with format: <state-index> <#features-in-state> <list-features>
        # each feature has format: <feature-index>:<value>
        for s in state_ids:
            feature_values = ((i, matrix[s][i]) for i in range(0, nfeatures))
            nonzero_features = [(i, v) for i, v in feature_values if v != 0]
            print("{} {} {}".format(s, len(nonzero_features), " ".join(
                "{}:{}".format(i, int(v)) for i, v in nonzero_features)), file=f)
