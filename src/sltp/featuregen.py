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

from sltp.matrices import NP_FEAT_VALUE_TYPE, cast_feature_value_to_numpy_value
from sltp.models import DLModel
from tarski.dl import PrimitiveConcept, UniversalConcept

from .features import parse_all_instances, compute_models, InstanceInformation
from .util.command import execute
from .returncodes import ExitCode


def run(config, data, rng):
    return extract_features(config, data.sample)


def convert_tuple(universe, name, values):
    return (name, ) + tuple(universe.value(i) for i in (values if isinstance(values, tuple) else (values, )))


def serialize_dl_model(model: DLModel, info: InstanceInformation):
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


def print_sample_info(sample, infos, model_cache, all_predicates, all_functions, all_objects, workspace):

    state_fn = os.path.join(workspace, "states.io")
    sample_fn = os.path.join(workspace, "sample.io")

    logging.info("Printing state data to {}".format(state_fn))
    all_atoms = set()  # We will accumulate all possible atoms here
    with open(state_fn, "w") as f:
        for expected_id, (id_, state) in enumerate(sample.states.items(), 0):
            # sample.states is an ordered dict. All ids are consecutive, so no need to print them,
            # but we check just in case
            assert expected_id == id_
            state_instance = sample.instance[id_]  # The instance to which the state belongs
            full_state = serialize_dl_model(model_cache.models[id_], infos[state_instance])

            # print one line per state with all state atoms, e.g. at,bob,shed   at,spanner,location;
            print("\t".join(",".join(atom) for atom in full_state), file=f)

            all_atoms.update(full_state)

    logging.info("Printing sample information to {}".format(sample_fn))
    with open(sample_fn, "w") as f:
        print("dummy-sample-name", file=f)  # First line: sample name

        print(" ".join("{}/{}".format(name, arity) for name, arity in sorted(all_predicates)),
              file=f)  # Second line: all predicate names

        print(" ".join("{}/{}".format(name, arity) for name, arity in sorted(all_functions)),
              file=f)  # # Third line: all function names

        print(" ".join(sorted(all_objects)), file=f)  # Fourth line: all object names

        # Fifth line: all possible atoms, i.e. grounded predicates (possibly from different problem instances):
        print("\t".join(",".join(atom) for atom in sorted(all_atoms)), file=f)


def transform_generator_output(config, matrix_filename, info_filename):
    logging.info("Transforming generator output to numpy format...")
    import numpy as np

    logging.info("Reading feature information from {}".format(info_filename))
    with open(info_filename, "r") as f:
        # Line  # 1: feature names
        names = f.readline().rstrip().split("\t")

        # Line #2: feature complexities
        complexities = [int(x) for x in f.readline().rstrip().split("\t")]

        # Line #3: feature types (0: boolean; 1: numeric)
        types = [bool(x) for x in f.readline().rstrip().split("\t")]

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

    bin_matrix = np.array(matrix, dtype=np.bool)
    np.save(config.feature_matrix_filename, matrix)
    np.save(config.bin_feature_matrix_filename, bin_matrix)
    np.save(config.feature_complexity_filename, np.array(complexities, dtype=NP_FEAT_VALUE_TYPE))
    np.save(config.feature_names_filename, np.array(names, dtype=np.unicode_))
    # np.savez(config.feature_matrix_filename, matrix)


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

    all_goal_predicates = set(itertools.chain.from_iterable(info.goal_predicates for info in infos))
    if any(all_goal_predicates != info.goal_predicates for info in infos):
        logging.warning("Not all instances in the training set use the same goal predicate")

    all_objects = set()
    all_predicates, all_functions = set(), set()
    for problem, lang, _ in parsed_problems:
        all_objects.update(c.symbol for c in lang.constants())
        all_predicates.update(set((p.symbol, p.arity) for p in lang.predicates if not p.builtin))
        all_functions.update(set((p.symbol, p.arity) for p in lang.functions if not p.builtin))

        # Add type predicates as well
        all_predicates.update(set((p.name, 1) for p in lang.sorts if not p.builtin and p != lang.Object))

    logging.info('Invoking C++ feature generation module'.format())

    # Write sample information
    print_sample_info(sample, infos, model_cache, all_predicates, all_functions,
                      all_objects, config.experiment_dir)

    # Invoke C++ feature generation module
    cmd = os.path.realpath(os.path.join(config.featuregen_location, "featuregen"))
    args = [str(config.max_concept_size), config.experiment_dir]
    generate_debug_scripts(config.experiment_dir, cmd, args)
    retcode = execute([cmd] + args)

    if retcode != 0:
        return ExitCode.FeatureGenerationUnknownError, dict()

    # Read off the output of the module and transform it into the numpy matrices to be consumed
    # by the next pipeline step
    transform_generator_output(config,
                               os.path.join(config.experiment_dir, "feature-matrix.io"),
                               os.path.join(config.experiment_dir, "feature-info.io"),)

    return ExitCode.Success, dict(enforced_feature_idxs=[])

