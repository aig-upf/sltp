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
import copy
import logging
import multiprocessing
import os
import sys
from signal import signal, SIGPIPE, SIG_DFL

from errors import CriticalPipelineError
from util import console
from util.bootstrap import setup_global_parser
from util.command import execute
from util.console import print_header, log_time
from util.naming import compute_instance_tag, compute_experiment_tag, compute_serialization_name, \
    compute_maxsat_filename, compute_info_filename, compute_maxsat_variables_filename
from util.serialization import deserialize, serialize

signal(SIGPIPE, SIG_DFL)

BASEDIR = os.path.dirname(os.path.realpath(__file__))
VERSION = "0.2"

BASE_DIR = os.path.dirname(os.path.abspath(__file__))
BENCHMARK_DIR = os.path.join(BASE_DIR, 'domains')
SAMPLE_DIR = os.path.join(BASE_DIR, 'samples')
EXPDATA_DIR = os.path.join(BASE_DIR, 'runs')


class InvalidConfigParameter(Exception):
    def __init__(self, msg=None):
        msg = msg or 'Invalid configuration parameter'
        super().__init__(msg)


def _get_step_index(steps, step_name):
    for index, step in enumerate(steps):
        if step.name == step_name:
            return index
    logging.critical('There is no step called "{}"'.format(step_name))


def get_step(steps, step_name):
    """*step_name* can be a step's name or number."""
    if step_name.isdigit():
        try:
            return steps[int(step_name) - 1]
        except IndexError:
            logging.critical('There is no step number {}'.format(step_name))
    return steps[_get_step_index(steps, step_name)]


def check_int_parameter(config, name, positive=False):
    try:
        config[name] = int(config[name])
        if positive and config[name] <= 0:
            raise ValueError()
    except ValueError:
        raise InvalidConfigParameter('Parameter "{}" must be a {}integer value'.format(
            name, "positive " if positive else ""))


class Step(object):

    def __init__(self, **kwargs):
        self.config = self.process_config(self.parse_config(**kwargs))

    def process_config(self, config):
        return config  # By default, we do nothing

    def get_required_attributes(self):
        raise NotImplementedError()

    def get_required_data(self):
        raise NotImplementedError()

    def parse_config(self, **kwargs):
        config = copy.deepcopy(kwargs)
        for attribute in self.get_required_attributes():
            if attribute not in kwargs:
                raise RuntimeError('Missing attribute "{}" in step "{}"'.format(attribute, self.__class__.__name__))
            config[attribute] = kwargs[attribute]

        return config

    def run(self):
        req_data = self.get_required_data()
        data = None
        if req_data:
            data = Bunch(load(self.config["experiment_dir"], self.get_required_data()))
        runner = StepRunner(self.description(), self.get_step_runner())
        result = runner.run(config=Bunch(self.config), data=data)
        return result

    def description(self):
        raise NotImplementedError()

    def get_step_runner(self):
        raise NotImplementedError()


def _run_planner(config, data):
    params = '--instance {} --domain {} --driver {} --disable-static-analysis --options="max_expansions={}"'\
        .format(config.instance, config.domain, config.driver, config.num_states)
    execute(command=[sys.executable, "run.py"] + params.split(' '),
            stdout=config.sample_file,
            cwd=config.planner_location
            )
    return dict()


class PlannerStep(Step):
    """ Run some planner on certain instance(s) to get the sample of transitions """

    VALID_DRIVERS = ("bfs", "ff", "iw2")

    def __init__(self, **kwargs):
        super().__init__(**kwargs)

    def get_required_attributes(self):
        return ["instance", "domain", "num_states", "planner_location", "driver"]

    def get_required_data(self):
        return []

    def process_config(self, config):
        if config["driver"] not in self.VALID_DRIVERS:
            raise InvalidConfigParameter('"driver" must be one of: {}'.format(self.VALID_DRIVERS))
        if not os.path.isfile(config["instance"]):
            raise InvalidConfigParameter('"instance" must be the path to an existing instance file')
        if not os.path.isfile(config["domain"]):
            raise InvalidConfigParameter('"domain" must be the path to an existing domain file')
        if not os.path.isdir(config["planner_location"]):
            raise InvalidConfigParameter('"planner_location" must be the path to the actual planner')
        check_int_parameter(config, "num_states", positive=True)
        config["instance_tag"] = compute_instance_tag(**config)
        config["experiment_tag"] = compute_experiment_tag(**config)
        config["experiment_dir"] = os.path.join(EXPDATA_DIR, config["experiment_tag"])
        config["sample_file"] = os.path.join(config["experiment_dir"], "samples.txt")

        # TODO This should prob be somewhere else:
        os.makedirs(config["experiment_dir"], exist_ok=True)

        return config

    def description(self):
        return "Sampling of the state space"

    def get_step_runner(self):
        return _run_planner


class ConceptGenerationStep(Step):
    """ Generate systematically a set of features of bounded complexity from the transition (state) sample """
    def __init__(self, **kwargs):
        super().__init__(**kwargs)

    def get_required_attributes(self):
        return ["sample_file", "domain", "experiment_dir", "max_concept_size"]

    def get_required_data(self):
        return []

    def process_config(self, config):
        check_int_parameter(config, "max_concept_size")

        config["concept_dir"] = os.path.join(config["experiment_dir"], 'terms')
        config["feature_stdout"] = os.path.join(config["experiment_dir"], 'feature-generation.stdout.txt')
        config["concept_generator"] = config.get("concept_generator", None)
        config["feature_generator"] = config.get("feature_generator", None)
        config["parameter_generator"] = config.get("parameter_generator", None)
        config["use_distance_features"] = config.get("use_distance_features", False)
        config["max_concept_grammar_iterations"] = config.get("max_concept_grammar_iterations", None)

        return config

    def description(self):
        return "Generation of concepts"

    def get_step_runner(self):
        import features
        return features.run


class FeatureMatrixGenerationStep(Step):
    """ Generate and output the feature and transition matrices for the problem  """
    def __init__(self, **kwargs):
        super().__init__(**kwargs)

    def get_required_attributes(self):
        return ["experiment_dir"]

    def process_config(self, config):
        config["feature_filename"] = compute_info_filename(config, "features.txt")
        config["feature_matrix_filename"] = compute_info_filename(config, "feature-matrix.dat")
        config["bin_feature_matrix_filename"] = compute_info_filename(config, "feature-matrix-bin.dat")
        config["feature_complexity_filename"] = compute_info_filename(config, "feature-complexity.dat")
        config["feature_names_filename"] = compute_info_filename(config, "feature-names.dat")
        config["transitions_filename"] = compute_info_filename(config, "transition-matrix.dat")
        config["goal_states_filename"] = compute_info_filename(config, "goal-states.dat")
        config["sat_transitions_filename"] = compute_info_filename(config, "sat_transitions.dat")
        config["sat_feature_matrix_filename"] = compute_info_filename(config, "sat_matrix.dat")
        config["feature_info_filename"] = compute_info_filename(config, "feature-info.dat")
        config["feature_denotation_filename"] = compute_info_filename(config, "feature-denotations.txt")
        return config

    def get_required_data(self):
        return ["features", "extensions", "states", "goal_states", "transitions"]

    def description(self):
        return "Generation of the feature and transition matrices"

    def get_step_runner(self):
        from matrices import generate_features
        return generate_features


class MaxsatProblemGenerationStep(Step):
    """ Generate the max-sat problem from a given set of generated features """
    def __init__(self, **kwargs):
        super().__init__(**kwargs)

    def get_required_attributes(self):
        return ["experiment_dir"]

    def process_config(self, config):
        config["cnf_filename"] = compute_maxsat_filename(config)
        config["maxsat_variables_file"] = compute_maxsat_variables_filename(config)
        config["relax_numeric_increase"] = config.get("relax_numeric_increase", False)
        return config

    def get_required_data(self):
        return ["goal_states", "transitions", "state_ids"]

    def description(self):
        return "Generation of the max-sat problem"

    def get_step_runner(self):
        import learn_actions
        return learn_actions.generate_maxsat_problem


class SatProblemGenerationStep(Step):
    """ Call Blai's SAT-encoding generator from a given set of generated features """
    def __init__(self, **kwargs):
        super().__init__(**kwargs)

    def get_required_attributes(self):
        return ["experiment_dir", "encoding_k", "encoding_m"]

    def process_config(self, config):
        config["sat_theory_prefix"] = os.path.join(config["experiment_dir"], "sat")
        return config

    def get_required_data(self):
        return ["features", "state_ids", "goal_states"]

    def description(self):
        return "Generation of the (alternative sat encoding) problem"

    def get_step_runner(self):
        from compact_encoding import encoder
        return encoder.encode


class SatProblemSolutionStep(Step):
    """ Call some SAT solver to solve Blai's SAT encoding """
    def __init__(self, **kwargs):
        super().__init__(**kwargs)

    def get_required_attributes(self):
        return ["experiment_dir", "encoding_k", "encoding_m", "sat_theory_prefix"]

    def process_config(self, config):
        config["sat_theory_filename"] = "{}_{}_{}_theory.cnf".format(
            config["sat_theory_prefix"], config["encoding_k"], config["encoding_m"])
        config["sat_solution_filename"] = "{}_{}_{}_model.cnf".format(
            config["sat_theory_prefix"], config["encoding_k"], config["encoding_m"])
        return config

    def get_required_data(self):
        return ["features", "state_ids", "goal_states"]

    def description(self):
        return "Solution of the SAT problem"

    def get_step_runner(self):
        from compact_encoding import solver
        return solver.run


class SatSolutionDecodingStep(Step):
    """ Decode the SAT solution from Blai's encoding """
    def __init__(self, **kwargs):
        super().__init__(**kwargs)

    def get_required_attributes(self):
        return ["experiment_dir", "encoding_k", "encoding_m", "sat_theory_prefix", "sat_solution_filename"]

    def process_config(self, config):
        return config

    def get_required_data(self):
        return ["features", "state_ids", "goal_states"]

    def description(self):
        return "Decoding of the SAT solution"

    def get_step_runner(self):
        from compact_encoding import encoder
        return encoder.decode


class MaxsatProblemSolutionStep(Step):
    """ Run some max-sat solver on the generated encoding """

    def __init__(self, **kwargs):
        super().__init__(**kwargs)

    def get_required_attributes(self):
        return []

    def get_required_data(self):
        return ["cnf_translator"]

    def description(self):
        return "Solution of the max-sat problem"

    def get_step_runner(self):
        import learn_actions
        return learn_actions.run_solver


def default_feature_namer(s):
    return str(s)


class ActionModelStep(Step):
    """ Generate an abstract action model from the solution of the max-sat encoding """

    def __init__(self, **kwargs):
        super().__init__(**kwargs)

    def get_required_attributes(self):
        return []

    def process_config(self, config):
        config["feature_namer"] = config.get("feature_namer", default_feature_namer) or default_feature_namer
        return config

    def get_required_data(self):
        return ["cnf_translator", "cnf_solution"]

    def description(self):
        return "Computation of the action model"

    def get_step_runner(self):
        import learn_actions
        return learn_actions.compute_action_model


PIPELINES = dict(
    maxsat=[
        PlannerStep,
        ConceptGenerationStep,
        FeatureMatrixGenerationStep,
        MaxsatProblemGenerationStep,
        MaxsatProblemSolutionStep,
        ActionModelStep,
    ],
    sat=[
        PlannerStep,
        ConceptGenerationStep,
        FeatureMatrixGenerationStep,
        SatProblemGenerationStep,
        SatProblemSolutionStep,
        SatSolutionDecodingStep,
    ],
    heuristic=[
        PlannerStep,
        ConceptGenerationStep,
        FeatureMatrixGenerationStep,
        # ...
    ],
)


def generate_pipeline(pipeline="maxsat", **kwargs):
    steps = []
    config = kwargs
    for klass in PIPELINES[pipeline]:
        step = klass(**config)
        config = step.config
        steps.append(step)
    return steps


class Experiment(object):
    def __init__(self, steps):
        self.args = None
        self.steps = steps

    def print_step_description(self):
        return "\t\t" + "\n\t\t".join("{}. {}".format(i, s.description()) for i, s in enumerate(self.steps, 1))

    def run(self):
        print_header("Generalized Action Model Learner, v.{}".format(VERSION))
        argparser = setup_global_parser(step_description=self.print_step_description())
        self.args = argparser.parse_args()
        if not self.args.steps and not self.args.run_all_steps:
            argparser.print_help()
            sys.exit(0)

        # If no steps were given on the commandline, run all exp steps.
        steps = [get_step(self.steps, name) for name in self.args.steps] or self.steps

        for step in steps:
            result = step.run()
            if result is not None:
                logging.error('Critical error while processing step "{}". Execution will be terminated. '
                              'Error message:'.format(step.description()))
                logging.error("\t{}".format(result))
                break


def save(basedir, output):
    if not output:
        return

    def serializer():
        return tuple(serialize(data, compute_serialization_name(basedir, name)) for name, data in output.items())

    log_time(serializer,
             'Serializing data elements "{}" to directory "{}"'.format(', '.join(output.keys()), basedir))


def _deserializer(basedir, items):
    return dict((k, deserialize(compute_serialization_name(basedir, k))) for k in items)


def load(basedir, items):
    def deserializer():
        return _deserializer(basedir, items)

    output = log_time(deserializer,
                      'Deserializing data elements "{}" from directory "{}"'.format(', '.join(items), basedir))
    return output


class StepRunner(object):
    def __init__(self, step_name, target):
        self.target = target
        self.step_name = step_name

    def run(self, **kwargs):
        pool = multiprocessing.Pool(processes=1)
        result = pool.apply_async(self._runner, (), kwargs)
        res = result.get()
        pool.close()
        pool.join()
        return res

    def _runner(self, config, data):
        # import tracemalloc
        # tracemalloc.start()
        # memutils.display_top(tracemalloc.take_snapshot())
        from util import performance

        result = None
        console.print_header("STARTING STEP: {}".format(self.step_name))
        start = performance.timer()
        try:
            output = self.target(config=config, data=data)
        except CriticalPipelineError as exc:
            output = dict()
            result = exc
        save(config.experiment_dir, output)
        performance.print_performance_stats(self.step_name, start)
        # profiling.start()
        return result


class Bunch(object):
    def __init__(self, adict):
        self.__dict__.update(adict)
