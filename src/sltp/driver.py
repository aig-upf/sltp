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
import resource
import sys
from signal import signal, SIGPIPE, SIG_DFL
import numpy as np

from .version import get_version
from .returncodes import ExitCode
from .errors import CriticalPipelineError
from .util import console
from .util.bootstrap import setup_global_parser
from .util.command import execute, create_experiment_workspace
from .util.naming import compute_instance_tag, compute_experiment_tag, compute_serialization_name, \
    compute_maxsat_filename, compute_info_filename, compute_maxsat_variables_filename, compute_sample_filenames, \
    compute_test_sample_filenames
from .util.serialization import deserialize, serialize
from .util import performance

signal(SIGPIPE, SIG_DFL)

BASEDIR = os.path.dirname(os.path.dirname(os.path.dirname(os.path.realpath(__file__))))
BENCHMARK_DIR = os.path.join(BASEDIR, 'domains')
SAMPLE_DIR = os.path.join(BASEDIR, 'samples')
EXPDATA_DIR = os.path.join(BASEDIR, 'runs')


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
        raise InvalidConfigParameter('Parameter "{}" must be a {} integer value'.format(
            name, "positive " if positive else ""))


class Step:

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

    def description(self):
        raise NotImplementedError()

    def get_step_runner(self):
        raise NotImplementedError()


def _run_planner(config, data, rng):
    # Run the planner on all training and test instances
    def run(d, i, o, w, n):
        params = '-i {} --domain {} --driver {} --options=max_expansions={},width.max={}' \
            .format(i, d, config.driver, n, w)
        execute(command=[sys.executable, "run.py"] + params.split(' '),
                stdout=o, cwd=config.planner_location)

    for i, o, w in zip(config.instances, config.sample_files, config.max_width):
        run(config.domain, i, o, w, config.num_states)

    for i, o in zip(config.test_instances, config.test_sample_files):
        run(config.test_domain, i, o, -1, config.num_tested_states)

    return ExitCode.Success, dict()


class PlannerStep(Step):
    """ Run some planner on certain instance(s) to get the sample of transitions """

    VALID_DRIVERS = ("bfs", "ff")

    def get_required_attributes(self):
        return ["instances", "domain", "num_states", "planner_location", "driver", "test_instances", "num_tested_states"]

    def get_required_data(self):
        return []

    def process_config(self, config):
        if config["driver"] not in self.VALID_DRIVERS:
            raise InvalidConfigParameter('"driver" must be one of: {}'.format(self.VALID_DRIVERS))
        if any(not os.path.isfile(i) for i in config["instances"]):
            raise InvalidConfigParameter('"instances" must be the path to existing instance files. Specified: "{}"'
                                         .format(config["instances"]))
        if not os.path.isfile(config["domain"]):
            raise InvalidConfigParameter('Specified domain file "{}" does not exist'.format(config["domain"]))
        if not os.path.isdir(config["planner_location"]):
            raise InvalidConfigParameter('Specified planner location "{}" does not exist'
                                         .format(config["planner_location"]))
        check_int_parameter(config, "num_states", positive=True)

        config["instance_tag"] = compute_instance_tag(**config)
        config["experiment_tag"] = compute_experiment_tag(**config)
        config["experiment_dir"] = os.path.join(EXPDATA_DIR, config["experiment_tag"])

        mw = config.get("max_width", [-1] * len(config["instances"]))
        if isinstance(mw, int):
            mw = [mw]

        if len(config["instances"]) != len(mw):
            if len(mw) == 1:
                mw = mw * len(config["instances"])
            else:
                raise InvalidConfigParameter('"max_width" should have same length as "instances"')

        config["max_width"] = mw

        config["sample_files"] = compute_sample_filenames(**config)
        config["test_sample_files"] = compute_test_sample_filenames(**config)

        # TODO This should prob be somewhere else
        create_experiment_workspace(config["experiment_dir"], rm_if_existed=False)

        return config

    def description(self):
        return "Sampling of the state space"

    def get_step_runner(self):
        return _run_planner


class TransitionSamplingStep(Step):
    """ Generate the sample of transitions from the set of solved planning instances """
    def get_required_attributes(self):
        return ["instances", "sample_files", "experiment_dir", "optimal_selection_strategy"]

    def get_required_data(self):
        return []

    def process_config(self, config):
        config["resampled_states_filename"] = os.path.join(config["experiment_dir"], 'resampled.txt')
        config["num_sampled_states"] = config.get("num_sampled_states", None)
        config["complete_only_wrt_optimal"] = config.get("complete_only_wrt_optimal", False)
        config["sampling"] = config.get("sampling", "all")
        config["sat_transitions_filename"] = compute_info_filename(config, "sat_transitions.dat")
        config["goal_states_filename"] = compute_info_filename(config, "goal-states.dat")
        config["unsolvable_states_filename"] = compute_info_filename(config, "unsolvable-states.dat")
        config["transitions_filename"] = compute_info_filename(config, "transition-matrix.dat")

        ns = config["num_sampled_states"]
        if ns is not None:
            if isinstance(ns, int):
                ns = [ns]

            if len(config["instances"]) != len(ns):
                if len(ns) == 1:
                    ns = ns * len(config["instances"])
                else:
                    raise InvalidConfigParameter('"num_sampled_states" should have same length as "instances"')
            config["num_sampled_states"] = ns

        if config["sampling"] == "random" and config["num_sampled_states"] is None:
            raise InvalidConfigParameter('sampling="random" requires that option "num_sampled_states" is set')

        if config["sampling"] == "random" and config["complete_only_wrt_optimal"]:
            raise InvalidConfigParameter('sampling="random" disallows the use of option "complete_only_wrt_optimal"')

        return config

    def description(self):
        return "Generation of the training sample"

    def get_step_runner(self):
        from . import sampling
        return sampling.run


class ConceptGenerationStep(Step):
    """ Generate systematically a set of features of bounded complexity from the transition (state) sample """
    def get_required_attributes(self):
        return ["domain", "experiment_dir", "max_concept_size"]

    def get_required_data(self):
        return ["sample"]

    def process_config(self, config):
        check_int_parameter(config, "max_concept_size")

        config["concept_dir"] = os.path.join(config["experiment_dir"], 'terms')
        config["concept_generator"] = config.get("concept_generator", None)
        config["feature_generator"] = config.get("feature_generator", None)
        config["enforce_features"] = config.get("enforce_features", None)
        config["parameter_generator"] = config.get("parameter_generator", None)
        config["distance_feature_max_complexity"] = config.get("distance_feature_max_complexity", None)
        config["max_concept_grammar_iterations"] = config.get("max_concept_grammar_iterations", None)
        config["concept_denotation_filename"] = compute_info_filename(config, "concept-denotations.txt")
        config["role_denotation_filename"] = compute_info_filename(config, "role-denotations.txt")

        return config

    def description(self):
        return "Generation of concepts"

    def get_step_runner(self):
        from . import features
        return features.run


def setup_feature_generation_filenames(config):
    config["feature_matrix_filename"] = compute_info_filename(config, "feature-matrix.dat")
    config["bin_feature_matrix_filename"] = compute_info_filename(config, "feature-matrix-bin.dat")
    config["feature_complexity_filename"] = compute_info_filename(config, "feature-complexity.dat")
    config["feature_names_filename"] = compute_info_filename(config, "feature-names.dat")


class FeatureMatrixGenerationStep(Step):
    """ Generate and output the feature and transition matrices for the problem  """
    def get_required_attributes(self):
        return ["experiment_dir", "prune_positive_features", "prune_infty_features"]

    def process_config(self, config):
        setup_feature_generation_filenames(config)
        config["feature_filename"] = compute_info_filename(config, "features.txt")
        config["sat_feature_matrix_filename"] = compute_info_filename(config, "sat_matrix.dat")
        config["feature_info_filename"] = compute_info_filename(config, "feature-info.dat")
        config["feature_denotation_filename"] = compute_info_filename(config, "feature-denotations.txt")
        return config

    def get_required_data(self):
        return ["features", "model_cache", "sample"]

    def description(self):
        return "Generation of the feature and transition matrices"

    def get_step_runner(self):
        from .matrices import generate_features
        return generate_features


class CPPFeatureGenerationStep(Step):
    """ Generate exhaustively a set of all features up to a given complexity from the transition (state) sample """
    def get_required_attributes(self):
        return ["domain", "experiment_dir", "max_concept_size"]

    def get_required_data(self):
        return ["sample"]

    def process_config(self, config):
        check_int_parameter(config, "max_concept_size")

        setup_feature_generation_filenames(config)

        config["concept_dir"] = os.path.join(config["experiment_dir"], 'terms')
        config["concept_generator"] = config.get("concept_generator", None)
        config["feature_generator"] = config.get("feature_generator", None)
        config["enforce_features"] = config.get("enforce_features", None)
        config["parameter_generator"] = config.get("parameter_generator", None)
        config["distance_feature_max_complexity"] = config.get("distance_feature_max_complexity", None)
        config["concept_denotation_filename"] = compute_info_filename(config, "concept-denotations.txt")
        config["sat_feature_matrix_filename"] = compute_info_filename(config, "sat_matrix.dat")
        config["feature_denotation_filename"] = compute_info_filename(config, "feature-denotations.txt")
        config["serialized_feature_filename"] = compute_info_filename(config, "serialized-features.io")

        if config["enforce_features"]:
            raise RuntimeError("Option enforce_features not allowed when using the C++ feature generator")

        return config

    def description(self):
        return "C++ feature generation module"

    def get_step_runner(self):
        from . import featuregen
        return featuregen.run


class MaxsatProblemGenerationStep(Step):
    """ Generate the max-sat problem from a given set of generated features """
    def get_required_attributes(self):
        return ["experiment_dir"]

    def process_config(self, config):
        config["cnf_filename"] = compute_maxsat_filename(config)
        config["maxsat_variables_file"] = compute_maxsat_variables_filename(config)
        return config

    def get_required_data(self):
        return ["sample", "enforced_feature_idxs", "in_goal_features"]

    def description(self):
        return "Generation of the max-sat problem"

    def get_step_runner(self):
        from . import learn_actions
        return learn_actions.generate_maxsat_problem


class MaxsatProblemSolutionStep(Step):
    """ Run some max-sat solver on the generated encoding """
    def get_required_attributes(self):
        return ['maxsat_solver', 'maxsat_timeout']

    def get_required_data(self):
        return []

    def description(self):
        return "Solution of the max-sat problem"

    def get_step_runner(self):
        from . import learn_actions
        return learn_actions.run_solver


class SatProblemGenerationStep(Step):
    """ Call Blai's SAT-encoding generator from a given set of generated features """
    def get_required_attributes(self):
        return ["experiment_dir", "encoding_k", "encoding_m"]

    def process_config(self, config):
        config["sat_theory_prefix"] = compute_info_filename(config, "sat")
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


class ActionModelStep(Step):
    """ Generate an abstract action model from the solution of the max-sat encoding """

    def get_required_attributes(self):
        return ['serialized_feature_filename']

    def process_config(self, config):
        config["feature_namer"] = config.get("feature_namer")
        config["qnp_abstraction_filename"] = compute_info_filename(config, "abstraction.qnp")
        return config

    def get_required_data(self):
        return ["cnf_translator", "cnf_solution", "sample"]

    def description(self):
        return "Computation of the action model"

    def get_step_runner(self):
        from . import learn_actions
        return learn_actions.compute_action_model


class ActionModelFromFeatureIndexesStep(Step):
    """ Generate an abstract action model from the solution of the max-sat encoding """

    def get_required_attributes(self):
        return []

    def process_config(self, config):
        config["feature_namer"] = config.get("feature_namer")
        config["qnp_abstraction_filename"] = compute_info_filename(config, "abstraction.qnp")
        return config

    def get_required_data(self):
        return ["selected_feature_idxs", "sample", "sat_feature_mapping"]

    def description(self):
        return "Computation of the action model from the feature indexes"

    def get_step_runner(self):
        from . import learn_actions
        return learn_actions.compute_action_model_from_feature_idxs


class QNPGenerationStep(Step):
    """ Generate a QNP encoding from the computed abstract action model """

    def get_required_attributes(self):
        return ["qnp_abstraction_filename"]

    def process_config(self, config):
        config["qnp_prefix"] = "encoding"
        return config

    def get_required_data(self):
        return []

    def description(self):
        return "Generation of the QNP encoding"

    def get_step_runner(self):
        from . import qnp
        return qnp.generate_encoding


class InhouseMaxsatSolverStep(Step):
    """ Generate and solve the maxsat theory with Blai's solver """
    def get_required_attributes(self):
        return ["experiment_dir"]

    def process_config(self, config):
        config["sat_theory_prefix"] = compute_info_filename(config, "sat")
        config["maxsat_solver_out"] = compute_info_filename(config, "maxsat.out")
        return config

    def get_required_data(self):
        return ["sample"]

    def description(self):
        return "Generation & (approximate) solving of the max-sat problem with Blai's solver"

    def get_step_runner(self):
        from . import approxmaxsat
        return approxmaxsat.run
    

def generate_pipeline(pipeline="maxsat", **kwargs):
    pipeline, config = generate_pipeline_from_list(PIPELINES[pipeline], **kwargs)
    return pipeline


def generate_pipeline_from_list(elements, **kwargs):
    steps = []
    config = kwargs
    for klass in elements:
        step = klass(**config)
        config = step.config
        steps.append(step)
    return steps, config


def save(basedir, output):
    if not output:
        return

    def serializer():
        return tuple(serialize(data, compute_serialization_name(basedir, name)) for name, data in output.items())

    console.log_time(serializer, logging.DEBUG,
                     'Serializing data elements "{}" to directory "{}"'.format(', '.join(output.keys()), basedir))


def _deserializer(basedir, items):
    return dict((k, deserialize(compute_serialization_name(basedir, k))) for k in items)


def load(basedir, items):
    def deserializer():
        return _deserializer(basedir, items)

    output = console.log_time(deserializer, logging.DEBUG,
                              'Deserializing data elements "{}" from directory "{}"'.format(', '.join(items), basedir))
    return output


class StepRunner:
    """ Run the given step """
    def __init__(self, stepnum, step_name, target, required_data):
        self.start = self.elapsed_time()
        self.target = target
        self.stepnum = stepnum
        self.step_name = step_name
        self.required_data = required_data
        self.loglevel = None

    def elapsed_time(self):
        info_children = resource.getrusage(resource.RUSAGE_CHILDREN)
        info_self = resource.getrusage(resource.RUSAGE_SELF)
        # print("({}) Self: {}".format(os.getpid(), info_self))
        # print("({}) Children: {}".format(os.getpid(), info_children))
        return info_children.ru_utime + info_children.ru_stime + info_self.ru_utime + info_self.ru_stime

    def used_memory(self):
        return performance.memory_usage()

    def setup(self, quiet):
        self.loglevel = logging.getLogger().getEffectiveLevel()
        if quiet:
            logging.getLogger().setLevel(logging.ERROR)
        else:
            print(console.header("(pid: {}) STARTING STEP #{}: {}".format(os.getpid(), self.stepnum, self.step_name)))

    def teardown(self, quiet):
        if quiet:
            logging.getLogger().setLevel(self.loglevel)
        else:
            current = self.elapsed_time()
            print(console.header("END OF STEP #{}: {}. {:.2f} CPU sec - {:.2f} MB".format(
                self.stepnum, self.step_name, current - self.start, self.used_memory())))

    def run(self, config):
        exitcode = self._run(config)
        self.teardown(config.quiet)
        return exitcode

    def _run(self, config):
        """ Run the StepRunner target.
            This method will also be the entry point of spawned subprocess in the case
            of SubprocessStepRunners
        """
        self.setup(config.quiet)

        data = Bunch(load(config.experiment_dir, self.required_data)) if self.required_data else None
        rng = np.random.RandomState(config.random_seed)  # ATM we simply create a RNG in each subprocess

        try:
            exitcode, output = self.target(config=config, data=data, rng=rng)
        except Exception as exception:
            # Flatten the exception so that we make sure it can be serialized,
            # and return it immediately so that it can be reported from the parent process
            logging.error("Critical error in the pipeline")
            import traceback
            traceback.print_exception(None, exception, exception.__traceback__)
            raise CriticalPipelineError("Error: {}".format(str(exception)))

        save(config.experiment_dir, output)
        # profiling.start()
        return exitcode


class SubprocessStepRunner(StepRunner):
    """ Run the given step by spawning a subprocess and waiting for its finalization """
    def _run(self, config):
        pool = multiprocessing.Pool(processes=1)
        result = pool.apply_async(StepRunner._run, (), {"self": self, "config": config})
        exitcode = result.get()
        pool.close()
        pool.join()
        return exitcode

    def used_memory(self):
        info_children = resource.getrusage(resource.RUSAGE_CHILDREN)
        # print("({}) Children: {}".format(os.getpid(), info_children))
        return info_children.ru_maxrss / 1024.0


class SATStateFactorizationStep(Step):
    """ Generate the SAT problem that learns a state factorization """
    def get_required_attributes(self):
        return ["experiment_dir"]

    def process_config(self, config):
        config["cnf_filename"] = compute_maxsat_filename(config)
        config["maxsat_variables_file"] = compute_maxsat_variables_filename(config)
        return config

    def get_required_data(self):
        return ["goal_states", "transitions", "state_ids", "enforced_feature_idxs", "optimal_transitions"]

    def description(self):
        return "SAT encoding to learn a state factorization "

    def get_step_runner(self):
        from . import factorization
        return factorization.learn_factorization


# class DFAGenerationStep(Step):
#     """ Generate the DFA from observation traces """
#
#     def get_required_attributes(self):
#         return ["experiment_dir"]
#
#     def process_config(self, config):
#         config["cnf_filename"] = compute_maxsat_filename(config)
#         config["maxsat_variables_file"] = compute_maxsat_variables_filename(config)
#         return config
#
#     def get_required_data(self):
#         return ["goal_states", "transitions", "state_ids", "enforced_feature_idxs", "optimal_transitions"]
#
#     def description(self):
#         return "DFA generation from observation traces"
#
#     def get_step_runner(self):
#         from . import factorization
#         return None
#         # return learn_actions.generate_maxsat_problem


class AbstractionTestingComputation(Step):
    """  """
    def get_required_attributes(self):
        return ["experiment_dir", "test_instances", "test_domain"]

    def process_config(self, config):
        if config["test_domain"] is not None:
            if not os.path.isfile(config["test_domain"]):
                raise InvalidConfigParameter('"test_domain" must be either None or the path to an existing domain file')
            if any(not os.path.isfile(i) for i in config["test_instances"]):
                raise InvalidConfigParameter('"test_instances" must point to existing files')

        return config

    def get_required_data(self):
        return ["abstract_actions", "selected_features"]

    def description(self):
        return "Testing of the computed abstraction on the testing instances"

    def get_step_runner(self):
        from . import tester
        return tester.run


class Experiment:
    def __init__(self, steps, parameters):
        self.args = None
        self.steps = steps

    def print_step_description(self):
        return "\t\t" + "\n\t\t".join("{}. {}".format(i, s.description()) for i, s in enumerate(self.steps, 1))

    def hello(self, args=None):
        print(console.header("SLTP v.{}".format(get_version())))
        argparser = setup_global_parser(step_description=self.print_step_description())
        self.args = argparser.parse_args(args)
        if not self.args.steps and not self.args.run_all_steps:
            argparser.print_help()
            sys.exit(0)

    def run(self, args=None):
        self.hello(args)
        # If no steps were given on the commandline, run all exp steps.
        steps = [get_step(self.steps, name) for name in self.args.steps] or self.steps

        for stepnum, step in enumerate(steps, start=1):
            try:
                run_and_check_output(step, stepnum, SubprocessStepRunner)
            except Exception as e:
                logging.error((step, _create_exception_msg(step, e)))
                return e
        return ExitCode.Success


def run_and_check_output(step, stepnum, runner_class, raise_on_error=True):
    runner = runner_class(stepnum, step.description(), step.get_step_runner(), step.get_required_data())
    exitcode = runner.run(config=Bunch(step.config))
    if raise_on_error and exitcode is not ExitCode.Success:
        raise RuntimeError(_create_exception_msg(step, exitcode))
    return exitcode


def _create_exception_msg(step, e):
    return 'Critical error processing step "{}". Error message: {}'.\
        format(step.description(), e)


def run_experiment(experiment, argv):
    retcode = experiment.run(argv)
    if retcode != ExitCode.Success:
        sys.exit(-1)
    

class Bunch:
    def __init__(self, adict):
        self.__dict__.update(adict)

    def to_dict(self):
        return self.__dict__.copy()


PIPELINES = dict(
    maxsat=[
        PlannerStep,
        TransitionSamplingStep,
        ConceptGenerationStep,
        FeatureMatrixGenerationStep,
        MaxsatProblemGenerationStep,
        MaxsatProblemSolutionStep,
        ActionModelStep,
        AbstractionTestingComputation,
        # QNPGenerationStep,
    ],
    maxsatcpp=[
        PlannerStep,
        TransitionSamplingStep,
        CPPFeatureGenerationStep,
        MaxsatProblemGenerationStep,
        MaxsatProblemSolutionStep,
        ActionModelStep,
        AbstractionTestingComputation,
        # QNPGenerationStep,
    ],
    maxsat_poly=[
        PlannerStep,
        TransitionSamplingStep,
        CPPFeatureGenerationStep,
        # ConceptGenerationStep,
        # FeatureMatrixGenerationStep,
        InhouseMaxsatSolverStep,  # Blai's polynomial ad-hoc maxsat algorithm
        ActionModelFromFeatureIndexesStep,
        AbstractionTestingComputation,
    ],
    sat=[
        PlannerStep,
        ConceptGenerationStep,
        FeatureMatrixGenerationStep,
        SatProblemGenerationStep,
        SatProblemSolutionStep,
        SatSolutionDecodingStep,
    ],
    # observations=[
    #     DFAGenerationStep,
    #     SATStateFactorizationStep,
    # ],
)
