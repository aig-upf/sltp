import itertools
import logging
import os

import numpy as np
from sltp.util import console

from .util.naming import compute_info_filename
from .util.tools import load_selected_features, process_features
from tarski.dl import compute_dl_vocabulary

from .validator import AbstractionValidator

from .sampling import log_sampled_states, print_transition_matrices
from .features import create_model_cache_from_samples, parse_all_instances, compute_instance_information, \
    report_use_goal_denotation
from .returncodes import ExitCode
from .util.misc import update_dict
from .driver import Experiment, generate_pipeline_from_list, PlannerStep, check_int_parameter, \
    InvalidConfigParameter, TransitionSamplingStep, SubprocessStepRunner, \
    run_and_check_output, save, MaxsatProblemSolutionStep, \
    load, save, Bunch, CPPFeatureGenerationStep, CPPMaxsatProblemGenerationStep, CPPActionModelStep, AbstractionTestingStep


def initial_sample_selection(sample, config, rng):
    # Get the initial sample of M states plus one goal state per instance, if any goal is available.
    # plus all states on one path to each goal.
    # plus all instance roots
    # all_state_ids_shuffled = list(sample.states.keys())
    # rng.shuffle(all_state_ids_shuffled)
    expanded_state_ids_shuffled = list(sample.expanded)
    rng.shuffle(expanded_state_ids_shuffled)
    initial_size = min(len(expanded_state_ids_shuffled), config.initial_sample_size-len(sample.roots))
    working_sample_idxs = set(expanded_state_ids_shuffled[:initial_size])
    working_sample_idxs.update(sample.get_one_goal_per_instance())
    working_sample_idxs.update(itertools.chain.from_iterable(sample.optimal_transitions))  # include all optimal tx
    working_sample_idxs.update(sample.roots)
    return working_sample_idxs, expanded_state_ids_shuffled


def full_learning(config, sample, rng):
    working_sample_idxs, expanded_state_ids_shuffled = initial_sample_selection(sample, config, rng)
    working_sample = sample.resample(working_sample_idxs)

    use_goal_denotation = report_use_goal_denotation(config.parameter_generator)
    parsed_problems = parse_all_instances(config.domain, config.instances)
    infos = [compute_instance_information(problem, use_goal_denotation) for problem, _, _ in parsed_problems]
    # We assume all problems languages are the same and simply pick the first one
    language = parsed_problems[0][0].language
    vocabulary = compute_dl_vocabulary(language)

    # Note that the validator will validate wrt the full sample
    _, model_cache = create_model_cache_from_samples(vocabulary, sample, config.domain, config.parameter_generator, infos)
    validator = AbstractionValidator(model_cache, sample, expanded_state_ids_shuffled)

    k, k_max, k_step = config.initial_concept_bound, config.max_concept_bound, config.concept_bound_step
    assert k <= k_max
    while True:
        print("Working sample idxs: {}".format(sorted(working_sample_idxs)))
        print("Working sample: {}".format(working_sample.info()))
        print("States in Working sample: {}".format(sorted(working_sample.remapping.keys())))
        res, k, abstraction = try_to_compute_abstraction_in_range(config, working_sample, k, k_max, k_step)
        if res == ExitCode.NoAbstractionUnderComplexityBound:
            logging.error("No abstraction possible for given sample set under max. complexity {}".format(k_max))
            return res

        # Otherwise, we learnt an abstraction with max. complexity k. Let's refine the working sample set with it!
        logging.info("Abstraction with k={} found for sample set of size {} ".format(k, working_sample.num_states()))

        # Look for flaws in the full sample set
        flaws = validator.find_flaws(abstraction, config.batch_refinement_size, check_completeness=False)
        if not flaws:
            break
        logging.info("{} flaws found in the computed abstraction: {}".format(len(flaws), sorted(flaws)))

        if set(flaws).intersection(working_sample_idxs):
            msg = "Set of flaws {} intersects with working sample!".format(flaws)
            logging.warning(msg)
            raise RuntimeError(msg)
        # Augment the sample set with the flaws found for the current abstraction
        working_sample_idxs.update(flaws)
        working_sample = sample.resample(working_sample_idxs)  # This will "close" the sample set with children states

    logging.info("The computed abstraction is sound & complete wrt all of the training instances")
    save(config.experiment_dir, dict(abstraction=abstraction))
    return ExitCode.Success


def try_to_compute_abstraction_in_range(config, sample, k_0, k_max, k_step):
    k, abstraction = k_0, None  # To avoid referenced-before-assigned warnings
    for k in range(k_0, k_max + 1, k_step):
        exitcode, abstraction = try_to_compute_abstraction(config, sample, k)
        if exitcode == ExitCode.Success:
            return exitcode, k, abstraction
        else:
            logging.info("No abstraction possible with max. concept complexity {}".format(k))
    return ExitCode.NoAbstractionUnderComplexityBound, k, abstraction


def setup_workspace(config, sample, k):
    """ Create the workspace and update the config for the incremental search with the given parameters """
    # Create a subdir, e.g. inc_k5_s20, to place all the data for this learning attempt
    subdir = os.path.join(config.experiment_dir, "inc_k{}_s{}".format(k, sample.num_states()))
    os.makedirs(subdir, exist_ok=True)
    save(subdir, dict(sample=sample))
    log_sampled_states(sample, os.path.join(subdir, "sample.txt"))
    logging.info("Search for abstraction: |S|={}, k=|{}|".format(sample.num_states(), k))
    logging.info("Data stored in {}".format(subdir))
    subconfig = update_dict(config.to_dict(), max_concept_size=k, experiment_dir=subdir)
    return subconfig


def teardown_workspace(domain, experiment_dir, serialized_feature_filename, clean_workspace, **kwargs):
    """ Clean up the given workspace """
    data = load(experiment_dir, ["abstraction"])
    # if clean_workspace:
    #     print(experiment_dir)
    #     shutil.rmtree(experiment_dir)  # TODO activate
    return data


def load_features(selected_features, domain, filename):
    indexes = set(elem['idx'] for elem in selected_features)
    features = load_selected_features(indexes, domain, filename)
    return process_features(features)


def try_to_compute_abstraction(config, sample, k):
    """ Try to learn a sound&complete abstract action model for the given sample,
    with max. concept complexity given by k
    """
    subconfig = setup_workspace(config, sample, k)

    steps, subconfig = generate_pipeline_from_list([
        CPPFeatureGenerationStep,
        # MaxsatProblemGenerationStep,
        CPPMaxsatProblemGenerationStep,
        MaxsatProblemSolutionStep,
        CPPActionModelStep,
        # InhouseMaxsatSolverStep,  # Blai's polynomial ad-hoc maxsat algorithm
        # ActionModelFromFeatureIndexesStep
    ], **subconfig)

    # Update some parameters and paths of the subconfig and save there transition info
    subconfig["sat_transitions_filename"] = compute_info_filename(subconfig, "sat_transitions.dat")
    subconfig["transitions_filename"] = compute_info_filename(subconfig, "transition-matrix.dat")
    print_transition_matrices(sample, Bunch(subconfig))

    for step in steps:
        exitcode = run_and_check_output(step, "[inc]", SubprocessStepRunner, raise_on_error=False)
        if exitcode != ExitCode.Success:
            return exitcode, None

    # All steps successfully executed, ergo we found an abstraction.
    data = teardown_workspace(**subconfig)
    return ExitCode.Success, data["abstraction"]


class IncrementalLearner:
    """ Generate concepts, compute feature set, compute A_F, all in one """
    def __init__(self, **kwargs):
        self.config = Bunch(self.process_config(kwargs))

    def get_required_attributes(self):
        return ["domain", "experiment_dir", "initial_concept_bound", "max_concept_bound", "concept_bound_step", "initial_sample_size"]

    def get_required_data(self):
        return ["sample"]

    def process_config(self, config):
        check_int_parameter(config, "initial_concept_bound")
        check_int_parameter(config, "concept_bound_step")
        check_int_parameter(config, "max_concept_bound")

        if config["initial_concept_bound"] > config["max_concept_bound"]:
            raise InvalidConfigParameter("initial_config_bound ({}) must be <= than max_concept_bound ({})".
                                         format(config["initial_concept_bound"], config["max_concept_bound"]))
        return config

    def run(self):
        data = load(self.config.experiment_dir, ["sample"])  # Load the initial sample
        rng = np.random.RandomState(self.config.random_seed)
        return full_learning(self.config, data["sample"], rng)


class IncrementalExperiment(Experiment):
    """ This class handcodes the sequence of steps to be run, as the logic behind incremental running is
        more complex.
    """
    def __init__(self, steps, parameters):
        super().__init__([], parameters)  # Simply ignore the steps
        self.parameters = parameters

    def run(self, steps=None):
        console.print_hello()
        # We ignore the specified steps and run the incremental algorithm
        # 1. Extract and resample the whole training set
        initial_steps, config = generate_pipeline_from_list([PlannerStep, TransitionSamplingStep], **self.parameters)
        exitcode = run_and_check_output(initial_steps[0], 0, SubprocessStepRunner)
        exitcode = run_and_check_output(initial_steps[1], 1, SubprocessStepRunner)

        learner = IncrementalLearner(**config)
        exitcode = learner.run()
        if exitcode == ExitCode.Success:
            testing_steps, config = generate_pipeline_from_list([AbstractionTestingStep], **config)
            exitcode = run_and_check_output(testing_steps[0], 999, SubprocessStepRunner)

        return exitcode
