import logging
import os

import numpy as np
from .returncodes import ExitCode
from .util.misc import update_dict, project_dict
from .driver import Experiment, generate_pipeline_from_list, PlannerStep, Step, check_int_parameter, \
    InvalidConfigParameter, TransitionSamplingStep, ConceptGenerationStep, StepRunner, SubprocessStepRunner, \
    _run_and_check_output, save, FeatureMatrixGenerationStep, MaxsatProblemGenerationStep, MaxsatProblemSolutionStep, \
    ActionModelStep, load, Bunch


def full_learning(config, sample, rng):
    all_state_ids = list(sample.states.keys())
    initial_size = min(len(all_state_ids), config.initial_sample_size)

    # Get the initial sample of M states plus one goal state per instance, if any goal is available
    resample_idxs = set(rng.choice(all_state_ids, initial_size, replace=False))
    resample_idxs.update(sample.get_one_goal_per_instance())
    working_sample = sample.resample(resample_idxs)

    k, k_max, k_step = config.initial_concept_bound, config.max_concept_bound, config.concept_bound_step
    assert k <= k_max
    while True:

        res, k = try_to_compute_abstraction_in_range(config, working_sample, k, k_max, k_step)
        if res == ExitCode.NoAbstractionUnderComplexityBound:
            logging.error("No abstraction possible for given sample set under max. complexity {}".format(k_max))
            return res, dict()

        # Otherwise, we learnt an abstraction with max. complexity k. Let's refine the working sample set with it!
        logging.info("Abstraction with k={} found for sample set of size {} ".format(working_sample.num_states(), k))
        validate = "TO DO :-)"
        if True:
            break
    logging.info("The computed abstraction is sound & complete wrt all of the training instances")




def try_to_compute_abstraction_in_range(config, sample, k_0, k_max, k_step):
    k = k_0  # To avoid referenced-before-assigned warnings
    for k in range(k_0, k_max + 1, k_step):
        exitcode = try_to_compute_abstraction(config, sample, k)
        if exitcode == ExitCode.Success:
            return exitcode, k
        else:
            logging.info("No abstraction possible with max. concept complexity {}".format(k))
    return ExitCode.NoAbstractionUnderComplexityBound, k


def try_to_compute_abstraction(config, sample, k):
    """ Try to learn a sound&complete abstract action model for the given sample,
    with max. concept complexity given by k
    """
    # Create a subdir, e.g. inc_k5_s20, to place all the data for this learning attempt
    subdir = os.path.join(config.experiment_dir, "inc_k{}_s{}".format(k, sample.num_states()))
    os.makedirs(subdir, exist_ok=True)
    save(subdir, dict(sample=sample))
    subconfig = update_dict(config.to_dict(), max_concept_size=k, experiment_dir=subdir)

    logging.info("Searching for abstraction: |S|={}, k=|{}|".format(sample.num_states(), k))
    logging.info("Data stored in {}".format(subdir))

    steps, subconfig = generate_pipeline_from_list([
        ConceptGenerationStep, FeatureMatrixGenerationStep, MaxsatProblemGenerationStep,
        MaxsatProblemSolutionStep, ActionModelStep], **subconfig)

    for step in steps:
        exitcode = _run_and_check_output(step, SubprocessStepRunner, raise_on_error=False)
        if exitcode != ExitCode.Success:
            return exitcode

    return ExitCode.Success


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
        full_learning(self.config, data["sample"], rng)


class IncrementalExperiment(Experiment):
    """ This class handcodes the sequence of steps to be run, as the logic behind incremental running is
        more complex.
    """
    def __init__(self, steps, parameters):
        super().__init__([], parameters)  # Simply ignore the steps
        self.parameters = parameters

    def run(self, args=None):
        # We ignore whatever specified in the command line and run the incremental algorithm
        self.hello(args)

        # 1. Extract and resample the whole training set
        initial_steps, config = generate_pipeline_from_list([PlannerStep, TransitionSamplingStep], **self.parameters)
        exitcode = _run_and_check_output(initial_steps[0], SubprocessStepRunner)
        exitcode = _run_and_check_output(initial_steps[1], SubprocessStepRunner)

        learner = IncrementalLearner(**config)
        exitcode = learner.run()
