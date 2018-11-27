import logging
import os
import shutil
from collections import defaultdict

import numpy as np
from tarski.dl import FeatureValueChange

from .sampling import log_sampled_states
from .features import create_denotation_processor
from .returncodes import ExitCode
from .util.misc import update_dict
from .driver import Experiment, generate_pipeline_from_list, PlannerStep, check_int_parameter, \
    InvalidConfigParameter, TransitionSamplingStep, ConceptGenerationStep, SubprocessStepRunner, \
    run_and_check_output, save, FeatureMatrixGenerationStep, MaxsatProblemGenerationStep, MaxsatProblemSolutionStep, \
    ActionModelStep, load, Bunch


def initial_sample_selection(sample, config, rng):
    # Get the initial sample of M states plus one goal state per instance, if any goal is available
    # all_state_ids_shuffled = list(sample.states.keys())
    # rng.shuffle(all_state_ids_shuffled)
    expanded_state_ids_shuffled = list(sample.expanded)
    rng.shuffle(expanded_state_ids_shuffled)
    initial_size = min(len(expanded_state_ids_shuffled), config.initial_sample_size)
    working_sample_idxs = set(expanded_state_ids_shuffled[:initial_size])
    working_sample_idxs.update(sample.get_one_goal_per_instance())
    return working_sample_idxs, expanded_state_ids_shuffled


def full_learning(config, sample, rng):
    working_sample_idxs, expanded_state_ids_shuffled = initial_sample_selection(sample, config, rng)
    working_sample = sample.resample(working_sample_idxs)

    # TODO This needs to be adjusted for sample sets coming from diff instances
    use_goal_denotation = config.parameter_generator is None
    _, processor = create_denotation_processor(config.domain, config.instances[0], use_goal_denotation)
    validator = AbstractionValidator(processor, sample, expanded_state_ids_shuffled)

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
        flaws = validator.find_flaws(abstraction, config.batch_refinement_size)
        if not flaws:
            break
        logging.info("{} flaws found in the computed abstraction: {}".format(len(flaws), sorted(flaws)))

        # Augment the sample set with the flaws found for the current abstraction
        working_sample_idxs.update(flaws)
        working_sample = sample.resample(working_sample_idxs)  # This will "close" the sample set with children states


    logging.info("The computed abstraction is sound & complete wrt all of the training instances")
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


def teardown_workspace(experiment_dir, clean_workspace, **kwargs):
    """ Clean up the given workspace """
    data = load(experiment_dir, ["abstract_actions", "selected_features", "features"])
    if clean_workspace:
        print(experiment_dir)
        # shutil.rmtree(experiment_dir)  # TODO activate
    return data


def try_to_compute_abstraction(config, sample, k):
    """ Try to learn a sound&complete abstract action model for the given sample,
    with max. concept complexity given by k
    """
    subconfig = setup_workspace(config, sample, k)

    steps, subconfig = generate_pipeline_from_list([
        ConceptGenerationStep, FeatureMatrixGenerationStep, MaxsatProblemGenerationStep,
        MaxsatProblemSolutionStep, ActionModelStep], **subconfig)

    for step in steps:
        exitcode = run_and_check_output(step, SubprocessStepRunner, raise_on_error=False)
        if exitcode != ExitCode.Success:
            return exitcode, None

    # All steps successfully executed, ergo we found an abstraction.
    data = teardown_workspace(**subconfig)
    return ExitCode.Success, data


class AbstractionValidator:
    def __init__(self, processor, sample, state_ids):
        """ """
        self.processor = processor
        self.sample = sample
        self.state_ids = state_ids

    def action_captures(self, models, s, sprime, action_effs, feature_set):
        """ Return true iff the abstract action captures (i.e. represents) transition (s, s') """
        for idx, f in feature_set:
            effect_of_action_on_f = action_effs[idx]
            x0 = self.get_possibly_cached_model(models, s).compute_feature_value(f, 0)
            x1 = self.get_possibly_cached_model(models, sprime).compute_feature_value(f, 0)
            diff = f.diff(x0, x1)
            if diff != effect_of_action_on_f:
                return False

        return True

    def action_set_captures(self, models, s, sprime, action_eff_set, feature_set):
        """ Return true iff some action in the given abstract action set captures transition (s, s') """
        return any(self.action_captures(models, s, sprime, effs, feature_set) for a, effs in action_eff_set.items())

    def get_possibly_cached_model(self, models, state_id):
        if state_id not in models:
            state = self.sample.states[state_id]
            models[state_id] = self.processor.compute_model_from_state(state)
        return models[state_id]

    def find_flaws(self, abstraction, max_flaws):
        """ """
        sample = self.sample
        unsound, not_represented = set(), set()
        abstract_actions = abstraction["abstract_actions"]
        pool = abstraction["features"]

        feature_idx = self.compute_feature_idx(abstract_actions)
        models = dict()  # We will store here the model corresponding to each state

        # Map the idxs of selected features to the actual features
        selected_feature_idxs = abstraction["selected_features"]
        selected_features = [(f["idx"], abstraction["features"][f["idx"]]) for f in selected_feature_idxs]
        # assert all(str(f) == f2["name"] for f, f2 in zip(selected_features, selected_feature_idxs))

        logging.info("Abstract action set:\n{}".format(abstract_actions))

        assert all(s in sample.expanded for s in self.state_ids)  # state_ids assumed to contain only expanded states.
        for sid in self.state_ids:
            model = self.get_possibly_cached_model(models, sid)
            # Check soundness
            for action in abstract_actions:
                # We need to cast the actual feature value into a bool to compare it with the abstraction bool value
                is_applicable = all(bool(model.compute_feature_value(pool[idx], 0)) is val for idx, val in action.preconditions)
                if is_applicable and not any(self.action_captures(models, sid, sprime, feature_idx[action], selected_features) for sprime in sample.transitions[sid]):
                    # The abstract action is not sound
                    print("Action {} not *sound* wrt state {}".format(action, sid))
                    unsound.add(sid)

            # Check completeness
            if not all(self.action_set_captures(models, sid, sprime, feature_idx, selected_features) for sprime in sample.transitions[sid]):
                # The abstract action is not complete
                print("Action set not *complete* wrt state {}".format(sid))
                not_represented.add(sid)

            if len(unsound) >= max_flaws:
                break

        return unsound if unsound else not_represented

    def compute_feature_idx(self, actions):
        index = dict()
        for a in actions:
            effects = [(eff.feature, eff.change) for eff in a.effects]
            index[a] = defaultdict(lambda: FeatureValueChange.NIL, effects)
        return index


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
        exitcode = run_and_check_output(initial_steps[0], SubprocessStepRunner)
        exitcode = run_and_check_output(initial_steps[1], SubprocessStepRunner)

        learner = IncrementalLearner(**config)
        exitcode = learner.run()
