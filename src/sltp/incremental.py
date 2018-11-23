import logging
import os

from sltp.util.misc import update_dict, project_dict
from .driver import Experiment, generate_pipeline_from_list, PlannerStep, Step, check_int_parameter, \
    InvalidConfigParameter, TransitionSamplingStep, ConceptGenerationStep, StepRunner, SubprocessStepRunner, \
    _run_and_check_output


class FullLearningStep(Step):
    """ Generate concepts, compute feature set, compute A_F, all in one """
    def __init__(self, **kwargs):
        super().__init__(**kwargs)

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

        config["concept_dir"] = os.path.join(config["experiment_dir"], 'terms')
        config["concept_generator"] = config.get("concept_generator", None)
        config["feature_generator"] = config.get("feature_generator", None)
        config["enforce_features"] = config.get("enforce_features", None)
        config["parameter_generator"] = config.get("parameter_generator", None)
        config["distance_feature_max_complexity"] = config.get("distance_feature_max_complexity", 0)
        config["max_concept_grammar_iterations"] = config.get("max_concept_grammar_iterations", None)

        return config

    def description(self):
        return "Full Learning Step"

    def get_step_runner(self):
        return full_learning


def full_learning(config, data, rng):
    # TransitionSamplingStep

    all_states = data.sample.states
    all_transitions = data.sample.transitions

    # Get the initial sample of M states
    working_set = project_dict(all_states, rng.choice(list(all_states.keys()), config.initial_sample_size, replace=False))

    k = 1
    while True:
        parameters = update_dict(config.to_dict(), states=working_set, max_concept_size=k)
        conceptgen = ConceptGenerationStep(**parameters)
        _run_and_check_output(conceptgen, SubprocessStepRunner)

        if "UNSAT":
            pass


        x = 1




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
        _run_and_check_output(initial_steps[0], SubprocessStepRunner)
        _run_and_check_output(initial_steps[1], SubprocessStepRunner)

        learner = FullLearningStep(**config)
        result = _run_and_check_output(learner, StepRunner)

        print("YES")

        # ConceptGenerationStep,
        # FeatureMatrixGenerationStep,
        # MaxsatProblemGenerationStep,
        # MaxsatProblemSolutionStep,
        # ActionModelStep,


