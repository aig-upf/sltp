import logging
import os

from .driver import Experiment, generate_pipeline_from_list, PlannerStep, Step, check_int_parameter, \
    InvalidConfigParameter


class FullLearningStep(Step):
    """ Generate concepts, compute feature set, compute A_F, all in one """
    def __init__(self, **kwargs):
        super().__init__(**kwargs)

    def get_required_attributes(self):
        return ["domain", "experiment_dir", "initial_concept_bound", "max_concept_bound", "concept_bound_step"]

    def get_required_data(self):
        return ["states"]

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
        return "Full Learning Step!"

    def get_step_runner(self):
        return full_learning


def full_learning(config, data, rng):
    # TransitionSamplingStep
    sample_set = config.states



def run_step(step):
    result = step.run()
    if result is not None:
        msg = 'Critical error while processing step "{}". Execution will be terminated.'.format(step.description())
        msg += 'Error message:\t{}'.format(result)
        logging.error(msg)
        raise RuntimeError(msg)
    return result


class IncrementalExperiment(Experiment):
    def __init__(self, steps, parameters):
        super().__init__([], parameters)  # Simply ignore the steps
        self.parameters = parameters

    def run(self, args=None):
        self.hello(args)
        # We ignore whatever specified in the command line and run the incremental algorithm

        # 1. Sample the state spaces
        # ppl = generate_pipeline_from_list(inc_pipeline, **self.parameters)
        # result = ppl[0].run()

        planner = PlannerStep(**self.parameters)
        run_step(planner)

        learner = FullLearningStep(**planner.config)
        run_step(learner)


