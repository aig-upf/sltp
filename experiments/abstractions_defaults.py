import os
import sys


def generate_experiment(domain_dir, domain, **kwargs):
    """ """

    sys.path.insert(0, '..')
    from driver import Experiment, generate_pipeline, BENCHMARK_DIR
    from learn_actions import OptimizationPolicy

    if "instance" not in kwargs:
        raise RuntimeError("Please specify domain and instance when generating an experiment")

    kwargs["domain"] = os.path.join(BENCHMARK_DIR, domain_dir, domain)
    kwargs["instance"] = os.path.join(BENCHMARK_DIR, domain_dir, kwargs["instance"])

    defaults = dict(
        pipeline="maxsat",
        # pipeline="sat",
        # pipeline="heuristic",

        # Location of the FS planner, used to do the state space sampling
        planner_location=os.getenv("FS_PATH", os.path.expanduser("~/projects/code/fs")),

        # Type of sampling procedure. Only breadth-first search implemented ATM
        driver="bfs",

        # Number of states to be expanded in the sampling procedure
        num_states=50,

        # Max. size of the generated concepts (mandatory)
        max_concept_size=10,

        # Max. number of iterations of the concept-generation grammar. Optional. Defaults to infinity,
        # in which case the limit is set by max_concept_size alone.
        max_concept_grammar_iterations=2,

        # Provide a special, handcrafted method to generate concepts, if desired.
        # This will override the standard concept generation procedure (default: None)
        concept_generator=None,

        # Or, alternatively, provide directly the features instead of the concepts (default: None)
        feature_generator=None,

        # Max. allowed complexity for distance features (default: 0)
        distance_feature_max_complexity=0,

        # Method to generate domain parameters (goal or otherwise). If None, goal predicates will
        # be used (default: None)
        parameter_generator=None,

        # Use the relaxed (weak) increase semantics
        relax_numeric_increase=False,

        # Optionally, use a method that gives handcrafted names to the features
        # (default: None, which will use their string representation)
        feature_namer=None,

        # What optimization criteria to use in the max-sat problem
        optimization=OptimizationPolicy.TOTAL_FEATURE_COMPLEXITY,
        # optimization=OptimizationPolicy.NUM_FEATURES

        # The number of features and actions for the SAT encoding
        encoding_k=10,
        encoding_m=10,
    )

    parameters = {**defaults, **kwargs}  # Copy defaults, overwrite with user-specified parameters

    steps = generate_pipeline(**parameters)
    exp = Experiment(steps)

    return exp
