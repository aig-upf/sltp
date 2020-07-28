import os
from sltp.driver import Experiment, generate_pipeline, BENCHMARK_DIR, BASEDIR
from sltp.learn_actions import OptimizationPolicy


def get_experiment_class(kwargs):
    exptype = kwargs.get('experiment_type', None)
    if exptype is not None:  # parameter 'experiment_type' takes precedence
        if exptype == 'incremental':
            from sltp.incremental import IncrementalExperiment
            return IncrementalExperiment
        else:
            return Experiment

    expclass = kwargs.get('experiment_class', None)
    return expclass if expclass is not None else Experiment


def generate_experiment(domain_dir, domain, **kwargs):
    """ """

    if "instances" not in kwargs:
        raise RuntimeError("Please specify domain and instance when generating an experiment")

    instances = kwargs["instances"] if isinstance(kwargs["instances"], (list, tuple)) else [kwargs["instances"]]
    kwargs["domain"] = os.path.join(BENCHMARK_DIR, domain_dir, domain)
    kwargs["instances"] = [os.path.join(BENCHMARK_DIR, domain_dir, i) for i in instances]

    if "test_domain" in kwargs:
        kwargs["test_domain"] = os.path.join(BENCHMARK_DIR, domain_dir, kwargs["test_domain"])
        kwargs["test_instances"] = \
            [os.path.join(BENCHMARK_DIR, domain_dir, i) for i in kwargs.get("test_instances", [])]
        kwargs["test_policy_instances"] = \
            [os.path.join(BENCHMARK_DIR, domain_dir, i) for i in kwargs.get("test_policy_instances", [])]
    else:
        kwargs["test_domain"] = None
        kwargs["test_instances"] = []
        kwargs["test_policy_instances"] = []

    defaults = dict(
        # pipeline="maxsat",
        pipeline="maxsatcpp",
        # pipeline="maxsatcpp_old",
        # pipeline="sat",
        # pipeline="maxsat_poly",

        # The directory where the experiment outputs will be left
        workspace=os.path.join(BASEDIR, 'workspace'),

        # Location of the FS planner, used to do the state space sampling
        # planner_location=os.getenv("FS_PATH", os.path.expanduser("~/projects/code/fs")),
        planner_location=os.getenv("FS_PATH", os.path.expanduser("~/software/github/fs-planner")),

        # Type of sampling procedure. Only breadth-first search implemented ATM
        driver="bfs",

        # Whether to do IW-like pruning of nodes with novelty higher than the specified (default: -1, no pruning at all)
        max_width=-1,

        # Whether we are happy to obtain completeness guarantees only with respect to those states
        # that lie in some arbitrary (one) optimal path. Default: False
        complete_only_wrt_optimal=False,

        # The selection strategy to be used when marking which transitions are considered as optimal.
        # - "arbitrary": marks as optimal the transitions in one single path btw initial state
        #               and (some) goal per instance, chosen arbitrarily.
        # - "complete": marks all transitions that are optimal between some _non-dead_ state and some goal.
        optimal_selection_strategy="arbitrary",

        # Type of sampling. Accepted options are:
        # - "all" (default): Use all expanded states
        # - "random": Use only a random sample of the expanded states, of size given by the option "num_sampled_states"
        # - "optimal": Use those expanded states on some optimal path (only one arbitrary optimal path)
        # Note: ATM random sampling is deactivated
        sampling="all",

        # Number of states to be expanded in the sampling procedure. Either a positive integer, or the string
        # "until_first_goal", or the string "all", both with obvious meanings
        num_states=50,

        # Number randomly sampled states among the set of expanded states. The default of None means
        # all expanded states will be used
        num_sampled_states=None,

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

        # A list of features the user wants to be in the abstraction (possibly along others).
        # Default is [], i.e. just run the generation process without enforcing anything
        enforce_features=[],

        # Max. allowed complexity for distance and conditional features (default: 0)
        distance_feature_max_complexity=0,
        cond_feature_max_complexity=0,

        # Whether to generate comparison features of the type F1 < F2
        comparison_features=False,

        # Method to generate domain parameters (goal or otherwise). If None, goal predicates will
        # be used (default: None)
        parameter_generator=None,

        # Whether to create goal-identifying features (e.g. of the form p_g AND not p_s for every unary predicate
        # apperaring in the goal)
        create_goal_features_automatically=False,

        # Generator of goal-like expressions
        goal_selector=None,

        # Optionally, use a method that gives handcrafted names to the features
        # (default: None, which will use their string representation)
        feature_namer=default_feature_namer,

        # What optimization criteria to use in the max-sat problem
        optimization=OptimizationPolicy.TOTAL_FEATURE_COMPLEXITY,
        # optimization=OptimizationPolicy.NUM_FEATURES

        # Set a random seed for reproducibility (default: 1)
        random_seed=1,

        # the max-sat solver to use. Accepted: openwbo, openwbo-inc, wpm3, maxino
        maxsat_solver='openwbo',
        maxsat_timeout=None,

        # The number of features and actions for the SAT encoding
        encoding_k=10,
        encoding_m=10,

        domain_dir=domain_dir,

        # The Experiment class to be used (e.g. standard, or incremental)
        experiment_class=get_experiment_class(kwargs),

        # The size of the initial batch of states for the incremental sampling approach
        initial_sample_size=100,

        # The max. number of states in the Flaw set when validating an incremental abstraction
        batch_refinement_size=10,

        # Whether to clean the (sub-) workspaces used to compute incremental abstractions after finishing.
        clean_workspace=True,

        # Reduce output to a minimum
        quiet=False,

        # Whether to take into acount states labeled as unsolvable by whatever planner is being used
        compute_unsolvable_states=False,

        # Whether to prune features that have denotation > 0 over *all states*
        prune_positive_features=True,

        # Whether to prune features that have infinity denotation *in some state*
        prune_infty_features=False,

        # Whether ad-hoc solve runs with verbose option
        wsat_solver_verbose=False,

        # Number of states to expand & test on the testing instances
        num_tested_states=50000,

        # Prune those states that are redundante wrt the feature pool before building the CNF theory
        prune_redundant_states=True,

        # Use the D2-Tree for the CNF encoding
        maxsat_encoding="d2tree",

        # Some debugging help to print the denotations of all features over all states (expensive!)
        print_all_denotations=False,

        # By default don't timeout the concept generation process
        concept_generation_timeout=-1,

        # A function to manually provide a transition-classification policy that we want to test
        transition_classification_policy=None,

        # In the transition-separation CNF encoding, whether to distinguish good transitions *only from*
        # unmarked transitions that start in an alive state
        use_only_alive_unmarked_transitions=True,

        # In the transition-separation CNF encoding, whether to distinguish good transitions *only from*
        # unmarked transitions that start in the same state as the good transition
        distinguish_transitions_locally=True,
    )

    parameters = {**defaults, **kwargs}  # Copy defaults, overwrite with user-specified parameters

    steps = generate_pipeline(**parameters)
    exp = parameters["experiment_class"](steps, parameters)
    return exp


def default_feature_namer(s):
    return str(s)
