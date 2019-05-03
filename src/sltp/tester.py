import logging
import sys

from .models import FeatureModel
from .util.tools import Abstraction, AbstractAction, generate_qualitative_effects_from_transition
from .features import generate_model_cache, create_model_factory, compute_static_atoms
from .returncodes import ExitCode
from .sampling import read_transitions_from_files
from .validator import AbstractionValidator
from . import PYPERPLAN_DIR


def run(config, data, rng):
    if config.test_domain is None:
        logging.info("No testing instances were specified")
        return ExitCode.Success, dict()
    logging.info("Testing learnt abstraction on instances: {}".format(config.test_instances))

    # First test whether the abstraction is sound on the test instances
    assert isinstance(data.abstraction, Abstraction)
    sample, _ = read_transitions_from_files(config.test_sample_files)
    model_cache = generate_model_cache(config.test_domain, config.test_instances, sample, config.parameter_generator)

    validator = AbstractionValidator(model_cache, sample, list(sample.expanded))
    flaws = validator.find_flaws(data.abstraction, 1, check_completeness=False)
    if flaws:
        logging.error("The computed abstraction is not sound".format())
        return ExitCode.AbstractionFailsOnTestInstances, dict()

    logging.info("The computed abstraction is sound & complete in all test instances!".format())

    pyperplan = _import_pyperplan()
    for instance in config.test_instances:
        logging.info('Testing computed abstract policy heuristic on instance "{}"'.format(instance))

        try:
            _run_pyperplan(pyperplan, config.test_domain, instance, data.abstraction, config.parameter_generator)
        except PolicySearchException as e:
            logging.warning("Testing of abstract policy failed with code: {}".format(e.code))
            return e.code, dict()

    return ExitCode.Success, dict()


def _import_pyperplan():
    sys.path.insert(0, PYPERPLAN_DIR)
    import pyperplan
    sys.path = sys.path[1:]
    return pyperplan


def _run_pyperplan(pyperplan, domain, instance, abstraction, parameter_generator):
    # Let's fake a pyperplan command-line arguments object so that
    # we don't have to worry about default parameters. We use "gbf" but will override that next.
    args = pyperplan.parse_args([domain, instance])

    # Parse the domain & instance and create a model generator
    problem, model_factory = create_model_factory(domain, instance, parameter_generator)
    static_atoms, _ = compute_static_atoms(problem)

    # Compute an actual policy object that returns the next action to be applied
    search_policy = create_abstract_policy(model_factory, static_atoms, abstraction)

    # And now we inject our desired search and heuristic functions
    args.forced_search = create_pyperplan_abstract_policy_based_search(pyperplan, search_policy)

    # And run pyperplan!
    pyperplan.main(args)


def translate_atom(atom):
    assert atom[0] == '(' and atom[-1] == ')'
    return atom[1:-1].split(' ')


def translate_state(state, static_atoms):
    """ Translate a pyperplan-like state into a list with the format required by SLTP's concept denotation processor """
    return [translate_atom(atom) for atom in state] + list(static_atoms)


def generate_model_from_state(model_factory, state, static_atoms):
    translated = translate_state(state, static_atoms)
    return FeatureModel(model_factory.create_model(translated))


def create_abstract_policy(model_factory, static_atoms, abstraction):
    assert isinstance(abstraction, Abstraction)

    def _policy(state, successor_generator):
        m0 = generate_model_from_state(model_factory, state, static_atoms)

        # Obtain the abstract action dictated by the policy
        abs_state, abs_action = abstraction.optimal_action(m0)
        if abs_action is None:
            raise PolicySearchException(ExitCode.AbstractPolicyNotCompleteOnTestInstances)

        assert isinstance(abs_action, AbstractAction)
        # Find a concrete representative for that abstract action - TODO Could be more efficient than a linear search
        for op, succ in successor_generator:
            m1 = generate_model_from_state(model_factory, succ, static_atoms)
            concrete_changeset = generate_qualitative_effects_from_transition(abstraction.features, m0, m1)

            if concrete_changeset == abs_action.effects:
                # We stick with the 1st concrete operator representing the abstract action dictated by the policy
                return op, succ

        logging.warning("No concrete action represents abstract action {} on concrete state {}".format(abs_action, state))
        raise PolicySearchException(ExitCode.AbstractPolicyNotSoundOnTestInstances)

    return _policy


def create_pyperplan_abstract_policy_based_search(pyperplan, search_policy):
    """ Creates a pyperplan like search object that uses the abstract generalized policy to perform the search """
    searchspace = pyperplan.search.searchspace

    def policy_based_search(task):  # The actual search object
        steps = 0  # Some bookkeping
        node = searchspace.make_root_node(task.initial_state)
        closed = set()

        while not task.goal_reached(node.state):
            closed.add(node.state)
            steps += 1
            if steps % 1000 == 0:
                logging.debug("Number of visited states so far in policy-based search: {}".format(steps))

            successor_generator = ((op, op.apply(node.state)) for op in task.operators if op.applicable(node.state))
            op, succ = search_policy(node.state, successor_generator)

            if succ in closed:  # loop detection
                logging.error("Policy incurred in a loop after {} steps. Repeated node: {}".format(steps, node.state))
                raise PolicySearchException(ExitCode.AbstractPolicyNonTerminatingOnTestInstances)

            node = searchspace.make_child_node(node, op, succ)

        logging.info("Goal found!")
        return node.extract_solution()

    return policy_based_search


class PolicySearchException(Exception):
    def __init__(self, code):
        super().__init__("")
        self.code = code
