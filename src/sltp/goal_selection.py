import itertools
import logging

from tarski.dl import Concept, Role

from .features import generate_model_cache
from .util.serialization import parse_concept


def procedural_selection(language, selector, sample, rng):
    logging.info("Searching for states maximizing the extension of user-provided procedural selector")

    maxargs = [[]] * len(sample.roots)
    maxvals = [-1] * len(sample.roots)

    for sid, state in sample.states.items():
        instance = sample.instance[sid]
        value = selector(language, state)
        if value > maxvals[instance]:
            maxargs[instance] = [sid]
            maxvals[instance] = value
        elif value == maxvals[instance]:
            maxargs[instance].append(sid)

    # TODO We could just sample a max. number of these if too many?
    # maxargs is sorted by increasing state ID, meaning that the states further away from the initial state
    # are the last ones in the list. These are the ones we want to pick
    goals_per_instance = 8
    goal_maximizing_states = set()
    for instgoals in maxargs:
        k = max(len(instgoals), goals_per_instance*10)
        goal_maximizing_states.update(rng.choice(instgoals[-k:], goals_per_instance, replace=False))
    # goal_maximizing_states = set(itertools.chain.from_iterable(l[-goals_per_instance:] for l in maxargs))
    logging.info("{} states marked as goal-maximizing, max. value: {}".format(len(goal_maximizing_states), maxvals))
    logging.debug("Goal-maximizing states: {}".format(goal_maximizing_states))
    assert goal_maximizing_states
    return goal_maximizing_states


def select_goal_maximizing_states(config, sample, rng):
    selector = config.goal_selector
    if selector is None:
        return []

    language, model_cache = generate_model_cache(config.domain, config.instances, sample, config.parameter_generator)

    if hasattr(selector, "procedural"):
        return procedural_selection(language, selector, sample, rng)

    goal_expression = selector(language)
    assert isinstance(goal_expression, str)
    concept = parse_concept(language, goal_expression)
    logging.info("Searching for states maximizing the extension of goal-like expression: {}".format(concept))
    assert isinstance(concept, (Concept, Role))

    root_extensions = []
    for sid in sample.instance_roots:
        root_extensions.append(concept.denotation(model_cache.get_term_model(sid)))

    maxargs = [[]] * len(sample.roots)
    # maxvals = [-1] * len(sample.roots)
    maxvals = [float("inf")] * len(sample.roots)
    for sid in sample.states.keys():
        instance = sample.instance[sid]
        root_extension = root_extensions[instance]
        extension = concept.denotation(model_cache.get_term_model(sid))
        # TODO This could be improved to take into account the precision / recall tradeoff?
        diff = (extension & ~root_extension).count()
        if diff < maxvals[instance]:
            maxargs[instance] = [sid]
            maxvals[instance] = diff
        elif diff == maxvals[instance]:
            maxargs[instance].append(sid)

    # TODO We could just sample a max. number of these if too many?
    goal_maximizing_states = set(itertools.chain.from_iterable(l[:5] for l in maxargs))
    logging.info("{} states marked as goal-maximizing".format(len(goal_maximizing_states)))
    logging.debug("Goal-maximizing states: {}".format(goal_maximizing_states))
    assert goal_maximizing_states
    return goal_maximizing_states

