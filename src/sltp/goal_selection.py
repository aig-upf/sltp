import itertools
import logging

from tarski.dl import Concept

from .features import generate_model_cache
from .util.serialization import parse_concept


def select_goal_maximizing_states(config, sample):
    selector = config.goal_selector
    if selector is None:
        return []

    language, model_cache = generate_model_cache(config.domain, config.instances, sample, config.parameter_generator)

    goal_expression = selector(language)
    assert isinstance(goal_expression, str)
    concept = parse_concept(language, goal_expression)
    logging.info("Searching for states maximizing the extension of goal-like expression: {}".format(concept))
    assert isinstance(concept, Concept)

    root_extensions = []
    for sid in sample.instance_roots:
        root_extensions.append(concept.denotation(model_cache.get_term_model(sid)))

    maxargs = [[]] * len(sample.roots)
    maxvals = [-1] * len(sample.roots)
    for sid in sample.states.keys():
        instance = sample.instance[sid]
        root_extension = root_extensions[instance]
        extension = concept.denotation(model_cache.get_term_model(sid))
        # TODO This could be improved to take into account the precision / recall tradeoff?
        diff = (extension & ~root_extension).count()
        if diff > maxvals[instance]:
            maxargs[instance] = [sid]
            maxvals[instance] = diff
        elif diff == maxvals[instance]:
            maxargs[instance].append(sid)

    # logging.info("Goal-maximizing state: {}".format(sample.states[maxargs[0]]))

    # TODO We could just sample a max. number of these if too many?
    goal_maximizing_states = set(itertools.chain.from_iterable(l[:5] for l in maxargs))
    logging.info("{} states marked as goal-maximizing".format(len(goal_maximizing_states)))
    return goal_maximizing_states

