import logging
from collections.__init__ import defaultdict, OrderedDict


def sample_generated_states(states, goal_states, transitions, config, rng):
    num_expanded_states = min(len(states), config.num_states)
    if config.num_sampled_states is None:
        return states, goal_states, transitions

    if config.num_sampled_states > num_expanded_states:
        logging.warning("Only {} were expanded, sampling all of them!".format(num_expanded_states))
        return states, goal_states, transitions

    # Although this might fail is some state was a dead-end? In that case we should perhaps revise the code below
    assert num_expanded_states == len(transitions)
    parents = compute_parents(transitions)

    selected = sample_expanded_and_goal_states(config, goal_states, num_expanded_states, parents, rng)
    sampled_expanded = set(selected)

    return remap_sample_expanded_states(sampled_expanded, states, goal_states, transitions)


def sample_expanded_and_goal_states(config, goal_states, num_expanded_states, parents, rng):
    expanded_states = list(range(0, num_expanded_states))
    # We will at least select all of those goal states that have been expanded plus one parent of those that not,
    # so that we maximize the number of goal states in the sample.
    # This is not the only possible strategy, and will be problematic if e.g. most states are goal states, but
    # for the moment being we're happy with it
    enforced = set()
    for x in expanded_states:
        if x in goal_states:
            enforced.add(x)
        elif parents[x]:
            enforced.add(next(iter(parents[x])))  # We simply pick one arbitrary parent of the non-expanded goal state
    #
    enforced = set()
    non_enforced = [i for i in expanded_states if i not in enforced]
    rng.shuffle(non_enforced)
    all_shuffled = list(enforced) + non_enforced
    return all_shuffled[:config.num_sampled_states]


def compute_parents(transitions):
    parents = defaultdict(set)
    for source, targets in transitions.items():
        for t in targets:
            parents[t].add(source)
    return parents


def remap_sample_expanded_states(sampled_expanded, states, goal_states, transitions):
    all_in_sample = set(sampled_expanded)
    for s in sampled_expanded:
        all_in_sample.update(transitions[s])

    # Now all_in_sample contains all states we want to sample plus their children state IDs

    ordered_sample = list(sorted(all_in_sample))
    idx = {x: i for i, x in enumerate(ordered_sample)}
    revidx = {}

    # Pick the selected elements from the data structures
    new_goal_states = {idx[x] for x in ordered_sample if x in goal_states}
    new_states = OrderedDict()
    for i, s in states.items():
        assert i == s[0]
        if i in idx:
            new_states[idx[i]] = (idx[i], s[1])

    new_transitions = defaultdict(set)
    for source, targets in transitions.items():
        if source in sampled_expanded:
            new_transitions[idx[source]] = {idx[t] for t in targets}

    return new_states, new_goal_states, new_transitions
