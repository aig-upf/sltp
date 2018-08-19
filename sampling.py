
from collections.__init__ import defaultdict, OrderedDict


def iw_sampling(goal_states, states, transitions):
    assert len(goal_states) == 1
    goal = next(iter(goal_states))
    parents = compute_parents(transitions)

    selected = {goal}
    current = goal
    while current != 0:
        # A small trick - node IDs are ordered by depth, so we can pick the min parent ID and know the resulting path
        # will be optimal
        current = min(parents[current])
        selected.add(current)
    return selected


def sample_generated_states(states, goal_states, transitions, config, rng):

    if config.num_sampled_states is None and config.max_width < 1:
        return states, goal_states, transitions

    if config.num_sampled_states is not None:
        selected = random_sample(config, goal_states, rng, states, transitions)
    else:
        assert config.max_width >= 1
        selected = iw_sampling(goal_states, states, transitions)

    return remap_sample_expanded_states(set(selected), states, goal_states, transitions)


def random_sample(config, goal_states, rng, states, transitions):
    num_expanded_states = min(len(states), config.num_states)
    if config.num_sampled_states > num_expanded_states:
        raise RuntimeError(
            "Only {} expanded statesm cannot sample {}".format(num_expanded_states, config.num_sampled_states))
    # Although this might fail is some state was a dead-end? In that case we should perhaps revise the code below
    assert num_expanded_states == len(transitions)
    parents = compute_parents(transitions)
    selected = sample_expanded_and_goal_states(config, goal_states, num_expanded_states, parents, rng)
    return selected


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
