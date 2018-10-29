import itertools
import json
import logging

from collections import defaultdict, OrderedDict

from .util.command import read_file
from .util.naming import filename_core


def mark_optimal(goal_states, root_states, parents):
    """ Collect all those states that lie on one arbitrary optimal path to the goal """
    roots = set(root_states)
    optimal_transitions = set()
    for goal in goal_states:
        previous = current = goal

        while current not in roots:
            # A small trick: node IDs are ordered by depth, so we can pick the min parent ID and know the resulting path
            # will be optimal
            current = min(parents[current])
            optimal_transitions.add((current, previous))  # We're traversing backwards
            previous = current

    return optimal_transitions


def print_atom(atom):
    assert len(atom) > 0
    if len(atom) == 1:
        return atom[0]
    return "{}({})".format(atom[0], ", ".join(atom[1:]))


def log_sampled_states(states, goals, transitions, expanded, optimal_transitions, root_states, unsolvable, resampled_states_filename):
    # We need to recompute the parenthood relation with the remapped states to log it!
    parents = compute_parents(transitions)
    optimal_s = set(x for x, _ in optimal_transitions)

    with open(resampled_states_filename, 'w') as f:
        for id_, state in states.values():
            state_parents = ", ".join(sorted(map(str, parents[id_])))
            state_children = ", ".join(sorted(map(str, transitions[id_])))
            atoms = ", ".join(print_atom(atom) for atom in state)
            is_goal = "*" if id_ in goals else ""
            is_expanded = "^" if expanded(id_) else ""
            is_root = "=" if id_ in root_states else ""
            is_optimal = "+" if id_ in optimal_s else ""
            is_unsolvable = "U" if id_ in unsolvable else ""
            print("#{}{}{}{}{}{} (parents: {}, children: {}):\n\t{}".
                  format(id_, is_root, is_goal, is_optimal, is_expanded, is_unsolvable, state_parents, state_children, atoms), file=f)
    logging.info('Resampled states logged at "{}"'.format(resampled_states_filename))


def sample_first_x_states(root_states, sample_sizes):
    sampled = set()
    assert len(root_states) == len(sample_sizes)
    for root, size in zip(root_states, sample_sizes):
        sampled.update(range(root, root+size))
    return sampled


def sample_generated_states(config, rng):
    logging.info('Loading state space samples...')
    states, goal_states, transitions, root_states, goals_by_instance, unsolvable = read_transitions(config.sample_files)
    if not goal_states:
        raise RuntimeError("No goal found in the sample - increase # expanded states!")

    # if config.num_sampled_states is None and config.max_width < 1 and not config.complete_only_wrt_optimal:
    #     return states, goal_states, transitions

    parents = compute_parents(transitions)

    optimal, optimal_transitions = set(), set()

    if config.sampling == "random":
        assert config.num_sampled_states is not None
        assert not config.complete_only_wrt_optimal
        selected = random_sample(config, goal_states, rng, states, transitions, parents)

    else:
        assert config.sampling in ("all", "optimal")

        # For each instance, we keep the first-reached goal, as a way of selecting an arbitrary optimal path.
        goal_selection = set(min(x) for x in goals_by_instance)

        optimal_transitions = mark_optimal(goal_selection, root_states, parents)
        states_in_some_optimal_transition = set(itertools.chain.from_iterable(optimal_transitions))
        logging.info("Resampling: states/goals/optimal tx: {} / {} / {}".format(len(states), len(goal_states), len(optimal_transitions)))

        if config.sampling == "optimal":  # Select only the optimal states
            assert False, "not sure this makes too much sense!"
            # selected = set(optimal)

        else:  # sampling = "all"
            if config.num_sampled_states is not None:
                selected = sample_first_x_states(root_states, config.num_sampled_states)
                selected.update(states_in_some_optimal_transition)
            else:
                selected = set(states)

    states, goals, transitions, optimal_transitions, root_states, unsolvable = \
        remap_sample_expanded_states(set(selected), states, goal_states, transitions, optimal_transitions, set(root_states), unsolvable)

    expanded = lambda s: len(transitions[s]) > 0
    total_tx = sum(len(transitions[sid]) for sid, _ in states.values())
    logging.info("Resampling (states / goals / unsolvable / tx / optimal tx): {} / {} / {} / {}".format(len(states), len(goals), len(unsolvable), total_tx, len(optimal_transitions)))
    log_sampled_states(states, goals, transitions, expanded, optimal_transitions, root_states, unsolvable, config.resampled_states_filename)
    return states, goals, transitions, optimal_transitions, root_states, unsolvable


def random_sample(config, goal_states, rng, states, transitions, parents):
    num_states = min(len(states), config.num_states)
    assert len(set(config.num_sampled_states)) == 1, "ATM only random sample with fixed sample size"
    sample_size = config.num_sampled_states[0]
    if sample_size > num_states:
        raise RuntimeError(
            "Only {} expanded statesm cannot sample {}".format(num_states, sample_size))
    # Although this might fail is some state was a dead-end? In that case we should perhaps revise the code below
    assert num_states == len(transitions)

    selected = sample_expanded_and_goal_states(config.num_sampled_states, goal_states, num_states, parents, rng)
    return selected


def sample_expanded_and_goal_states(sample_size, goal_states, num_states, parents, rng):
    expanded_states = list(range(0, num_states))
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
    return all_shuffled[:sample_size]


def compute_parents(transitions):
    parents = defaultdict(set)
    for source, targets in transitions.items():
        for t in targets:
            parents[t].add(source)
    return parents


def remap_sample_expanded_states(sampled_expanded, states, goal_states, transitions, optimal_transitions, root_states, unsolvable):
    # all_in_sample will contain all states we want to sample plus their children state IDs
    all_in_sample = set(sampled_expanded)
    for s in sampled_expanded:
        all_in_sample.update(transitions[s])

    ordered_sample = list(sorted(all_in_sample))
    idx = {x: i for i, x in enumerate(ordered_sample)}

    # Pick the selected elements from the data structures
    new_goal_states = {idx[x] for x in ordered_sample if x in goal_states}
    new_unsolvable = {idx[x] for x in ordered_sample if x in unsolvable}
    new_root = {idx[x] for x in ordered_sample if x in root_states}
    new_optimal_transitions = {(idx[x], idx[y]) for x, y in optimal_transitions}

    new_states = OrderedDict()
    for i, s in states.items():
        assert i == s[0]
        if i in idx:
            new_states[idx[i]] = (idx[i], s[1])

    new_transitions = defaultdict(set)
    for source, targets in transitions.items():
        if source in sampled_expanded:
            new_transitions[idx[source]] = {idx[t] for t in targets}

    return new_states, new_goal_states, new_transitions, new_optimal_transitions, new_root, new_unsolvable


def normalize_atom_name(name):
    tmp = name.replace('()', '').replace(')', '').replace('(', ',')
    if "=" in tmp:  # We have a functional atom
        tmp = "=," + tmp.replace("=", ',')  # Mark the string putting the "=" as the first position

    return tmp.split(',')


def remap_state_ids(states, goals, transitions, unsolvable, remap):

    new_goals = {remap(x) for x in goals}
    new_unsolvable = {remap(x) for x in unsolvable}
    new_states = OrderedDict()
    for i, s in states.items():
        assert i == s[0]
        new_states[remap(i)] = (remap(i), s[1])

    new_transitions = defaultdict(set)
    for source, targets in transitions.items():
        new_transitions[remap(source)] = {remap(t) for t in targets}

    return new_states, new_goals, new_transitions, new_unsolvable


def read_transitions(filenames):
    all_samples = [read_single_sample_file(f) for f in filenames]
    assert len(all_samples) > 0
    states, goals, transitions, unsolvable = all_samples[0]
    assert states[0][0] == 0
    root_states = [0]
    goals_by_instance = [set(goals)]
    for s, g, tx, unsolv in all_samples[1:]:
        starting_state_id = len(states)
        s, g, tx, unsolv = remap_state_ids(s, g, tx, unsolv, remap=lambda state: state + starting_state_id)
        assert next(iter(s)) == starting_state_id
        root_states.append(starting_state_id)
        states.update(s)
        transitions.update(tx)
        goals.update(g)
        goals_by_instance.append(g)
        unsolvable.update(unsolv)

    return states, goals, transitions, root_states, goals_by_instance, unsolvable


def read_single_sample_file(filename):
    states_by_id = {}
    states_by_str = {}
    transitions = defaultdict(set)
    transitions_inv = defaultdict(set)
    seen = set()
    goal_states = set()
    unsolvable_states  = set()

    def register_transition(state):
        transitions[state['parent']].add(state['id'])
        transitions_inv[state['id']].add(state['parent'])

    def register_state(state):
        data = (state['id'], state['normalized_atoms'])
        states_by_str[state['atoms_string']] = data
        states_by_id[state['id']] = data
        seen.add(state['id'])
        if j['goal']:
            goal_states.add(state['id'])
        if 'unsolvable' in j and j['unsolvable']:
            unsolvable_states.add(state['id'])

    raw_file = [line.replace(' ', '') for line in read_file(filename) if line[0:6] == '{"id":']
    for raw_line in raw_file:
        j = json.loads(raw_line)
        j['normalized_atoms'] = [normalize_atom_name(atom) for atom in j['atoms']]
        j['atoms_string'] = str(j['normalized_atoms'])

        if j['id'] in seen:
            # We hit a repeated state in the search, so we simply need to record the transition
            register_transition(j)
            continue

        # The state must be _really_ new
        assert j['atoms_string'] not in states_by_str
        register_state(j)

        if j['parent'] != j['id']:  # i.e. if not in the root node, which has 0 as its "fake" parent
            # BELOW CHECK NO LONGER CORRECT, AS REPEATED NODES ARE ALSO RECORDED WHEN ENCOUNTERED IN DIFF TRANSITIONS
            # assert json.loads(raw_file[j['parent']])['id'] == j['parent']  # Just a check
            register_transition(j)

    # check soundness
    for src in transitions:
        assert src in states_by_id
        for dst in transitions[src]:
            assert dst in states_by_id

    assert sum([len(t) for t in transitions.values()]) == sum([len(t) for t in transitions_inv.values()])

    logging.info('%s: #lines-raw-file=%d, #state-by-str=%d, #states-by-id=%d, #transition-entries=%d, #transitions=%d' %
                 (filename_core(filename), len(raw_file), len(states_by_str), len(states_by_id), len(transitions),
                  sum([len(targets) for targets in transitions.values()])))

    ordered = OrderedDict()  # Make sure we return an ordered dictionary
    for id_ in sorted(states_by_id.keys()):
        ordered[id_] = states_by_id[id_]
    return ordered, goal_states, transitions, unsolvable_states
