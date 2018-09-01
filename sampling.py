import json
import logging

from collections import defaultdict, OrderedDict

from util.command import read_file
from util.naming import filename_core


def mark_optimal(goal_states, root_states, parents):
    """ Collect all those states that lie on one arbitrary optimal path to the goal """
    optimal = set()
    for goal in goal_states:
        optimal.add(goal)
        current = goal
        while current not in root_states:
            # A small trick: node IDs are ordered by depth, so we can pick the min parent ID and know the resulting path
            # will be optimal
            current = min(parents[current])
            optimal.add(current)
    return optimal


def print_atom(atom):
    assert len(atom) > 0
    if len(atom) == 1:
        return atom[0]
    return "{}({})".format(atom[0], ", ".join(atom[1:]))


def log_sampled_states(states, goals, transitions, expanded, optimal, resampled_states_filename):
    # We need to recompute the parenthood relation with the remapped states to log it!
    parents = compute_parents(transitions)

    with open(resampled_states_filename, 'w') as f:
        for id_, state in states.values():
            state_parents = ", ".join(sorted(map(str, parents[id_])))
            state_children = ", ".join(sorted(map(str, transitions[id_])))
            atoms = ", ".join(print_atom(atom) for atom in state)
            is_goal = "*" if id_ in goals else ""
            is_expanded = "^" if expanded(id_) else ""
            is_optimal = "+" if id_ in optimal and not is_goal else ""
            print("#{}{}{}{} (parents: {}, children: {}):\n\t{}".
                  format(id_, is_goal, is_optimal, is_expanded, state_parents, state_children, atoms), file=f)
    logging.info('Resampled states logged at "{}"'.format(resampled_states_filename))


def sample_first_x_states(root_states, num_states):
    sampled = set()
    for root in root_states:
        sampled.update(range(root, root+num_states))
    return sampled


def sample_generated_states(config, rng):
    logging.info('Loading state space samples...')
    states, goal_states, transitions, root_states = read_transitions(config.sample_files)
    if not goal_states:
        raise RuntimeError("No goal found in the sample - increase # expanded states!")

    # if config.num_sampled_states is None and config.max_width < 1 and not config.complete_only_wrt_optimal:
    #     return states, goal_states, transitions

    parents = compute_parents(transitions)

    if config.sampling == "random":
        assert config.num_sampled_states is not None
        assert not config.complete_only_wrt_optimal
        selected = random_sample(config, goal_states, rng, states, transitions, parents)
        optimal = set()

    else:
        assert config.sampling in ("all", "optimal")

        optimal = mark_optimal(goal_states, root_states, parents)

        if config.sampling == "optimal":  # Select only the optimal states
            selected = set(optimal)

        else:  # sampling = "all"
            if config.num_sampled_states is not None:
                selected = sample_first_x_states(root_states, config.num_sampled_states)
                selected.update(list(optimal)[0:9])
            else:
                selected = set(states)

    states, goals, transitions, optimal = remap_sample_expanded_states(set(selected), states, goal_states, transitions, optimal)

    expanded = lambda s: len(transitions[s]) > 0
    log_sampled_states(states, goals, transitions, expanded, optimal, config.resampled_states_filename)
    logging.info("Total sampled states / goals / optimal: {} / {} / {}".format(len(states), len(goals), len(optimal)))
    return states, goals, transitions, optimal


def random_sample(config, goal_states, rng, states, transitions, parents):
    num_states = min(len(states), config.num_states)
    if config.num_sampled_states > num_states:
        raise RuntimeError(
            "Only {} expanded statesm cannot sample {}".format(num_states, config.num_sampled_states))
    # Although this might fail is some state was a dead-end? In that case we should perhaps revise the code below
    assert num_states == len(transitions)

    selected = sample_expanded_and_goal_states(config, goal_states, num_states, parents, rng)
    return selected


def sample_expanded_and_goal_states(config, goal_states, num_states, parents, rng):
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
    return all_shuffled[:config.num_sampled_states]


def compute_parents(transitions):
    parents = defaultdict(set)
    for source, targets in transitions.items():
        for t in targets:
            parents[t].add(source)
    return parents


def remap_sample_expanded_states(sampled_expanded, states, goal_states, transitions, optimal):
    # all_in_sample will contain all states we want to sample plus their children state IDs
    all_in_sample = set(sampled_expanded)
    for s in sampled_expanded:
        all_in_sample.update(transitions[s])

    ordered_sample = list(sorted(all_in_sample))
    idx = {x: i for i, x in enumerate(ordered_sample)}

    # Pick the selected elements from the data structures
    new_goal_states = {idx[x] for x in ordered_sample if x in goal_states}
    new_optimal = {idx[x] for x in ordered_sample if x in optimal}

    new_states = OrderedDict()
    for i, s in states.items():
        assert i == s[0]
        if i in idx:
            new_states[idx[i]] = (idx[i], s[1])

    new_transitions = defaultdict(set)
    for source, targets in transitions.items():
        if source in sampled_expanded:
            new_transitions[idx[source]] = {idx[t] for t in targets}

    return new_states, new_goal_states, new_transitions, new_optimal


def normalize_atom_name(name):
    tmp = name.replace('()', '').replace(')', '').replace('(', ',')
    if "=" in tmp:  # We have a functional atom
        tmp = "=," + tmp.replace("=", ',')  # Mark the string putting the "=" as the first position

    return tmp.split(',')


def remap_state_ids(states, goals, transitions, remap):

    new_goals = {remap(x) for x in goals}
    new_states = OrderedDict()
    for i, s in states.items():
        assert i == s[0]
        new_states[remap(i)] = (remap(i), s[1])

    new_transitions = defaultdict(set)
    for source, targets in transitions.items():
        new_transitions[remap(source)] = {remap(t) for t in targets}

    return new_states, new_goals, new_transitions


def read_transitions(filenames):
    all_samples = [read_single_sample_file(f) for f in filenames]
    assert len(all_samples) > 0
    states, goals, transitions = all_samples[0]
    assert states[0][0] == 0
    root_states = {0}
    for s, g, tx in all_samples[1:]:
        starting_state_id = len(states)
        s, g, tx = remap_state_ids(s, g, tx, remap=lambda state: state + starting_state_id)
        assert next(iter(s)) == starting_state_id
        root_states.add(starting_state_id)
        states.update(s)
        transitions.update(tx)
        goals.update(g)

    return states, goals, transitions, root_states


def read_single_sample_file(filename):
    states_by_id = {}
    states_by_str = {}
    transitions = defaultdict(set)
    transitions_inv = defaultdict(set)
    seen = set()
    goal_states = set()

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
    return ordered, goal_states, transitions
