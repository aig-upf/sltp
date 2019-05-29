import itertools
import json
import logging

from collections import defaultdict, OrderedDict

from .util.command import read_file
from .util.naming import filename_core
from .returncodes import ExitCode


class TransitionSample:
    """ """
    def __init__(self):
        self.states = OrderedDict()
        self.transitions = defaultdict(set)
        self.parents = dict()
        # self.problem = dict()
        self.roots = set()  # The set of all roots
        self.instance_roots = []  # The root of each instance
        self.instance_goals = []  # The goals of each instance
        self.goals = set()
        self.unsolvable = set()
        self.optimal_transitions = set()
        self.expanded = set()
        self.instance = dict()  # A mapping between states and the problem instances they came from

    def add_transitions(self, states, transitions, instance_id):
        # TODO Could check for ID overlappings, but would be expensive
        self.states.update(states)
        self.transitions.update(transitions)
        self.parents.update(compute_parents(transitions))
        for s in states:
            self.instance[s] = instance_id
        # TODO: This is not correct, will fail whenever we have states in the sample that have indeed been expanded
        # TODO: but have no children:
        self.expanded.update(s for s in states if len(transitions[s]) > 0)

    def mark_as_root(self, state):
        self.roots.add(state)
        self.instance_roots.append(state)

    def mark_as_goals(self, goals):
        self.goals.update(goals)
        self.instance_goals.append(goals.copy())

    def mark_as_unsolvable(self, states):
        self.unsolvable.update(states)

    def num_states(self):
        return len(self.states)

    def num_transitions(self):
        return sum(len(x) for x in self.transitions.values())

    def mark_as_optimal(self, optimal):
        self.optimal_transitions.update(optimal)

    def compute_optimal_states(self):
        """ Return those states that lie in some transition marked as optimal """
        return set(itertools.chain.from_iterable(self.optimal_transitions))

    def info(self):
        return "roots: {}, states: {}, transitions: {} ({} optimal), goals: {}, unsolvable: {}".format(
            len(self.roots), len(self.states), self.num_transitions(), len(self.optimal_transitions),
            len(self.goals), len(self.unsolvable))

    def __str__(self):
        return "TransitionsSample[{}]".format(self.info())

    def get_sorted_state_ids(self):
        return sorted(self.states.keys())

    def resample(self, selected):
        """ Resample (i.e. project) the current sample into a new sample that will contain only states specified
        in `selected` plus their children. """
        # all_in_sample will contain all states we want to sample plus their children state IDs
        all_in_sample = set(selected)
        for s in selected:  # Add the children
            all_in_sample.update(self.transitions[s])

        # Sort the selected states and compute the remapping function
        ordered_sample = list(sorted(all_in_sample))
        remapping = {x: i for i, x in enumerate(ordered_sample)}

        # Pick the selected elements from the data structures
        goals = {remapping[x] for x in ordered_sample if x in self.goals}
        unsolvable = {remapping[x] for x in ordered_sample if x in self.unsolvable}
        roots = {remapping[x] for x in ordered_sample if x in self.roots}
        optimal = {(remapping[x], remapping[y]) for x, y in self.optimal_transitions
                   if x in remapping and y in remapping}

        instance = dict()
        states = OrderedDict()
        for i, s in self.states.items():
            if i in remapping:
                states[remapping[i]] = s
                instance[remapping[i]] = self.instance[i]

        transitions = defaultdict(set)
        for source, targets in self.transitions.items():
            if source in selected:
                transitions[remapping[source]] = {remapping[t] for t in targets}

        resampled = TransitionSample()
        resampled.add_transitions(states, transitions, 0)
        resampled.instance = instance
        resampled.mark_as_goals(goals)
        resampled.mark_as_optimal(optimal)
        resampled.mark_as_unsolvable(unsolvable)
        resampled.remapping = remapping
        _ = [resampled.mark_as_root(r) for r in roots]
        return resampled

    def get_one_goal_per_instance(self):
        return set(min(x) for x in self.instance_goals if x)


def mark_optimal(goal_states, root_states, parents):
    """ Collect all those states that lie on one arbitrary optimal path to the goal """
    optimal_transitions = set()
    for goal in goal_states:
        previous = current = goal

        while current not in root_states:
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


def log_sampled_states(sample, filename):
    # We need to recompute the parenthood relation with the remapped states to log it!
    optimal_s = set(x for x, _ in sample.optimal_transitions)

    with open(filename, 'w') as f:
        for id_, state in sample.states.items():
            parents = sample.parents[id_] if id_ in sample.parents else []
            state_parents = ", ".join(sorted(map(str, parents)))
            state_children = ", ".join(sorted(map(str, sample.transitions[id_])))
            atoms = ", ".join(print_atom(atom) for atom in state)
            is_goal = "*" if id_ in sample.goals else ""
            is_expanded = "^" if id_ in sample.expanded else ""
            is_root = "=" if id_ in sample.roots else ""
            is_optimal = "+" if id_ in optimal_s else ""
            is_unsolvable = "U" if id_ in sample.unsolvable else ""
            print("#{}{}{}{}{}{} (parents: {}, children: {}):\n\t{}".
                  format(id_, is_root, is_goal, is_optimal, is_expanded, is_unsolvable, state_parents, state_children, atoms), file=f)
    logging.info('Resampled states logged at "{}"'.format(filename))


def sample_first_x_states(root_states, sample_sizes):
    """ Sample the first sample_sizes[i] states of instance i, i.e., the integers going from root_states[i] to
    root_states[i] + sample_sizes[i] """
    sampled = set()
    assert len(root_states) == len(sample_sizes)
    for root, size in zip(root_states, sample_sizes):
        sampled.update(range(root, root+size))
    return sampled


def sample_generated_states(config, rng):
    logging.info('Loading state space samples...')
    sample, goals_by_instance = read_transitions_from_files(config.sample_files)

    if not sample.goals:
        raise RuntimeError("No goal found in the sample - increase # expanded states!")

    mark_optimal_transitions(config.optimal_selection_strategy, sample, goals_by_instance)
    logging.info("Entire sample: {}".format(sample.info()))

    # if config.num_sampled_states is None and config.max_width < 1 and not config.complete_only_wrt_optimal:
    #     return states, goal_states, transitions
    # Let's deactivate random sampling temporarily, as we're not using it
    # optimal, optimal_transitions = set(), set()
    # if config.sampling == "random":
    #     assert config.num_sampled_states is not None
    #     assert not config.complete_only_wrt_optimal
    #     selected = random_sample(config, goal_states, rng, states, transitions, parents)
    #
    # else:
    # assert config.sampling in ("all", "optimal")

    if config.num_sampled_states is not None:
        # Resample the full sample and extract only a few specified states
        selected = sample_first_x_states(sample.instance_roots, config.num_sampled_states)
        states_in_some_optimal_transition = sample.compute_optimal_states()
        selected.update(states_in_some_optimal_transition)
        sample = sample.resample(set(selected))
        logging.info("Sample after resampling: {}".format(sample.info()))

    log_sampled_states(sample, config.resampled_states_filename)
    return sample


def mark_optimal_transitions(selection_strategy, sample: TransitionSample, goals_by_instance):
    """ Marks which transitions are optimal in a transition system according to some selection criterium,F
    such as marking *all* optimal transitions, or marking just one *single* optimal path.
     """
    if selection_strategy == "arbitrary":
        # For each instance, we keep the first-reached goal, as a way of selecting an arbitrary optimal path.
        goals = sample.get_one_goal_per_instance()
        optimal = mark_optimal(goals, sample.roots, sample.parents)
        sample.mark_as_optimal(optimal)
        return

    if selection_strategy == "complete":
        assert 0, "To implement"
    raise RuntimeError("Unknown optimal selection strategy")


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
        tmp = tmp.replace("=", ',')

    return tmp.split(',')


def remap_state_ids(states, goals, transitions, unsolvable, remap):

    new_goals = {remap(x) for x in goals}
    new_unsolvable = {remap(x) for x in unsolvable}
    new_states = OrderedDict()
    for i, s in states.items():
        new_states[remap(i)] = s

    new_transitions = defaultdict(set)
    for source, targets in transitions.items():
        new_transitions[remap(source)] = {remap(t) for t in targets}

    return new_states, new_goals, new_transitions, new_unsolvable


def read_transitions_from_files(filenames):
    assert len(filenames) > 0

    goals_by_instance = []
    sample = TransitionSample()
    for instance_id, filename in enumerate(filenames, 0):
        starting_state_id = sample.num_states()
        s, g, tx, unsolv = read_single_sample_file(filename)
        assert next(iter(s.keys())) == 0  # Make sure state IDs in the sample file start by 0
        s, g, tx, unsolv = remap_state_ids(s, g, tx, unsolv, remap=lambda state: state + starting_state_id)
        assert next(iter(s)) == starting_state_id

        sample.add_transitions(s, tx, instance_id)
        sample.mark_as_root(starting_state_id)
        sample.mark_as_goals(g)
        sample.mark_as_unsolvable(unsolv)

        goals_by_instance.append(g)

    return sample, goals_by_instance


def read_single_sample_file(filename):
    states_by_id = {}
    # states_by_str = {}
    transitions = defaultdict(set)
    transitions_inv = defaultdict(set)
    seen = set()
    goal_states = set()
    unsolvable_states = set()

    def register_transition(state):
        transitions[state['parent']].add(state['id'])
        transitions_inv[state['id']].add(state['parent'])

    def register_state(state):
        # states_by_str[state['atoms_string']] = data
        states_by_id[state['id']] = state['normalized_atoms']
        seen.add(state['id'])
        if j['goal']:
            goal_states.add(state['id'])
        if 'unsolvable' in j and j['unsolvable']:
            unsolvable_states.add(state['id'])

    raw_file = [line.replace(' ', '') for line in read_file(filename)]
    for raw_line in raw_file:
        j = json.loads(raw_line)
        j['normalized_atoms'] = tuple(normalize_atom_name(atom) for atom in j['atoms'])
        # j['atoms_string'] = str(j['normalized_atoms'])

        if j['id'] in seen:
            # We hit a repeated state in the search, so we simply need to record the transition
            register_transition(j)
            continue

        # The state must be _really_ new
        # assert j['atoms_string'] not in states_by_str
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

    logging.info('%s: #lines-raw-file=%d, #states-by-id=%d, #transition-entries=%d, #transitions=%d' %
                 (filename_core(filename), len(raw_file), len(states_by_id), len(transitions),
                  sum([len(targets) for targets in transitions.values()])))

    ordered = OrderedDict()  # Make sure we return an ordered dictionary
    for id_ in sorted(states_by_id.keys()):
        ordered[id_] = states_by_id[id_]
    return ordered, goal_states, transitions, unsolvable_states


def run(config, data, rng):
    assert not data
    sample = sample_generated_states(config, rng)
    return ExitCode.Success, dict(sample=sample)
