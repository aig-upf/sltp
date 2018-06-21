"""
A few utility functions to parse the output of the FS planner and extract the relevant information
about expanded states.
"""

import json
import logging
from collections import OrderedDict, defaultdict

from util.command import read_file


def normalize_atom_name(name):
    return name.replace('()', '').rstrip(')').replace('(', ',').split(',')


def read_transitions(transitions_filename):
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

    raw_file = [line.replace(' ', '') for line in read_file(transitions_filename) if line[0:6] == '{"id":']
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

        if j['parent'] != j['id']:  # i.e. if not in the root node, which is not the target of any transition
            assert json.loads(raw_file[j['parent']])['id'] == j['parent']  # Just a check
            register_transition(j)

    # check soundness
    for src in transitions:
        assert src in states_by_id
        for dst in transitions[src]:
            assert dst in states_by_id

    assert sum([len(t) for t in transitions.values()]) == sum([len(t) for t in transitions_inv.values()])

    logging.info('#lines-raw-file=%d, #state-by-str=%d, #states-by-id=%d, #transition-entries=%d, #transitions=%d' % (
        len(raw_file), len(states_by_str), len(states_by_id), len(transitions),
        sum([len(targets) for targets in transitions.values()])))

    ordered = OrderedDict()  # Make sure we return an ordered dictionary
    for id_ in sorted(states_by_id.keys()):
        ordered[id_] = states_by_id[id_]
    return ordered, goal_states, transitions
