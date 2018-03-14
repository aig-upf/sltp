"""
A few utility functions to parse the output of the FS planner and extract the relevant information
about expanded states.
"""

import json
from collections import OrderedDict

from utils import read_file

# read transitions
g_states_by_str = {}
g_transitions = {}


def normalize_atom_name(name):
    return name.replace('()', '').rstrip(')').replace('(', ',').split(',')


def read_transitions(transitions_filename):
    states_by_id = {}
    
    raw_file = [line.replace(' ', '') for line in read_file(transitions_filename) if line[0:6] == '{"id":']
    for raw_line in raw_file:
        j = json.loads(raw_line)
        j_atoms = [normalize_atom_name(atom) for atom in j['atoms']]
        j_atoms_str = str(j_atoms)

        # insert state into hash with (normalized) id
        if j_atoms_str not in g_states_by_str:
            j_id = int(j['id'])
            g_states_by_str[j_atoms_str] = states_by_id[j_id] = (j_id, j_atoms)
        else:
            j_id = g_states_by_str[j_atoms_str][0]

        # insert (normalized) transition into hash
        if j['parent'] != j['id']:
            j_pid = int(j['parent'])
            jp = json.loads(raw_file[j_pid])
            assert jp['id'] == j_pid
            jp_atoms = [normalize_atom_name(atom) for atom in jp['atoms']]

            jp_atoms_str = str(jp_atoms)

            assert jp_atoms_str in g_states_by_str
            jp_id = g_states_by_str[jp_atoms_str][0]
            if jp_id not in g_transitions: g_transitions[jp_id] = []
            g_transitions[jp_id].append(j_id)

    # check soundness
    for src in g_transitions:
        assert src in states_by_id
        for dst in g_transitions[src]:
            assert dst in states_by_id

    print('#lines-raw-file=%d, #state-by-str=%d, #states-by-id=%d, #transition-entries=%d, #transitions=%d' % (
        len(raw_file), len(g_states_by_str), len(states_by_id), len(g_transitions),
        sum([len(g_transitions[src]) for src in g_transitions])))

    ordered = OrderedDict()  # Make sure we return an ordered dictionary
    for id_ in sorted(states_by_id.keys()):
        ordered[id_] = states_by_id[id_]
    return ordered
