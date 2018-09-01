#!/usr/bin/python

import sys


def read_transition_file(path):
    transitions = []
    adj_list = {}

    with open(path, 'r') as f:
        for line in f.readlines():
            state = line.split()[0]
            adj_list['s' + state] = []
            for transition in line.split()[1:]:
                transitions.append(("s" + state, "s" + transition))
                adj_list['s' + state].append("s" + transition)
    return transitions, adj_list


def read_features_file(path):
    features_per_state = {}
    with open(path, 'r') as f:
        for line in f.readlines():
            state = 's' + line.split()[0]
            features_per_state[state] = []
            num_features = len(line.split()[1:])
            for feature in line.split()[1:]:
                features_per_state[state].append(int(feature))
    return features_per_state, num_features


def read_goal_states(path):
    goal_states = []
    try:
        f = open(path, 'r')
    except IOError:
        print("Could not read file:", path)
        sys.exit()

    with f:
        for line in f.readlines():
            for state in line.split():
                goal_states.append('s' + state)
    return set(x for x in goal_states)


def read_complexity_file(path):
    feature_complexity = []
    names = []
    try:
        f = open(path, 'r')
    except IOError:
        print("Could not read file:", path)
        sys.exit()

    with f:
        for line in f.readlines():
            splitted = line.split()
            names.append(' '.join(splitted[:-1]))
            feature_complexity.append(int(splitted[-1]))
    return feature_complexity, names
