import logging


def print_state_set(goal_states, filename):
    with open(filename, 'w') as f:
        print(" ".join(str(s) for s in goal_states), file=f)


def print_transition_matrix(state_ids, transitions, transitions_filename):
    num_transitions = sum(len(targets) for targets in transitions.values())
    logging.info("Printing transition matrix with {} states and {} transitions to '{}'".
                 format(len(state_ids), num_transitions, transitions_filename))
    with open(transitions_filename, 'w') as f:
        for s, succ in transitions.items():
            print("{} {}".format(s, " ".join("{}".format(sprime) for sprime in succ)), file=f)


def print_sat_transition_matrix(state_ids, transitions, transitions_filename):
    num_transitions = sum(len(targets) for targets in transitions.values())
    num_states = len(transitions.keys())
    logging.info("Printing SAT transition matrix with {} states, {} expanded states and {} transitions to '{}'".
                 format(len(state_ids), num_states, num_transitions, transitions_filename))
    with open(transitions_filename, 'w') as f:
        # first line: <#states> <#transitions>
        print("{} {}".format(num_states, num_transitions), file=f)

        # second line: <#expanded states>
        print("{}".format(num_states), file=f)

        for s, succ in transitions.items():
            print("{} {} {}".format(s, len(succ), " ".join("{}".format(sprime) for sprime in sorted(succ))), file=f)

