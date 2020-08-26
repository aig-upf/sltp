import logging
import itertools


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


def print_sat_transition_matrix(state_ids, transitions, alive, vstar, optimal_transitions, transitions_filename):
    num_transitions = sum(len(targets) for targets in transitions.values())
    num_states = len(transitions.keys())
    num_optimal_transitions = len(optimal_transitions)
    logging.info(f"Printing SAT transition matrix with {len(state_ids)} states,"
                 f" {num_states} expanded states and {num_transitions} transitions to '{transitions_filename}'")

    # State Ids should start at 0 and be contiguous
    assert state_ids == list(range(0, len(state_ids)))

    with open(transitions_filename, 'w') as f:
        # first line: <#states> <#transitions> <#marked-transitions>
        print(f"{num_states} {num_transitions} {num_optimal_transitions}", file=f)

        # second line: all marked transitions, in format: <src0> <dst0> <src1> <dst1> ...
        print(" ".join(map(str, itertools.chain(*optimal_transitions))), file=f)

        # third line: <#expanded states>
        print(f"{num_states}", file=f)

        # Next lines: transitions, one per source state, in format: source_id, num_successors, succ_1, succ_2, ...
        for s, successors in transitions.items():
            sorted_succs = " ".join(map(str, sorted(successors)))
            print(f"{s} {len(successors)} {sorted_succs}", file=f)

        # Next: A list of all alive states
        print(f'{len(alive)} ' + " ".join(map(str, alive)), file=f)

        # Next: A space-separated list of V^*(s) values, one per each state s, where -1 denotes infinity
        print(' '.join(str(vstar.get(s, -1)) for s in state_ids), file=f)
