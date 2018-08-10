import sys


def hill_climbing(s0, adjacencies, heuristic, goal_states):
    """ """

    # counts the number of loops (only for printing)
    s = s0
    plan = [s]
    while s not in goal_states:

        h_s = heuristic(s)
        improvement_found = False
        for succ in adjacencies[s]:
            h_succ = heuristic(succ)
            if h_succ < h_s:
                s = succ
                plan.append(s)
                improvement_found = True
                break

        if not improvement_found:
            # Error!
            successor_values = [(succ, heuristic(succ)) for succ in adjacencies[s]]
            print("", file=sys.stderr)
            print("h({})={}, but:".format(s, h_s), file=sys.stderr)
            print("\t" + "\n\t".join("h({})={}".format(succ, hsucc) for succ, hsucc in successor_values), file=sys.stderr)
            sys.stderr.flush()
            raise RuntimeError("Your model doesn't work")


    assert s in goal_states
    print("\nGenius: ")
    print(", ".join(plan))
