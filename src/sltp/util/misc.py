from sltp.extensions import UniverseIndex


def update_dict(d, **kwargs):
    res = d.copy()
    res.update(kwargs)
    return res


def project_dict(obj, idxs):
    """ Project a given dictionary/set/list to only a subset of its elements, given by idxs """
    if isinstance(obj, dict):
        return {k: obj[k] for k in idxs}
    elif isinstance(obj, set):
        return {k for k in idxs if k in obj}
    elif isinstance(obj, set):
        return [o for o in obj if o in idxs]
    raise RuntimeError("Cannot project unexpected data structure '{}'".format(type(obj)))


def try_number(x):
    """ Return the number-like version of x, if x is a number, or x, otherwise"""
    try:
        return int(x)
    except ValueError:
        try:
            return float(x)
        except ValueError:
            return x


# The commented code below was used for creating a universe index from a sample
# of states, but now we create that from the actual PDDL model instead, which
# should be more reliable
# def compute_universe_from_state_set(states):
#     """ Iterate through all states and collect all possible PDDL objects """
#     universe = UniverseIndex()
#
#     for sid, state in states.items():
#         for atom in state:
#             assert atom
#             if atom[0] == "=":  # Functional atom
#                 assert len(atom) <= 4, "Cannot deal with arity>1 functions yet"
#                 [universe.add(try_number(obj)) for obj in atom[2:]]
#             else:
#                 assert len(atom) in (1, 2, 3), "Cannot deal with arity>2 predicates yet"
#                 [universe.add(try_number(obj)) for obj in atom[1:]]
#
#     universe.finish()  # No more objects possible
#     return universe


def compute_universe_from_pddl_model(language):
    """ Compute a Universe Index from the PDDL model """
    universe = UniverseIndex()
    _ = [universe.add(try_number(c.symbol)) for c in language.constants()]
    universe.finish()  # No more objects possible
    return universe
