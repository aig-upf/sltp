import sys

from collections.__init__ import deque, defaultdict


def check_all_equal(iterator):
    iterator = iter(iterator)
    try:
        first = next(iterator)
    except StopIteration:
        return True
    return all(first == rest for rest in iterator)


def transitive_closure(elements):
    closure = set(elements)
    while True:
        closure_until_now = closure | set((x, w) for x, y in closure for q, w in closure if q == y)

        if len(closure_until_now) == len(closure):
            break

        closure = closure_until_now

    return closure


def filter_subnodes(elem, t):
    return list(filter(lambda x: type(x) == t, elem.flatten()))


def compute_min_distance(c1s, relation, c2s):
    """  """
    # Cover first a couple of base cases to enhance performance
    if c1s & c2s:
        return 0
    if not c1s or not c2s or not relation:
        return sys.maxsize

    # index the adjacency relation
    adjacencies = defaultdict(list)
    [adjacencies[s].append(t) for s, t in relation]

    # initialize table of distances
    min_distances = defaultdict(lambda: sys.maxsize)
    queue = deque()
    for obj in c1s:
        queue.appendleft(obj)
        min_distances[obj] = 0

    # Perform breadth-first search
    while queue:
        s = queue.pop()
        dist_through_s = min_distances[s]+1
        for t in adjacencies[s]:
            if dist_through_s < min_distances[t]:
                min_distances[t] = min_distances[s]+1
                queue.appendleft(t)

    # Return the minimum distance
    return min(min_distances[o] for o in c2s)