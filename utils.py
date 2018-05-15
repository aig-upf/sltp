import logging
import os
import sys
from collections import deque, defaultdict
import argparse


def check_all_equal(iterator):
    iterator = iter(iterator)
    try:
        first = next(iterator)
    except StopIteration:
        return True
    return all(first == rest for rest in iterator)


def breadth_first_search(src, objects, adj):
    # initialize table
    table = dict()
    for u in objects:
        table[u] = (0, float('inf'), None)  # (u.color = BLACK, u.d = inf, u.p = nil)
    table[src] = (1, 0, None)  # (src.color = GRAY, src.d = 0, src.p = nil)

    # perform bfs starting from src
    queue = deque()
    while queue:
        u = queue.pop()
        for v in adj[u]:
            if table[v][0] == 0:
                table[v] = (1, table[u][1] + 1, u)  # (v.color = GRAY, v.d = u.d + 1, v.p = u)
                queue.appendleft(v)
        table[u][0] = 2  # u.color = BLACK
    return table


def compute_transitive_closure(relation):
    # compute adjacency lists from extension
    adj = defaultdict(list)
    for (u, v) in relation:
        adj[u].append(v)

    # objects in extension
    objects_in_ext = set()
    for (u, v) in relation:
        objects_in_ext.add(u)
        objects_in_ext.add(v)

    # perform a breadth-first search from each object
    bfs_table = dict()
    for u in objects_in_ext:
        table = breadth_first_search(u, objects_in_ext, adj)
        bfs_table[u] = table

    # compute transitive closure (u,v) is in closure iff u == v or v.p != nil in bfs-tree(u)
    closure = []
    for u in objects_in_ext:
        for v in objects_in_ext:
            if u == v:
                closure.append((u, v))
            elif (u in bfs_table) and (v in bfs_table[u]):
                assert len(bfs_table[u][v]) == 3
                if bfs_table[u][v][2] is not None:
                    closure.append((u, v))

    return closure


def transitive_closure(elements):
    closure = set(elements)
    while True:
        closure_until_now = closure | set((x, w) for x, y in closure for q, w in closure if q == y)

        if len(closure_until_now) == len(closure):
            break

        closure = closure_until_now

    return closure


# read file line by line
def read_file(filename):
    with open(filename) as f:
        for line in f:
            yield line.rstrip('\n')


def filter_subnodes(elem, t):
    return list(filter(lambda x: type(x) == t, elem.flatten()))


def parse_arguments(args):
    parser = argparse.ArgumentParser(description="Learn generalized features and concepts")
    parser.add_argument('-k', help='Number of iterations to derive concepts and roles', action='store', default=0,
                        type=int)
    parser.add_argument('transitions', help='Name of file containing transitions (output from planner)')
    parser.add_argument('-d', '--domain', required=True, help='The PDDL domain filename')
    parser.add_argument('--debug', action='store_true', help='Print additional debugging information')
    return parser.parse_args(args)


def configure_logging(args):
    level = logging.DEBUG if args.debug else logging.INFO
    filename = os.path.basename(args.transitions)
    args.result_filename = '.'.join(filename.split('.')[:-1]) + ".{}it".format(args.k)
    filename = os.path.join('logs', args.result_filename + '.log')
    logging.basicConfig(filename=filename, filemode='w', level=level)


def bootstrap(arguments):
    args = parse_arguments(arguments)
    configure_logging(args)
    return args


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
