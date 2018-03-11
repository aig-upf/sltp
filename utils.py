from collections import deque, defaultdict


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
