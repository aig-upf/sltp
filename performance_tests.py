from timeit import timeit

setup = """
from random import sample, shuffle
universe = list(range(1000))
shuffle(universe)
a = tuple(sample(universe, 50))
b = tuple(sample(universe, 50))
setu = set(universe)
seta = set(a)
setb = set(b)
"""

forin = setup + """
def forin():
    return len(list(elem for elem in b if elem in a))
"""

setin = setup + """
def setin():
    return len(seta & setb)
"""

setcm = setup + """
def setcm():
    return len([x for x in seta if x in setb])
"""

# print(timeit("forin()", forin, number=10000))
# print(timeit("setin()", setin, number=10000))
# print(timeit("setcm()", setcm, number=10000))

# Compute the complement of A
neg_lst = setup + """
def neg_lst():
    return len(list(elem for elem in universe if elem not in a))
"""

neg_set = setup + """
def neg_set():
    return len(setu - seta)
"""

neg_cmp = setup + """
def neg_cmp():
    return len([x for x in setu if x not in seta])
"""

times = 1000
print(timeit("neg_lst()", neg_lst, number=times))
print(timeit("neg_set()", neg_set, number=times))
print(timeit("neg_cmp()", neg_cmp, number=times))
