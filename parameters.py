"""
This hacky module is used to identify, on a per-domain (and possibly per-instance) basis,
which of the objects in a certain problem are also parameters, i.e. can be considered to
be objects of all instances in the domain.

It is not too elegant, but spares us the need to manually modify every domain/instance we want to
test our approach on.
"""


def add_domain_goal_parameters(domain_name, language):

    if domain_name == "blocks":
        # We simply add block "a" as a domain constant
        language.constant("a", "object")
        # language.constant("b", "object")

    elif domain_name == "gripper-strips":
        # language.constant("ball1", "object")
        language.constant("roomb", "object")

    elif domain_name == "grid-visit-all":
        language.constant("loc-x0-y0", "place")
        language.constant("loc-x0-y1", "place")
    else:
        print('WARNING: Domain name "{}" not recognized, no goal parameters added'.format(domain_name))
    return language
