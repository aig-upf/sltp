"""
This hacky module is used to identify, on a per-domain (and possibly per-instance) basis,
which of the objects in a certain problem are also parameters, i.e. can be considered to
be objects of all instances in the domain.

It is not too elegant, but spares us the need to manually modify every domain/instance we want to
test our approach on.
"""
import logging


def add_domain_goal_parameters(domain_name, language):
    generic_constants = []

    if domain_name == "blocks":
        # We simply add block "a" as a domain constant
        generic_constants = [language.constant("a", "object")]
        # language.constant("b", "object")

    elif domain_name == "gripper-strips":
        # language.constant("ball1", "object")
        generic_constants = [language.constant("roomb", "object")]

    elif domain_name == "gridworld":
        # language.constant(2, "coordinate")  # x-goal coordinate
        # language.constant(3, "coordinate")  # x-goal coordinate
        # language.constant(10, "coordinate")  # grid limits!!
        generic_constants = [language.constant(1, "coordinate")]  # grid limits!!

        # [language.constant(i, "coordinate") for i in range(1, 11)]

    elif domain_name == "gridworld-strips":
        # language.constant("c1", "coordinate")  # grid limits!!
        generic_constants = []

    elif domain_name == "grid-visit-all":
        generic_constants = [language.constant("loc-x0-y0", "place"),
                             language.constant("loc-x0-y1", "place")]

    else:
        logging.critical('Domain name "{}" not recognized, no goal parameters added'.format(domain_name))

    return generic_constants
