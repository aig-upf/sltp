import logging


from tarski.io import FstripsReader
from tarski.syntax.transform.errors import TransformationError
from tarski.syntax.transform.simplifications import transform_to_ground_atoms


def parse_pddl(domain_file, instance_file=None):
    """ Parse the given PDDL domain and instance files, the latter only if provided """
    reader = FstripsReader()
    reader.parse_domain(domain_file)
    problem = reader.problem

    # The generic constants are those which are parsed before parsing the instance file
    generic_constants = problem.language.constants()

    if instance_file is not None:
        reader.parse_instance(instance_file)

    return problem, problem.language, generic_constants


def compute_goal_denotation(problem, use_goal_denotation):
    """ Compute the goal-version of the relevant concepts from the problem goal specification """
    if not use_goal_denotation:
        return []

    try:
        return transform_to_ground_atoms(problem.goal)
    except TransformationError:
        logging.error("Cannot create goal concepts when problem goal is not a conjunction of ground atoms")
        raise

