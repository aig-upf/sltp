#! /usr/bin/env python3
# -*- coding: utf-8 -*-
import os

from tarski.dl import AndConcept, NominalConcept, PrimitiveConcept, NotConcept, ExistsConcept, InverseRole, \
    PrimitiveRole, StarRole, ForallConcept, GoalRole, EqualConcept


def main():
    import sys
    sys.path.insert(0, '..')
    from driver import Experiment, generate_pipeline, BENCHMARK_DIR
    from learn_actions import OptimizationPolicy

    domain_dir = "blocks"
    # domain = "domain.pddl"
    # instance = "probBLOCKS-3-0.pddl"
    domain = "domain_fstrips.pddl"
    instance = "instance_5_fstrips.pddl"

    steps = generate_pipeline(pipeline="maxsat",
                              domain=os.path.join(BENCHMARK_DIR, domain_dir, domain),
                              instance=os.path.join(BENCHMARK_DIR, domain_dir, instance),

                              # Location of the FS planner, used to do the state space sampling
                              planner_location=os.getenv("FS_PATH", os.path.expanduser("~/projects/code/fs")),

                              # Type of sampling procedure. Only breadth-first search implemented ATM
                              driver="bfs",

                              # Number of states to be expanded in the sampling procedure
                              num_states=5000,

                              max_concept_size=10,

                              # Provide a special, handcrafted method to generate concepts, if desired.
                              # This will override the standard concept generation procedure (default: None)
                              # concept_generator=generate_chosen_concepts,

                              # Max. allowed complexity for distance features (default: 0)
                              distance_feature_max_complexity=0,

                              # Method to generate domain parameters (goal or otherwise). If None, goal predicates will
                              # be used (default: None)
                              # parameter_generator=add_domain_parameters,

                              # What optimization criteria to use in the max-sat problem
                              optimization=OptimizationPolicy.TOTAL_FEATURE_COMPLEXITY,
                              # optimization=OptimizationPolicy.NUM_FEATURES
                              )
    exp = Experiment(steps)
    exp.run()


def generate_chosen_concepts(lang):
    """ Return the set concepts from Hector & Blai's IJCAI paper which are sound and complete for clear(X) goals """
    obj_t = lang.Object
    holding_p = lang.get_predicate("holding")

    ontable = PrimitiveRole(lang.get("ontable"))
    on_r = PrimitiveRole(lang.get("on"))
    on_g = GoalRole(lang.get("on"))
    on_g_closure = StarRole(on_g)

    y = EqualConcept(on_g, on_r, "object")
    x = ForallConcept(on_g_closure, y)
    well_placed_block = ForallConcept(on_g, x)
    misplaced_block = NotConcept(well_placed_block, obj_t)
    misplaced_and_ontable = AndConcept(misplaced_block, ontable, "object")

    concepts = [misplaced_block, misplaced_and_ontable]
    return [], concepts, []  # atoms, concepts, roles


def add_domain_parameters(language):
    # We simply add block "a" as a domain constant
    return [language.constant("a", "object")]
    # language.constant("b", "object")


if __name__ == "__main__":
    main()
