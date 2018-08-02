#! /usr/bin/env python3
# -*- coding: utf-8 -*-
import os

from tarski.dl import AndConcept, NominalConcept, PrimitiveConcept, NotConcept, ExistsConcept, InverseRole, \
    PrimitiveRole, StarRole


def main():
    import sys
    sys.path.insert(0, '..')
    from driver import Experiment, generate_full_pipeline, BENCHMARK_DIR
    from learn_actions import OptimizationPolicy

    domain_dir = "blocks"
    steps = generate_full_pipeline(domain=os.path.join(BENCHMARK_DIR, domain_dir, "domain.pddl"),
                                   instance=os.path.join(BENCHMARK_DIR, domain_dir, "probBLOCKS-4-0.pddl"),

                                   # Location of the FS planner, used to do the state space sampling
                                   planner_location=os.getenv("FS_PATH", os.path.expanduser("~/projects/code/fs")),

                                   # Type of sampling procedure. Only breadth-first search implemented ATM
                                   driver="bfs",

                                   # Number of states to be expanded in the sampling procedure
                                   num_states=90,

                                   # Number of iterations of the concept generation grammar.
                                   concept_depth=2,

                                   # Provide a special, handcrafted method to generate concepts, if desired.
                                   # This will override the standard concept generation procedure (default: None)
                                   # concept_generator=build_paper_concepts,

                                   # Whether to use distance features (default: False)
                                   use_distance_features=True,

                                   # Whether to use concepts based on goal predicate denotations (default: False)
                                   use_goal_features=True,

                                   # What optimization criteria to use in the max-sat problem
                                   optimization=OptimizationPolicy.TOTAL_FEATURE_COMPLEXITY,
                                   # optimization=OptimizationPolicy.NUM_FEATURES
                                   )
    exp = Experiment(steps)
    exp.run()


def build_ijcai_paper_concepts(lang):
    """ Return the concepts from Hector & Blai's IJCAI paper """
    obj_t = lang.Object
    x_nominal = NominalConcept("a", obj_t)
    holding_p = lang.get_predicate("holding")
    on_p = lang.get_predicate("on")
    on_r = PrimitiveRole(on_p)
    inv_on = InverseRole(on_r)
    on_star = StarRole(on_r)
    inv_on_star = StarRole(inv_on)

    holding = PrimitiveConcept(holding_p)
    not_held = NotConcept(holding, obj_t)
    x_held = AndConcept(x_nominal, holding, "object")  # H: Block x is being held

    other_block = NotConcept(x_nominal, obj_t)
    other_block_held = AndConcept(other_block, holding, "object")

    there_is_a_block_below_x = ExistsConcept(inv_on, x_nominal)

    blocks_above_x = ExistsConcept(on_star, x_nominal)
    blocks_below_x = ExistsConcept(inv_on_star, x_nominal)

    not_above_x = NotConcept(blocks_above_x, obj_t)
    not_below_x = NotConcept(blocks_below_x, obj_t)
    not_x = NotConcept(x_nominal, obj_t)

    blocks_not_in_x_tower = AndConcept(AndConcept(not_above_x, not_below_x, "object"), not_x, "object")
    blocks_not_in_x_tower_or_held = AndConcept(blocks_not_in_x_tower, not_held, "object")

    concepts = [x_held, other_block_held, there_is_a_block_below_x, blocks_above_x, blocks_not_in_x_tower_or_held]
    return [], concepts, []  # atoms, concepts, roles


if __name__ == "__main__":
    main()
