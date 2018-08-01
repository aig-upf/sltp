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
                                   instance=os.path.join(BENCHMARK_DIR, domain_dir, "instance_5_clear_x.pddl"),
                                   driver="bfs",
                                   planner_location=os.getenv("FS_PATH", os.path.expanduser("~/projects/code/fs")),
                                   num_states=100,
                                   concept_depth=1,
                                   concept_generator=build_paper_concepts,
                                   use_distance_features=False,
                                   optimization=OptimizationPolicy.TOTAL_FEATURE_COMPLEXITY,
                                   # optimization=OptimizationPolicy.NUM_FEATURES
                                   )
    exp = Experiment(steps)
    exp.run()


def build_paper_concepts(lang):
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
