from collections import defaultdict

from common import bwnamer, add_bw_domain_parameters_2, ipc_instances
from sltp.util.misc import update_dict


def experiments():

    base = dict(
        domain_dir="blocks",
        domain="domain.pddl",
        test_domain="domain.pddl",
        complete_only_wrt_optimal=True,
    )

    exps = dict()

    exps["on_x_y"] = update_dict(
        base,
        instances=["inst_on_x_y_14.pddl", "inst_on_x_y_16.pddl", ],
        num_states=40000, max_width=[-1],
        num_sampled_states=[500],
        # Note: Testing here is not simple, as we'd want to test only when X and Y are on different towers
        # test_instances=["inst_on_x_y_16.pddl"], num_tested_states=10000,
        # test_policy_instances=ipc_instances(),
        complete_only_wrt_optimal=True,
        max_concept_size=8,
        concept_generator=None,
        parameter_generator=add_bw_domain_parameters_2,
        feature_namer=bwnamer,
    )

    exps["on_x_y_gc"] = update_dict(exps["on_x_y"], parameter_generator=None)

    # A small incremental experiment mostly for testing purposes
    exps["small_inc"] = update_dict(
        base,
        experiment_type='incremental',
        instances=["probBLOCKS-4-0.pddl", ],
        test_instances=["probBLOCKS-5-0.pddl", ],
        num_states=10000,
        num_sampled_states=None,  # Take all expanded states into account
        initial_sample_size=5,
        initial_concept_bound=5, max_concept_bound=12, concept_bound_step=1,
        batch_refinement_size=1,
        clean_workspace=False,
    )

    exps["one_tower"] = update_dict(
        base,
        instances=["probBLOCKS-4-0.pddl"],
        test_instances=[
            "probBLOCKS-5-0.pddl",
            "probBLOCKS-5-1.pddl",
            "probBLOCKS-6-0.pddl",
            "probBLOCKS-6-1.pddl",
            # "probBLOCKS-10-0.pddl",
            # "probBLOCKS-10-1.pddl",
            # "probBLOCKS-10-2.pddl",
            # "probBLOCKS-14-0.pddl",
            # "probBLOCKS-14-1.pddl",
        ],
        test_policy_instances=ipc_instances(),
        num_states=2000,
        num_tested_states=50000,
        num_sampled_states=300,
        complete_only_wrt_optimal=True,
        max_concept_size=6,
        concept_generator=None,
        parameter_generator=None,
        feature_namer=bwnamer,
    )

    # Goal: build one single given tower of blocks from a random initial configuration
    exps["one_tower_inc"] = update_dict(
        exps["one_tower"],
        experiment_type='incremental',
        # instances=["probBLOCKS-4-0.pddl", "probBLOCKS-7-0.pddl", "probBLOCKS-8-1.pddl", "probBLOCKS-9-2.pddl", ],
        instances=["probBLOCKS-4-0.pddl", "probBLOCKS-4-1.pddl", "probBLOCKS-5-0.pddl", "probBLOCKS-5-1.pddl", "probBLOCKS-6-0.pddl", "probBLOCKS-6-1.pddl"],
        test_instances=[
            # "probBLOCKS-4-1.pddl",
            # "probBLOCKS-6-1.pddl",
            # "probBLOCKS-10-1.pddl",
        ],
        num_states=40000,
        num_sampled_states=None,  # Take all expanded states into account
        initial_sample_size=100, batch_refinement_size=20,
        initial_concept_bound=6, max_concept_bound=10, concept_bound_step=1,
        clean_workspace=False,
        # quiet=True,
    )

    exps["one_tower_inc_all"] = update_dict(
        exps["one_tower_inc"],
        num_states="all",
    )

    exps["one_tower_inc_all_g"] = update_dict(
        exps["one_tower_inc"],
        num_states="all",

        # goal_selector=goal_selector,
        goal_selector=maximize_num_blocks_on_same_tower,
        create_goal_features_automatically=True,
    )

    exps["one_tower_goal_selection"] = update_dict(
        exps["one_tower_inc"],
        instances=["probBLOCKS-7-0.pddl"],
        num_states=40000,
        goal_selector=maximize_num_blocks_on_same_tower,
        create_goal_features_automatically=False,
        maxsat_encoding="basic",  # d2tree seems to be too expensive!
    )

    exps["one_tower_test"] = update_dict(
        exps["one_tower"],
        instances=["probBLOCKS-4-0.pddl", "probBLOCKS-4-1.pddl", "probBLOCKS-5-0.pddl", "probBLOCKS-5-1.pddl", "probBLOCKS-6-0.pddl", "probBLOCKS-6-1.pddl"],
        test_instances=[
            # "probBLOCKS-10-0.pddl",
            # "probBLOCKS-10-1.pddl",
            # "probBLOCKS-10-2.pddl",
            # "probBLOCKS-11-0.pddl",
            # "probBLOCKS-11-1.pddl",
            # "probBLOCKS-11-2.pddl",
        ],
        num_states=40000,
        num_sampled_states=1000,
        feature_generator=feature_testing,
    )

    return exps


def feature_testing(lang):
    return [
        "Atom[handempty]",
        "Bool[And(Equal(on_g,on),Forall(on,<empty>))]",
        "Bool[And(Equal(Star(on_g),Star(Inverse(on))),Forall(on_g,<empty>))]",
        "Num[Exists(Star(on_g),clear)]",
        "Num[Exists(Star(on_g),And(Equal(Star(on_g),Star(on)),clear))]",

    ]


def genfeatures(lang):
    return [
        "Bool[holding]",
        "Bool[And(Equal(on_g,on),Forall(on,<empty>))]",
        "Num[Equal(Star(on_g),Star(Inverse(on)))]",
        "Num[Exists(on,Equal(Star(on_g),Star(Inverse(on))))]",
        "Num[Exists(Star(on_g),clear)]",
        "Num[Exists(Star(on_g),And(Equal(Star(on_g),Star(on)),clear))]",
    ]


def goal_selector(lang):
    # return "RoleDifference(on_g,on)"
    return "Equal(on_g,on)"


def maximize_num_blocks_on_same_tower(lang, state):
    atoms = [a for a in state if a[0] != 'handempty' and a[0] != 'clear']
    inv_on = {a[2]: a[1] for a in atoms if a[0] == 'on'}  # on^-1 relation
    tower_counts = defaultdict(int)

    for a in atoms:
        if a[0] == 'ontable':
            base = theoneabove = a[1]
            while theoneabove is not None:
                tower_counts[base] += 1
                theoneabove = inv_on.get(theoneabove, None)

    return max(tower_counts.values())


maximize_num_blocks_on_same_tower.procedural = True
