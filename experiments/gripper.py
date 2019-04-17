
from tarski.dl import ConceptCardinalityFeature, NominalConcept, ExistsConcept, PrimitiveRole, UniversalConcept, \
    NotConcept, AndConcept, PrimitiveConcept, EmpiricalBinaryConcept

from sltp.util.misc import update_dict


def experiments():
    domain = "domain.pddl"
    base = dict(
        domain_dir="gripper",
        domain=domain,
    )

    exps = dict()

    exps["sample_small"] = update_dict(base,
        instances="sample-small.pddl",
        test_domain=domain, test_instances=["prob03.pddl", "prob04.pddl"],
        num_states=100, max_concept_size=10, max_concept_grammar_iterations=3,
        concept_generator=None, parameter_generator=add_domain_parameters,
        feature_namer=feature_namer,)

    exps["prob01"] = update_dict(base,
        instances="prob01.pddl",
        num_states=300, num_sampled_states=None, random_seed=12,
        max_concept_size=10, max_concept_grammar_iterations=3,
        concept_generator=None, parameter_generator=add_domain_parameters,
        feature_namer=feature_namer,)

    exps["prob01_goalc"] = update_dict(base,
        instances=["prob01.pddl", "sample02.pddl", ],
        num_states=300, num_sampled_states=None, random_seed=12,
        max_concept_size=10, max_concept_grammar_iterations=3,
        concept_generator=None, parameter_generator=None,
        feature_namer=feature_namer,)

    #
    exps["prob01_rnd"] = update_dict(base,
        instances="prob01.pddl",
        num_states=2000, num_sampled_states=50, random_seed=12,
        max_concept_size=10, max_concept_grammar_iterations=3,
        concept_generator=None, parameter_generator=add_domain_parameters,
        feature_namer=feature_namer,)

    #
    exps["aaai_prob01"] = update_dict(base,
        instances=["prob01.pddl", "sample02.pddl"],
        test_domain=domain, test_instances=["prob03.pddl"],
        num_tested_states=5000,
        num_states=2000, max_width=[-1],
        num_sampled_states=100,
        complete_only_wrt_optimal=True,
        max_concept_size=8, max_concept_grammar_iterations=3,
        concept_generator=None, parameter_generator=add_domain_parameters,
        feature_namer=feature_namer,)

    # Same but using goal-concepts instead of goal parameters:
    exps["aaai_prob01_gc"] = update_dict(exps["aaai_prob01"], parameter_generator=None)

    exps["aaai_prob01_no_marking"] = update_dict(exps["aaai_prob01"], complete_only_wrt_optimal=False)

    exps["aaai_prob01_blai"] = update_dict(
        exps["aaai_prob01"], pipeline="maxsat_poly",
        # max_concept_size=4, max_concept_grammar_iterations=2,
        # num_states=100,
    )

    exps["aaai_prob01_blai_std"] = update_dict(  # Same config as Blai, but with standard pipeline
        exps["aaai_prob01_blai"], pipeline="maxsat")

    exps["aaai_prob01_deb"] = update_dict(exps["aaai_prob01"], feature_generator=debug_aaai_features)

    return exps


def add_domain_parameters(language):
    return [language.constant("roomb", "object")]


def feature_namer(feature):
    s = str(feature)
    return {
        "card[Exists(at,Not({roomb}))]": "nballs-A",
        "card[Exists(at,{roomb})]": "nballs-B",
        "card[Exists(carry,<universe>)]": "ncarried",
        "bool[And(at-robby,{roomb})]": "robot-at-B",
        "card[Exists(at,Not(at-robby))]": "nballs-in-rooms-with-no-robot",
        "card[free]": "nfree-grippers",
        "bool[Exists(at-robby,{roomb})]": "robby-is-at-B",
        "card[Exists(at,Exists(Inverse(at-robby),<universe>))]": "nballs-in-room-with-some-robot",
        "card[And(Exists(gripper,Exists(at-robby,{roomb})),free)]": "nfree-grippers-at-B",
        "card[Exists(at-robby,{roomb})]": "nrobots-at-B",
        "card[Exists(gripper,Exists(at-robby,{roomb}))]": "ngrippers-at-B",
        "card[Exists(carry,Exists(gripper,Exists(at-robby,{roomb})))]": "nballs-carried-in-B",
        "card[Exists(at,And(Forall(Inverse(at-robby),<empty>), Not({roomb})))]": "nballs-in-some-room-notB-without-any-robot",
        "bool[And(Exists(Inverse(at),<universe>), And({roomb}, Not(at-robby)))]": "some-ball-in-B-but-robbot-not-in-B",
    }.get(s, s)


def debug_aaai_features(lang):
    obj_t = lang.Object
    at = PrimitiveRole(lang.get("at"))
    carry = PrimitiveRole(lang.get("carry"))
    at_robby = PrimitiveConcept(lang.get("at-robby"))
    free = PrimitiveConcept(lang.get("free"))
    x_param = NominalConcept("roomb", obj_t)

    f1 = AndConcept(at_robby, x_param, "object")
    f2 = ExistsConcept(at, NotConcept(x_param, obj_t))
    f3 = ExistsConcept(carry, UniversalConcept("object"))

    return [
        EmpiricalBinaryConcept(ConceptCardinalityFeature(f1)),
        ConceptCardinalityFeature(f2),
        ConceptCardinalityFeature(f3),
        ConceptCardinalityFeature(free),
    ]
