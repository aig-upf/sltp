
import itertools

from sltp.util.misc import extend_namer_to_all_features
from tarski.dl import NominalConcept, PrimitiveRole, InverseRole, StarRole, PrimitiveConcept, NotConcept, AndConcept, \
    ExistsConcept, ConceptCardinalityFeature


def build_ijcai_paper_bw_concepts(lang):
    """ Return the set concepts from Hector & Blai's IJCAI paper which are sound and complete for clear(X) goals """
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


def build_alt_feature_set(lang):
    _, ijcai, _ = build_ijcai_paper_bw_concepts(lang)
    modified = ijcai[:-1]
    return [], modified, []


def build_on_x_y_feature_set(lang):
    obj_t = lang.Object
    x_nominal = NominalConcept("a", obj_t)
    y_nominal = NominalConcept("b", obj_t)
    holding_p = lang.get_predicate("holding")
    on_p = lang.get_predicate("on")
    on_r = PrimitiveRole(on_p)
    inv_on = InverseRole(on_r)
    on_star = StarRole(on_r)
    inv_on_star = StarRole(inv_on)

    holding = PrimitiveConcept(holding_p)
    not_held = NotConcept(holding, obj_t)
    x_held = AndConcept(x_nominal, holding, "object")  # H: Block x is being held
    y_held = AndConcept(y_nominal, holding, "object")

    other_block = AndConcept(NotConcept(x_nominal, obj_t), NotConcept(y_nominal, obj_t), "object")
    other_block_held = AndConcept(other_block, holding, "object")

    there_is_a_block_below_x = ExistsConcept(inv_on, x_nominal)
    there_is_a_block_below_y = ExistsConcept(inv_on, y_nominal)

    blocks_above_x = ExistsConcept(on_star, x_nominal)
    blocks_below_x = ExistsConcept(inv_on_star, x_nominal)
    blocks_above_y = ExistsConcept(on_star, y_nominal)
    blocks_below_y = ExistsConcept(inv_on_star, y_nominal)

    on_x_y = AndConcept(ExistsConcept(on_r, y_nominal), x_nominal, "object")
    on_y_x = AndConcept(ExistsConcept(on_r, x_nominal), y_nominal, "object")
    above_x_y = AndConcept(blocks_above_y, x_nominal, "object")
    above_y_x = AndConcept(blocks_above_x, y_nominal, "object")

    not_above_x = NotConcept(blocks_above_x, obj_t)
    not_below_x = NotConcept(blocks_below_x, obj_t)
    not_x = NotConcept(x_nominal, obj_t)
    not_y = NotConcept(y_nominal, obj_t)

    blocks_not_in_x_tower = AndConcept(AndConcept(not_above_x, not_below_x, "object"), not_x, "object")
    blocks_not_in_x_tower_or_held = AndConcept(blocks_not_in_x_tower, not_held, "object")

    blocks_not_in_y_tower = AndConcept(AndConcept(NotConcept(blocks_above_y, obj_t),
                                                  NotConcept(blocks_below_y, obj_t), "object"),
                                       not_y, "object")

    num_other_blocks = AndConcept(blocks_not_in_x_tower_or_held, blocks_not_in_y_tower, "object")

    concepts = [x_held, y_held, other_block_held,
                above_y_x, above_x_y, on_x_y, on_y_x,
                there_is_a_block_below_x, blocks_above_x,
                there_is_a_block_below_y, blocks_above_y,
                num_other_blocks]
    return [], concepts, []  # atoms, concepts, roles


def generate_features_n_ab(lang):
    obj_t = lang.Object
    x_nominal = NominalConcept("a", obj_t)
    y_nominal = NominalConcept("b", obj_t)
    holding_p = lang.get_predicate("holding")
    on_p = lang.get_predicate("on")
    on_r = PrimitiveRole(on_p)
    on_star = StarRole(on_r)

    holding = PrimitiveConcept(holding_p)
    x_held = AndConcept(x_nominal, holding, "object")  # H: Block x is being held
    y_held = AndConcept(y_nominal, holding, "object")

    other_block = AndConcept(NotConcept(x_nominal, obj_t), NotConcept(y_nominal, obj_t), "object")
    other_block_held = AndConcept(other_block, holding, "object")

    blocks_above_x = ExistsConcept(on_star, x_nominal)
    blocks_above_y = ExistsConcept(on_star, y_nominal)

    on_x_y = AndConcept(ExistsConcept(on_r, y_nominal), x_nominal, "object")
    on_y_x = AndConcept(ExistsConcept(on_r, x_nominal), y_nominal, "object")

    # concepts = [x_held, y_held, other_block_held,
    #             above_y_x, above_x_y, on_x_y, on_y_x,
    #             there_is_a_block_below_x, blocks_above_x,
    #             there_is_a_block_below_y, blocks_above_y,
    #             num_other_blocks]

    return [
        ConceptCardinalityFeature(on_x_y),
        ConceptCardinalityFeature(on_y_x),
        ConceptCardinalityFeature(blocks_above_x),
        ConceptCardinalityFeature(blocks_above_y),
        ConceptCardinalityFeature(x_held),
        ConceptCardinalityFeature(y_held),
        ConceptCardinalityFeature(other_block_held),
        ConceptCardinalityFeature(holding)
    ]


def get_on_x_y_feature(lang):
    obj_t = lang.Object
    x_nominal = NominalConcept("a", obj_t)
    y_nominal = NominalConcept("b", obj_t)
    on_r = PrimitiveRole(lang.get("on"))
    on_x_y = AndConcept(ExistsConcept(on_r, y_nominal), x_nominal, "object")
    return [ConceptCardinalityFeature(on_x_y)]


def bwnamer(feature):
    s = str(feature)
    base = {
        "And(clear,Nominal(a))": "clear(a)",
        "And(clear,Nominal(b))": "clear(b)",
        "holding": "holding(·)",
        "And(Nominal(a),holding)": "holding(a)",
        "And(holding,Nominal(a))": "holding(a)",
        "And(holding,Nominal(b))": "holding(b)",
        "And(Exists(on,Nominal(b)),Nominal(a))": "on(a,b)",
        "And(Exists(Inverse(on),Nominal(a)),Nominal(b))": "on(a,b)",
        "And(Exists(Star(on),Nominal(b)),Nominal(a))": "above(a,b)",
        "And(Not(Nominal(a)),holding)": "H",
        "Exists(Inverse(on),Nominal(a))": "Z",
        "Exists(Star(on),Nominal(a))": "n(a)",
        "Exists(Star(on),Nominal(b))": "n(b)",
        "And(ontable,Nominal(a))": "ontable(a)",
        "And(Forall(on,Nominal(b)),Nominal(a))": "a_on_b_ontable_or_held",
        "And(And(And(Not(Exists(Star(on),Nominal(a))),Not(Exists(Star(Inverse(on)),Nominal(a)))),Not(Nominal(a))),Not(holding))": "m(a)",
        "And(And(Forall(Star(on),Not(Nominal(a))),Forall(Star(Inverse(on)),Not(Nominal(a)))),And(Not(holding),Not(Nominal(a))))": "m(a)",
        "Exists(Star(on),Exists(on,Nominal(b)))": "n-at-least-2-above-b",
    }
    return extend_namer_to_all_features(base).get(s, s)


def features_clear_x(lang):
    obj_t = lang.Object
    x_nominal = NominalConcept("a", obj_t)
    y_nominal = NominalConcept("b", obj_t)
    holding_p = lang.get_predicate("holding")
    on_p = lang.get_predicate("on")
    on_r = PrimitiveRole(on_p)
    on_star = StarRole(on_r)

    holding = PrimitiveConcept(holding_p)
    x_held = AndConcept(x_nominal, holding, "object")  # H: Block x is being held
    y_held = AndConcept(y_nominal, holding, "object")

    other_block = AndConcept(NotConcept(x_nominal, obj_t), NotConcept(y_nominal, obj_t), "object")
    other_block_held = AndConcept(other_block, holding, "object")

    blocks_above_x = ExistsConcept(on_star, x_nominal)
    blocks_above_y = ExistsConcept(on_star, y_nominal)

    on_x_y = AndConcept(ExistsConcept(on_r, y_nominal), x_nominal, "object")
    on_y_x = AndConcept(ExistsConcept(on_r, x_nominal), y_nominal, "object")

    return [
        ConceptCardinalityFeature(on_x_y),
        # ConceptCardinalityFeature(on_y_x),
        ConceptCardinalityFeature(blocks_above_x),
        ConceptCardinalityFeature(blocks_above_y),
        ConceptCardinalityFeature(x_held),
        # ConceptCardinalityFeature(y_held),
        # ConceptCardinalityFeature(other_block_held),
        ConceptCardinalityFeature(holding)
    ]


def add_bw_domain_parameters(language):
    # We simply add block "a" as a domain constant
    return [language.constant("a", "object")]
    # language.constant("b", "object")


def add_bw_domain_parameters_2(language):
    return [language.constant("a", "object"), language.constant("b", "object")]


def no_parameter(lang):
    return []


def ipc_instances():
    tpl = "probBLOCKS-{}-{}.pddl"
    return [tpl.format(b, i) for b, i in itertools.product(range(5, 11), range(0, 3))] +\
           [tpl.format(b, i) for b, i in itertools.product(range(11, 16), range(0, 2))]
