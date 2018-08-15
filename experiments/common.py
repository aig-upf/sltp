from tarski.dl import NominalConcept, PrimitiveRole, InverseRole, StarRole, PrimitiveConcept, NotConcept, AndConcept, \
    ExistsConcept


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


def ijcai_paper_bw_feature_namer(feature):
    s = str(feature)
    return {
        "bool[And({a}, holding)]": "X",
        "bool[And(Not({a}), holding)]": "H",
        "bool[Exists(Inverse(on),{a})]": "Z",
        "card[Exists(Star(on),{a})]": "n(x)",
        "card[And(And(And(Not(Exists(Star(on),{a})), Not(Exists(Star(Inverse(on)),{a}))), Not({a})), Not(holding))]": "m(x)",
        "card[And(And(Forall(Star(on),Not({a})), Forall(Star(Inverse(on)),Not({a}))), And(Not(holding), Not({a})))]": "m(x)",
    }.get(s, s)


def add_bw_domain_parameters(language):
    # We simply add block "a" as a domain constant
    return [language.constant("a", "object")]
    # language.constant("b", "object")
