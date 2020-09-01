from .misc import extend_namer_to_all_features


def no_parameter(lang):
    return []


def gripper_names(feature):
    s = str(feature)
    base = {
        "Exists(at,Not(Nominal(roomb)))": "nballs-A",
        "Exists(at,Nominal(roomb))": "nballs-B",
        "Exists(carry,<universe>)": "ncarried",
        "And(at-robby,Nominal(roomb))": "robot-at-B",
        "Exists(at,at-robby)": "nballs-in-room-with-robot",
        "Exists(at,Not(at-robby))": "nballs-in-rooms-with-no-robot",
        "free": "nfree-grippers",
        "Exists(at,Exists(Inverse(at-robby),<universe>))": "nballs-in-room-with-some-robot",
        "And(Exists(gripper,Exists(at-robby,{roomb})),free)": "nfree-grippers-at-B",
        "Exists(at-robby,{roomb})": "nrobots-at-B",
        "Exists(gripper,Exists(at-robby,{roomb}))": "ngrippers-at-B",
        "Exists(carry,Exists(gripper,Exists(at-robby,{roomb})))": "nballs-carried-in-B",
        "Exists(at,And(Forall(Inverse(at-robby),<empty>), Not({roomb})))":
            "nballs-in-some-room-notB-without-any-robot",
        "And(Exists(Inverse(at),<universe>), And({roomb}, Not(at-robby)))": "some-ball-in-B-but-robot-not-in-B",
        "And(Forall(Inverse(at),<empty>),room)": "num-empty-rooms",
        "Exists(at,And(at-robby,Nominal(roomb)))": "num-balls-at-B-when-robot-at-B-as-well",
        "Not(And(Forall(carry,<empty>),Forall(at,at-robby)))": "num-balls-either-carried-or-not-in-same-room-as-robot",
        # "Not(And(Not(And(at-robby,Nominal(roomb))),Forall(at,And(at-robby,Nominal(roomb)))))": "",
        # "Not(And(Not(And(Forall(at,at-robby),ball)),Not(And(at-robby,Nominal(roomb)))))": "",
        # "Not(And(Forall(at-robby,And(Not(Nominal(roomb)),Exists(Inverse(at),<universe>))),Forall(carry,<empty>)))":
        #     ""
        "And(Exists(carry,<universe>),Exists(at_g,at-robby))": "if-robot-at-B-then-num-carried-balls-else-emptyset",
    }
    return extend_namer_to_all_features(base).get(s, s)


def gripper_parameters(language):
    return [language.constant("roomb", "object")]


def spanner_names(feature):
    s = str(feature)
    base = {
        "And(tightened_g,Not(tightened))": "n-untightened-nuts",
        "Exists(Inverse(carrying),<universe>)": "n-carried-spanners",
        "Forall(Inverse(link),<empty>)": "first-cell",  # Neat!
        "Exists(at,Forall(Inverse(link),<empty>))": "n-things-on-first-cell",
        "And(Exists(at,Exists(Inverse(at),man)),Not(man))": "n-spanners-in-same-cell-as-man",
        "And(Exists(at,Exists(Inverse(at),man)),spanner)":  "n-spanners-in-same-cell-as-man",
        "Exists(at,Exists(link,Exists(Inverse(at),<universe>)))": "",
        "loose": "n-untightened-nuts",
        "Exists(at,Exists(link,Exists(Inverse(at),man)))": "n-spanners-on-cell-left-to-man",
        "Exists(Inverse(at),spanner)": "locations-with-a-spanner",
        "Exists(carrying,useable)": "bob-is-carrying-a-usable-spanner",
        "tightened": "num-tightened-nuts",
        "Exists(Star(link),Exists(Inverse(at),man))": "n-unreachable-locations",
        "Exists(at,Exists(Star(link),Exists(Inverse(at),man)))": "n-unreachable-spanners",
        "LessThan{Num[Exists(Inverse(carrying),<universe>)]}{Num[loose]}": "not-carrying-enough-spanners",
        "Exists(at,Forall(Inverse(at),man))": "bob-in-empty-loc",
        "Exists(Star(link),Exists(Inverse(at),spanner))": "num-locs-from-which-some-spanner-is-reachable",
        "Exists(Inverse(Star(link)),Exists(Inverse(at),<universe>))": "num-locs-reachable-from-a-loc-with-things",
        "And(Forall(at,<empty>),useable)": "picked-up-spanners",
    }
    d = extend_namer_to_all_features(base)
    return d.get(s, s)


def blocksworld_names(feature):
    s = str(feature)
    base = {
        "And(clear,Nominal(a))": "clear(a)",
        "And(clear,Nominal(b))": "clear(b)",
        "holding": "holding(·)",
        "And(Nominal(a),holding)": "holding(a)",
        "And(holding,Nominal(a))": "holding(a)",
        "And(holding,Nominal(b))": "holding(b)",
        "And(Exists(on,Nominal(b)),Nominal(a))": "on(a,b)",
        "And(Exists(loc,Nominal(b)),Nominal(a))": "on(a,b)",  # FSTRIPS
        "And(Exists(Inverse(on),Nominal(a)),Nominal(b))": "on(a,b)",
        "And(Exists(Star(on),Nominal(b)),Nominal(a))": "above(a,b)",
        "And(Not(Nominal(a)),holding)": "H",
        "Exists(Inverse(on),Nominal(a))": "Z",
        "Exists(Star(on),Nominal(a))": "n(a)",
        "Exists(Star(on),Nominal(b))": "n(b)",
        "Exists(Star(loc),Nominal(a))": "n(a)",  # FSTRIPS
        "Exists(Star(loc),Nominal(b))": "n(b)",  # FSTRIPS
        "And(ontable,Nominal(a))": "ontable(a)",
        "And(Forall(on,Nominal(b)),Nominal(a))": "a_on_b_ontable_or_held",
        "And(And(And(Not(Exists(Star(on),Nominal(a))),Not(Exists(Star(Inverse(on)),Nominal(a)))),Not(Nominal(a))),Not(holding))": "m(a)",
        "And(And(Forall(Star(on),Not(Nominal(a))),Forall(Star(Inverse(on)),Not(Nominal(a)))),And(Not(holding),Not(Nominal(a))))": "m(a)",
        "Exists(Star(on),Exists(on,Nominal(b)))": "n-at-least-2-above-b",
        "Not(clear)": "num-unclear",
        "And(clear,ontable)": "n-single-blocks",
        "Atom[handempty]": "handempty",
        # "superficially-well-placed": all blocks below are the same as in goal
        "And(Equal(Star(on_g),Star(on)),clear)": "n-clear-and-superficially-well-placed-blocks",
        "Equal(Inverse(loc_g),Inverse(Star(loc)))": "n-blocks-below-their-hat",  # FSTRIPS
        "Exists(Star(loc),Exists(loc_g,Not(Equal(Inverse(loc_g),Inverse(loc)))))": "n-x-with-misplaced-block-below",  # FSTRIPS
        "clear": "num-clear",
        "And(Equal(loc_g,loc),Forall(Star(loc),Equal(loc_g,loc)))": "n-well-placed",  # FSTRIPS
        "Equal(Star(loc_g),Star(loc))": "n-superficially-well-placed",  # FSTRIPS
        "Equal(loc_g,loc)": "n-ontarget",  # FSTRIPS
        "Equal(Inverse(loc_g),Inverse(loc))": "n-right-under-target",  # FSTRIPS
    }
    return extend_namer_to_all_features(base).get(s, s)


def blocksworld_parameters_for_clear(language):
    # We simply add block "a" as a domain constant
    return [language.constant("a", "object")]


def blocksworld_parameters_for_on(language):
    return [language.constant("a", "object"), language.constant("b", "object")]


def reward_names(feature):
    s = str(feature)
    base = {
        "reward": "num-rewards",
        "Dist[at, adjacent, reward]": "dist-to-closest-reward",
        "Dist[at, Restrict(adjacent,unblocked), reward]": "unblocked-dist-to-closest-reward",
    }

    return extend_namer_to_all_features(base).get(s, s)
