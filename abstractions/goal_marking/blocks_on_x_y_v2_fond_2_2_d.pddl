(define (domain FOND_blocks_2_2)
    (:requirements :non-deterministic)
    (:types counter depth bit)
    (:constants n_a n_b - counter d0 d1 d2 - depth b0 b1 b2 - bit)
    (:predicates
        (zero ?c - counter)
        (stack-depth ?d - depth)
        (stack-idx ?c - counter ?d - depth)
        (in-stack ?c - counter)
        (bitvalue ?d - depth ?b - bit)
        (bool[handempty])
        (holding_a)
        (on_a_b_2)
    )

    (:action PUSH_d0_b0
        :parameters (?c - counter)
        :precondition (and (stack-depth d0) (bitvalue d0 b0))
        :effect (and (not (stack-depth d0)) (stack-depth d1) (in-stack ?c) (stack-idx ?c d1) (not (bitvalue d0 b0)) (bitvalue d1 b2) (bitvalue d1 b1) (bitvalue d1 b0))
    )
    (:action PUSH_d0_b1
        :parameters (?c - counter)
        :precondition (and (stack-depth d0) (bitvalue d0 b1) (not (bitvalue d0 b0)))
        :effect (and (not (stack-depth d0)) (stack-depth d1) (in-stack ?c) (stack-idx ?c d1) (not (bitvalue d0 b1)) (bitvalue d0 b0) (bitvalue d1 b2) (bitvalue d1 b1) (bitvalue d1 b0))
    )
    (:action PUSH_d0_b2
        :parameters (?c - counter)
        :precondition (and (stack-depth d0) (bitvalue d0 b2) (not (bitvalue d0 b1)) (not (bitvalue d0 b0)))
        :effect (and (not (stack-depth d0)) (stack-depth d1) (in-stack ?c) (stack-idx ?c d1) (not (bitvalue d0 b2)) (bitvalue d0 b1) (bitvalue d0 b0) (bitvalue d1 b2) (bitvalue d1 b1) (bitvalue d1 b0))
    )
    (:action PUSH_d1_b0
        :parameters (?c - counter)
        :precondition (and (stack-depth d1) (not (in-stack ?c)) (bitvalue d1 b0))
        :effect (and (not (stack-depth d1)) (stack-depth d2) (in-stack ?c) (stack-idx ?c d2) (not (bitvalue d1 b0)) (bitvalue d2 b2) (bitvalue d2 b1) (bitvalue d2 b0))
    )
    (:action PUSH_d1_b1
        :parameters (?c - counter)
        :precondition (and (stack-depth d1) (not (in-stack ?c)) (bitvalue d1 b1) (not (bitvalue d1 b0)))
        :effect (and (not (stack-depth d1)) (stack-depth d2) (in-stack ?c) (stack-idx ?c d2) (not (bitvalue d1 b1)) (bitvalue d1 b0) (bitvalue d2 b2) (bitvalue d2 b1) (bitvalue d2 b0))
    )
    (:action PUSH_d1_b2
        :parameters (?c - counter)
        :precondition (and (stack-depth d1) (not (in-stack ?c)) (bitvalue d1 b2) (not (bitvalue d1 b1)) (not (bitvalue d1 b0)))
        :effect (and (not (stack-depth d1)) (stack-depth d2) (in-stack ?c) (stack-idx ?c d2) (not (bitvalue d1 b2)) (bitvalue d1 b1) (bitvalue d1 b0) (bitvalue d2 b2) (bitvalue d2 b1) (bitvalue d2 b0))
    )
    (:action POP_d1
        :parameters (?c - counter)
        :precondition (and (stack-depth d1) (in-stack ?c) (stack-idx ?c d1))
        :effect (and (not (stack-depth d1)) (stack-depth d0) (not (in-stack ?c)) (not (stack-idx ?c d1)))
    )
    (:action POP_d2
        :parameters (?c - counter)
        :precondition (and (stack-depth d2) (in-stack ?c) (stack-idx ?c d2))
        :effect (and (not (stack-depth d2)) (stack-depth d1) (not (in-stack ?c)) (not (stack-idx ?c d2)))
    )
    (:action action_1
        :precondition (and (not (bool[handempty])) (not (holding_a)) (zero n_a) (not (on_a_b_2)))
        :effect (and (bool[handempty]))
    )
    (:action action_2
        :precondition (and (not (bool[handempty])) (not (holding_a)) (not (zero n_a)) (not (zero n_b)) (not (on_a_b_2)))
        :effect (and (bool[handempty]))
    )
    (:action action_3
        :precondition (and (bool[handempty]) (not (holding_a)) (zero n_a) (zero n_b) (not (on_a_b_2)))
        :effect (and (not (bool[handempty])) (holding_a))
    )
    (:action action_4
        :precondition (and (not (bool[handempty])) (holding_a) (zero n_a) (not (zero n_b)) (not (on_a_b_2)))
        :effect (and (bool[handempty]) (not (holding_a)))
    )
    (:action action_5_d1
        :precondition (and (bool[handempty]) (not (holding_a)) (not (zero n_a)) (not (zero n_b)) (not (on_a_b_2)) (in-stack n_a) (stack-idx n_a d1))
        :effect (and (not (bool[handempty])) (oneof (zero n_a) (not (zero n_a))) (bitvalue d1 b0) (bitvalue d1 b1) (bitvalue d1 b2) (bitvalue d2 b0) (bitvalue d2 b1) (bitvalue d2 b2))
    )
    (:action action_5_d2
        :precondition (and (bool[handempty]) (not (holding_a)) (not (zero n_a)) (not (zero n_b)) (not (on_a_b_2)) (in-stack n_a) (stack-idx n_a d2))
        :effect (and (not (bool[handempty])) (oneof (zero n_a) (not (zero n_a))) (bitvalue d2 b0) (bitvalue d2 b1) (bitvalue d2 b2))
    )
    (:action action_6_d1
        :precondition (and (bool[handempty]) (not (holding_a)) (zero n_a) (not (zero n_b)) (not (on_a_b_2)) (in-stack n_b) (stack-idx n_b d1))
        :effect (and (not (bool[handempty])) (oneof (zero n_b) (not (zero n_b))) (bitvalue d1 b0) (bitvalue d1 b1) (bitvalue d1 b2) (bitvalue d2 b0) (bitvalue d2 b1) (bitvalue d2 b2))
    )
    (:action action_6_d2
        :precondition (and (bool[handempty]) (not (holding_a)) (zero n_a) (not (zero n_b)) (not (on_a_b_2)) (in-stack n_b) (stack-idx n_b d2))
        :effect (and (not (bool[handempty])) (oneof (zero n_b) (not (zero n_b))) (bitvalue d2 b0) (bitvalue d2 b1) (bitvalue d2 b2))
    )
    (:action action_7
        :precondition (and (not (bool[handempty])) (holding_a) (zero n_a) (zero n_b) (not (on_a_b_2)) (not (in-stack n_b)))
        :effect (and (bool[handempty]) (not (holding_a)) (not (zero n_b)) (on_a_b_2))
    )
)

