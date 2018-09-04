(define (domain FOND_blocks_2_2)
    (:requirements :non-deterministic)
    (:types counter depth bit)
    (:constants n_a - counter d0 d1 d2 - depth b0 b1 b2 - bit)
    (:predicates
        (zero ?c - counter)
        (stack-depth ?d - depth)
        (stack-idx ?c - counter ?d - depth)
        (in-stack ?c - counter)
        (bitvalue ?d - depth ?b - bit)
        (bool[holding])
        (holding_a)
    )

    (:action Push_n_a_d0_b0
        :precondition (and (stack-depth d0) (bitvalue d0 b0))
        :effect (and (not (stack-depth d0)) (stack-depth d1) (in-stack n_a) (stack-idx n_a d1) (not (bitvalue d0 b0)) (bitvalue d1 b0) (bitvalue d1 b1) (bitvalue d1 b2))
    )
    (:action Push_n_a_d0_b1
        :precondition (and (stack-depth d0) (bitvalue d0 b1) (not (bitvalue d0 b0)))
        :effect (and (not (stack-depth d0)) (stack-depth d1) (in-stack n_a) (stack-idx n_a d1) (not (bitvalue d0 b1)) (bitvalue d0 b0) (bitvalue d1 b0) (bitvalue d1 b1) (bitvalue d1 b2))
    )
    (:action Push_n_a_d0_b2
        :precondition (and (stack-depth d0) (bitvalue d0 b2) (not (bitvalue d0 b1)) (not (bitvalue d0 b0)))
        :effect (and (not (stack-depth d0)) (stack-depth d1) (in-stack n_a) (stack-idx n_a d1) (not (bitvalue d0 b2)) (bitvalue d0 b1) (bitvalue d0 b0) (bitvalue d1 b0) (bitvalue d1 b1) (bitvalue d1 b2))
    )
    (:action Push_n_a_d1_b0
        :precondition (and (stack-depth d1) (bitvalue d1 b0))
        :effect (and (not (stack-depth d1)) (stack-depth d2) (in-stack n_a) (stack-idx n_a d2) (not (bitvalue d1 b0)) (bitvalue d2 b0) (bitvalue d2 b1) (bitvalue d2 b2))
    )
    (:action Push_n_a_d1_b1
        :precondition (and (stack-depth d1) (bitvalue d1 b1) (not (bitvalue d1 b0)))
        :effect (and (not (stack-depth d1)) (stack-depth d2) (in-stack n_a) (stack-idx n_a d2) (not (bitvalue d1 b1)) (bitvalue d1 b0) (bitvalue d2 b0) (bitvalue d2 b1) (bitvalue d2 b2))
    )
    (:action Push_n_a_d1_b2
        :precondition (and (stack-depth d1) (bitvalue d1 b2) (not (bitvalue d1 b1)) (not (bitvalue d1 b0)))
        :effect (and (not (stack-depth d1)) (stack-depth d2) (in-stack n_a) (stack-idx n_a d2) (not (bitvalue d1 b2)) (bitvalue d1 b1) (bitvalue d1 b0) (bitvalue d2 b0) (bitvalue d2 b1) (bitvalue d2 b2))
    )
    (:action Pop_n_a_d1
        :precondition (and (stack-depth d1) (in-stack n_a) (stack-idx n_a d1))
        :effect (and (stack-depth d0) (not (stack-depth d1)) (not (in-stack n_a)) (not (stack-idx n_a d1)))
    )
    (:action Pop_n_a_d2
        :precondition (and (stack-depth d2) (in-stack n_a) (stack-idx n_a d2))
        :effect (and (stack-depth d1) (not (stack-depth d2)) (not (in-stack n_a)) (not (stack-idx n_a d2)))
    )
    (:action action_1
        :precondition (and (bool[holding]) (not (holding_a)))
        :effect (and (not (bool[holding])))
    )
    (:action action_2_d1
        :precondition (and (not (bool[holding])) (not (holding_a)) (not (zero n_a)) (in-stack n_a) (stack-idx n_a d1))
        :effect (and (bool[holding]) (oneof (zero n_a) (not (zero n_a))) (bitvalue d1 b0) (bitvalue d1 b1) (bitvalue d1 b2) (bitvalue d2 b0) (bitvalue d2 b1) (bitvalue d2 b2))
    )
    (:action action_2_d2
        :precondition (and (not (bool[holding])) (not (holding_a)) (not (zero n_a)) (in-stack n_a) (stack-idx n_a d2))
        :effect (and (bool[holding]) (oneof (zero n_a) (not (zero n_a))) (bitvalue d2 b0) (bitvalue d2 b1) (bitvalue d2 b2))
    )
)

