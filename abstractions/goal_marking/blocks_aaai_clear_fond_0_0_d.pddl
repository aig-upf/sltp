(define (domain FOND_blocks_0_0)
    (:requirements :non-deterministic)
    (:types counter bit)
    (:constants n_a - counter b0 - bit)
    (:predicates
        (zero ?c - counter)
        (q ?c - counter)
        (bitvalue ?b - bit)
        (bool[holding])
        (holding_a)
    )

    (:action Set_n_a_b0
        :precondition (and (not (q n_a)) (bitvalue b0))
        :effect (and (q n_a) (not (bitvalue b0)))
    )
    (:action Unset_n_a
        :precondition (and (q n_a))
        :effect (and (not (q n_a)))
    )
    (:action action_1
        :precondition (and (bool[holding]) (not (holding_a)))
        :effect (and (not (bool[holding])))
    )
    (:action action_2
        :precondition (and (not (bool[holding])) (not (holding_a)) (not (zero n_a)) (q n_a))
        :effect (and (bool[holding]) (oneof (zero n_a) (not (zero n_a))))
    )
)

