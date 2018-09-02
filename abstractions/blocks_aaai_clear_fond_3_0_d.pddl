(define (domain FOND_blocks_3_0)
    (:requirements :non-deterministic)
    (:types counter bit)
    (:constants n_a - counter b0 b1 b2 b3 - bit)
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
    (:action Set_n_a_b1
        :precondition (and (not (q n_a)) (bitvalue b1) (not (bitvalue b0)))
        :effect (and (q n_a) (not (bitvalue b1)) (bitvalue b0))
    )
    (:action Set_n_a_b2
        :precondition (and (not (q n_a)) (bitvalue b2) (not (bitvalue b1)) (not (bitvalue b0)))
        :effect (and (q n_a) (not (bitvalue b2)) (bitvalue b1) (bitvalue b0))
    )
    (:action Set_n_a_b3
        :precondition (and (not (q n_a)) (bitvalue b3) (not (bitvalue b2)) (not (bitvalue b1)) (not (bitvalue b0)))
        :effect (and (q n_a) (not (bitvalue b3)) (bitvalue b2) (bitvalue b1) (bitvalue b0))
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

