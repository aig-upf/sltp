(define (domain FOND_blocks_1_0)
    (:requirements :non-deterministic)
    (:types counter bit)
    (:constants n_a n_b - counter b0 b1 - bit)
    (:predicates
        (zero ?c - counter)
        (q ?c - counter)
        (bitvalue ?b - bit)
        (bool[holding])
        (ontable_a)
        (holding_a)
        (a_on_b_ontable_or_held)
    )

    (:action Set_n_a_b0
        :precondition (and (not (q n_a)) (not (q n_b)) (bitvalue b0))
        :effect (and (q n_a) (not (bitvalue b0)))
    )
    (:action Set_n_a_b1
        :precondition (and (not (q n_a)) (not (q n_b)) (bitvalue b1) (not (bitvalue b0)))
        :effect (and (q n_a) (not (bitvalue b1)) (bitvalue b0))
    )
    (:action Set_n_b_b0
        :precondition (and (not (q n_a)) (not (q n_b)) (bitvalue b0))
        :effect (and (q n_b) (not (bitvalue b0)))
    )
    (:action Set_n_b_b1
        :precondition (and (not (q n_a)) (not (q n_b)) (bitvalue b1) (not (bitvalue b0)))
        :effect (and (q n_b) (not (bitvalue b1)) (bitvalue b0))
    )
    (:action Unset_n_a
        :precondition (and (q n_a))
        :effect (and (not (q n_a)))
    )
    (:action Unset_n_b
        :precondition (and (q n_b))
        :effect (and (not (q n_b)))
    )
    (:action action_1
        :precondition (and (a_on_b_ontable_or_held) (not (bool[holding])) (not (holding_a)) (zero n_a) (not (zero n_b)) (ontable_a) (q n_b))
        :effect (and (bool[holding]) (oneof (zero n_b) (not (zero n_b))))
    )
    (:action action_2
        :precondition (and (a_on_b_ontable_or_held) (not (bool[holding])) (not (holding_a)) (zero n_b) (ontable_a))
        :effect (and (bool[holding]))
    )
    (:action action_3
        :precondition (and (a_on_b_ontable_or_held) (not (bool[holding])) (not (holding_a)) (zero n_a) (not (zero n_b)) (not (ontable_a)))
        :effect (and (bool[holding]))
    )
    (:action action_4
        :precondition (and (a_on_b_ontable_or_held) (bool[holding]) (not (holding_a)) (zero n_a) (not (zero n_b)) (not (ontable_a)) (not (q n_a)) (not (q n_b)))
        :effect (and (not (bool[holding])) (not (zero n_a)) (not (zero n_b)))
    )
    (:action action_5
        :precondition (and (a_on_b_ontable_or_held) (not (bool[holding])) (not (holding_a)) (zero n_a) (zero n_b) (ontable_a))
        :effect (and (bool[holding]) (holding_a) (not (ontable_a)))
    )
    (:action action_6
        :precondition (and (a_on_b_ontable_or_held) (not (bool[holding])) (not (holding_a)) (not (zero n_a)) (zero n_b) (ontable_a) (q n_a))
        :effect (and (bool[holding]) (oneof (zero n_a) (not (zero n_a))))
    )
    (:action action_7
        :precondition (and (a_on_b_ontable_or_held) (bool[holding]) (not (holding_a)) (zero n_a) (not (zero n_b)) (not (ontable_a)))
        :effect (and (not (bool[holding])))
    )
    (:action action_8
        :precondition (and (a_on_b_ontable_or_held) (bool[holding]) (not (holding_a)) (zero n_a) (ontable_a))
        :effect (and (not (bool[holding])))
    )
    (:action action_9
        :precondition (and (a_on_b_ontable_or_held) (bool[holding]) (not (holding_a)) (not (zero n_a)) (zero n_b) (ontable_a))
        :effect (and (not (bool[holding])))
    )
    (:action action_10
        :precondition (and (a_on_b_ontable_or_held) (bool[holding]) (holding_a) (zero n_a) (not (zero n_b)) (not (ontable_a)))
        :effect (and (not (bool[holding])) (not (holding_a)) (ontable_a))
    )
    (:action action_11
        :precondition (and (not (a_on_b_ontable_or_held)) (not (bool[holding])) (not (holding_a)) (zero n_a) (zero n_b) (not (ontable_a)))
        :effect (and (a_on_b_ontable_or_held) (bool[holding]) (holding_a))
    )
    (:action action_12
        :precondition (and (a_on_b_ontable_or_held) (bool[holding]) (holding_a) (zero n_a) (zero n_b) (not (ontable_a)) (not (q n_b)))
        :effect (and (not (bool[holding])) (not (holding_a)) (not (zero n_b)))
    )
)

