(define (domain FOND_blocks-clear-x_3)
    (:types counter)
    (:constants nabove_A nother_A)
    (:predicates
        (zero ?c - counter)
        (hold_A)
        (hold-other_A)
        (some-below_A)
        (Q_nabove_A)
        (Q_nother_A)
        (bit_0)
        (bit_1)
        (bit_2)
        (bit_3)
    )

    (:action Set_nabove_A,0
        (:precondition (and (not (Q_nabove_A)) (not (Q_nother_A)) (bit_0)))
        (:effect (and (Q_nabove_A) (not (bit_0))))
    )
    (:action Set_nabove_A,1
        (:precondition (and (not (Q_nabove_A)) (not (Q_nother_A)) (bit_1) (not (bit_0))))
        (:effect (and (Q_nabove_A) (not (bit_1)) (bit_0)))
    )
    (:action Set_nabove_A,2
        (:precondition (and (not (Q_nabove_A)) (not (Q_nother_A)) (bit_2) (not (bit_1)) (not (bit_0))))
        (:effect (and (Q_nabove_A) (not (bit_2)) (bit_1) (bit_0)))
    )
    (:action Set_nabove_A,3
        (:precondition (and (not (Q_nabove_A)) (not (Q_nother_A)) (bit_3) (not (bit_2)) (not (bit_1)) (not (bit_0))))
        (:effect (and (Q_nabove_A) (not (bit_3)) (bit_2) (bit_1) (bit_0)))
    )
    (:action Set_nother_A,0
        (:precondition (and (not (Q_nabove_A)) (not (Q_nother_A)) (bit_0)))
        (:effect (and (Q_nother_A) (not (bit_0))))
    )
    (:action Set_nother_A,1
        (:precondition (and (not (Q_nabove_A)) (not (Q_nother_A)) (bit_1) (not (bit_0))))
        (:effect (and (Q_nother_A) (not (bit_1)) (bit_0)))
    )
    (:action Set_nother_A,2
        (:precondition (and (not (Q_nabove_A)) (not (Q_nother_A)) (bit_2) (not (bit_1)) (not (bit_0))))
        (:effect (and (Q_nother_A) (not (bit_2)) (bit_1) (bit_0)))
    )
    (:action Set_nother_A,3
        (:precondition (and (not (Q_nabove_A)) (not (Q_nother_A)) (bit_3) (not (bit_2)) (not (bit_1)) (not (bit_0))))
        (:effect (and (Q_nother_A) (not (bit_3)) (bit_2) (bit_1) (bit_0)))
    )
    (:action Unset_nabove_A
        (:precondition (and (Q_nabove_A)))
        (:effect (and (not (Q_nabove_A))))
    )
    (:action Unset_nother_A
        (:precondition (and (Q_nother_A)))
        (:effect (and (not (Q_nother_A))))
    )
    (:action Pick-x-some-below
        (:precondition (and (zero (nabove_A)) (not (hold_A)) (not (hold-other_A)) (some-below_A)))
        (:effect (and (not (zero (nother_A))) (hold_A) (not (some-below_A))))
    )
    (:action Pick-x-none-below
        (:precondition (and (zero (nabove_A)) (not (hold_A)) (not (hold-other_A)) (not (some-below_A))))
        (:effect (and (hold_A)))
    )
    (:action Pick-above-x
        (:precondition (and (not (zero (nabove_A))) (not (hold_A)) (not (hold-other_A)) (Q_nabove_A)))
        (:effect (and (oneof (zero (nabove_A)) (not (zero (nabove_A)))) (hold-other_A)))
    )
    (:action Pick-other
        (:precondition (and (not (zero (nother_A))) (not (hold_A)) (not (hold-other_A)) (Q_nabove_A)))
        (:effect (and (oneof (zero (nother_A)) (not (zero (nother_A)))) (hold-other_A)))
    )
    (:action Put-x-on-table
        (:precondition (and (hold_A)))
        (:effect (and (not (hold_A))))
    )
    (:action Put-x-above-some
        (:precondition (and (not (zero (nother_A))) (hold_A) (Q_nabove_A)))
        (:effect (and (oneof (zero (nother_A)) (not (zero (nother_A)))) (not (hold_A)) (some-below_A)))
    )
    (:action Put-aside
        (:precondition (and (hold-other_A)))
        (:effect (and (not (zero (nother_A))) (not (hold-other_A))))
    )
    (:action Put-above-x
        (:precondition (and (hold-other_A)))
        (:effect (and (not (zero (nabove_A))) (not (hold-other_A))))
    )
)

