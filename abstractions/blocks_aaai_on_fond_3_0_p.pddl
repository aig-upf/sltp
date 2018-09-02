(define (problem FOND_blocks_3_0_p)
    (:domain FOND_blocks_3_0)
    (:init (zero n_a) (holding_a) (bool[holding]) (a_on_b_ontable_or_held) (bitvalue b0) (bitvalue b1) (bitvalue b2) (bitvalue b3))
    (:goal (and (not (zero n_b)) (not (ontable_a)) (a_on_b_ontable_or_held) (not (holding_a))))
)

