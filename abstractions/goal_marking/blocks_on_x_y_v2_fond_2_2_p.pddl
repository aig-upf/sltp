(define (problem FOND_blocks_2_2_p)
    (:domain FOND_blocks_2_2)
    (:init (bool[handempty]) (stack-depth d0) (bitvalue d0 b0) (bitvalue d0 b1) (bitvalue d0 b2))
    (:goal (and (not (holding_a)) (zero n_a) (not (zero n_b)) (on_a_b_2) (stack-depth d0)))
)

