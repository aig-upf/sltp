(define (problem FOND_blocks_2_2_p)
    (:domain FOND_blocks_2_2)
    (:init (stack-depth d0) (bitvalue d0 b0) (bitvalue d0 b1) (bitvalue d0 b2))
    (:goal (and (not (holding_a)) (zero n_a) (stack-depth d0)))
)

