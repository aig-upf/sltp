(define (problem FOND_gripper_2_2_p)
    (:domain FOND_gripper_2_2)
    (:init (zero ncarried) (stack-depth d0) (bitvalue d0 b0) (bitvalue d0 b1) (bitvalue d0 b2))
    (:goal (and (zero nballs-A) (not (zero nfree-grippers)) (zero ncarried) (stack-depth d0)))
)

