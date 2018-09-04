(define (problem FOND_gripper_0_2_p)
    (:domain FOND_gripper_0_2)
    (:init (zero ncarried) (stack-depth d0) (bitvalue d0 b0))
    (:goal (and (zero nballs-A) (not (zero nfree-grippers)) (zero ncarried) (stack-depth d0)))
)

