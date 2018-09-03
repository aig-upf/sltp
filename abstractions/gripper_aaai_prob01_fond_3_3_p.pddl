(define (problem FOND_gripper_3_3_p)
    (:domain FOND_gripper_3_3)
    (:init (zero ncarried) (stack-depth d0) (bitvalue d0 b0) (bitvalue d0 b1) (bitvalue d0 b2) (bitvalue d0 b3))
    (:goal (and (zero nballs-A) (not (zero nfree-grippers)) (zero ncarried) (stack-depth d0)))
)

