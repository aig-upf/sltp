(define (domain FOND_gripper_0_2)
    (:requirements :non-deterministic)
    (:types counter depth bit)
    (:constants nfree-grippers ncarried nballs-A - counter d0 d1 d2 - depth b0 - bit)
    (:predicates
        (zero ?c - counter)
        (stack-depth ?d - depth)
        (stack-idx ?c - counter ?d - depth)
        (in-stack ?c - counter)
        (bitvalue ?d - depth ?b - bit)
        (robot-at-B)
    )

    (:action Push_nfree-grippers_d0_b0
        :precondition (and (stack-depth d0) (bitvalue d0 b0))
        :effect (and (not (stack-depth d0)) (stack-depth d1) (in-stack nfree-grippers) (stack-idx nfree-grippers d1) (not (bitvalue d0 b0)) (bitvalue d1 b0))
    )
    (:action Push_ncarried_d0_b0
        :precondition (and (stack-depth d0) (not (in-stack ncarried)) (bitvalue d0 b0))
        :effect (and (not (stack-depth d0)) (stack-depth d1) (in-stack ncarried) (stack-idx ncarried d1) (not (bitvalue d0 b0)) (bitvalue d1 b0))
    )
    (:action Push_nballs-A_d0_b0
        :precondition (and (stack-depth d0) (not (in-stack nballs-A)) (bitvalue d0 b0))
        :effect (and (not (stack-depth d0)) (stack-depth d1) (in-stack nballs-A) (stack-idx nballs-A d1) (not (bitvalue d0 b0)) (bitvalue d1 b0))
    )
    (:action Push_nfree-grippers_d1_b0
        :precondition (and (stack-depth d1) (bitvalue d1 b0))
        :effect (and (not (stack-depth d1)) (stack-depth d2) (in-stack nfree-grippers) (stack-idx nfree-grippers d2) (not (bitvalue d1 b0)) (bitvalue d2 b0))
    )
    (:action Push_ncarried_d1_b0
        :precondition (and (stack-depth d1) (not (in-stack ncarried)) (bitvalue d1 b0))
        :effect (and (not (stack-depth d1)) (stack-depth d2) (in-stack ncarried) (stack-idx ncarried d2) (not (bitvalue d1 b0)) (bitvalue d2 b0))
    )
    (:action Push_nballs-A_d1_b0
        :precondition (and (stack-depth d1) (not (in-stack nballs-A)) (bitvalue d1 b0))
        :effect (and (not (stack-depth d1)) (stack-depth d2) (in-stack nballs-A) (stack-idx nballs-A d2) (not (bitvalue d1 b0)) (bitvalue d2 b0))
    )
    (:action Pop_nfree-grippers_d1
        :precondition (and (stack-depth d1) (in-stack nfree-grippers) (stack-idx nfree-grippers d1))
        :effect (and (stack-depth d0) (not (stack-depth d1)) (not (in-stack nfree-grippers)) (not (stack-idx nfree-grippers d1)))
    )
    (:action Pop_ncarried_d1
        :precondition (and (stack-depth d1) (in-stack ncarried) (stack-idx ncarried d1))
        :effect (and (stack-depth d0) (not (stack-depth d1)) (not (in-stack ncarried)) (not (stack-idx ncarried d1)))
    )
    (:action Pop_nballs-A_d1
        :precondition (and (stack-depth d1) (in-stack nballs-A) (stack-idx nballs-A d1))
        :effect (and (stack-depth d0) (not (stack-depth d1)) (not (in-stack nballs-A)) (not (stack-idx nballs-A d1)))
    )
    (:action Pop_nfree-grippers_d2
        :precondition (and (stack-depth d2) (in-stack nfree-grippers) (stack-idx nfree-grippers d2))
        :effect (and (stack-depth d1) (not (stack-depth d2)) (not (in-stack nfree-grippers)) (not (stack-idx nfree-grippers d2)))
    )
    (:action Pop_ncarried_d2
        :precondition (and (stack-depth d2) (in-stack ncarried) (stack-idx ncarried d2))
        :effect (and (stack-depth d1) (not (stack-depth d2)) (not (in-stack ncarried)) (not (stack-idx ncarried d2)))
    )
    (:action Pop_nballs-A_d2
        :precondition (and (stack-depth d2) (in-stack nballs-A) (stack-idx nballs-A d2))
        :effect (and (stack-depth d1) (not (stack-depth d2)) (not (in-stack nballs-A)) (not (stack-idx nballs-A d2)))
    )
    (:action action_1_d1
        :precondition (and (not (zero ncarried)) (robot-at-B) (in-stack ncarried) (not (in-stack nfree-grippers)) (stack-idx ncarried d1))
        :effect (and (oneof (zero ncarried) (not (zero ncarried))) (not (zero nfree-grippers)) (bitvalue d1 b0) (bitvalue d2 b0))
    )
    (:action action_1_d2
        :precondition (and (not (zero ncarried)) (robot-at-B) (in-stack ncarried) (not (in-stack nfree-grippers)) (stack-idx ncarried d2))
        :effect (and (oneof (zero ncarried) (not (zero ncarried))) (not (zero nfree-grippers)) (bitvalue d2 b0))
    )
    (:action action_2
        :precondition (and (zero nballs-A) (not (zero ncarried)) (not (zero nfree-grippers)) (not (robot-at-B)))
        :effect (and (robot-at-B))
    )
    (:action action_3
        :precondition (and (not (zero ncarried)) (zero nfree-grippers) (not (robot-at-B)))
        :effect (and (robot-at-B))
    )
    (:action action_4_0_d1
        :precondition (and (not (zero nballs-A)) (not (zero nfree-grippers)) (not (robot-at-B)) (in-stack nballs-A) (not (in-stack ncarried)) (stack-idx nballs-A d1))
        :effect (and (oneof (zero nballs-A) (not (zero nballs-A))) (not (zero ncarried)) (oneof (zero nfree-grippers) (not (zero nfree-grippers))) (bitvalue d1 b0) (bitvalue d2 b0))
    )
    (:action action_4_0_d2
        :precondition (and (not (zero nballs-A)) (not (zero nfree-grippers)) (not (robot-at-B)) (in-stack nballs-A) (not (in-stack ncarried)) (stack-idx nballs-A d2))
        :effect (and (oneof (zero nballs-A) (not (zero nballs-A))) (not (zero ncarried)) (oneof (zero nfree-grippers) (not (zero nfree-grippers))) (bitvalue d2 b0))
    )
    (:action action_4_1_d1
        :precondition (and (not (zero nballs-A)) (not (zero nfree-grippers)) (not (robot-at-B)) (in-stack nfree-grippers) (not (in-stack ncarried)) (stack-idx nfree-grippers d1))
        :effect (and (oneof (zero nballs-A) (not (zero nballs-A))) (not (zero ncarried)) (oneof (zero nfree-grippers) (not (zero nfree-grippers))) (bitvalue d1 b0) (bitvalue d2 b0))
    )
    (:action action_4_1_d2
        :precondition (and (not (zero nballs-A)) (not (zero nfree-grippers)) (not (robot-at-B)) (in-stack nfree-grippers) (not (in-stack ncarried)) (stack-idx nfree-grippers d2))
        :effect (and (oneof (zero nballs-A) (not (zero nballs-A))) (not (zero ncarried)) (oneof (zero nfree-grippers) (not (zero nfree-grippers))) (bitvalue d2 b0))
    )
    (:action action_5
        :precondition (and (zero ncarried) (not (zero nfree-grippers)) (robot-at-B))
        :effect (and (not (robot-at-B)))
    )
)

