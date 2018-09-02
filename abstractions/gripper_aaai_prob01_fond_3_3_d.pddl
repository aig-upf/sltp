(define (domain FOND_gripper_3_3)
    (:requirements :non-deterministic)
    (:types counter depth bit)
    (:constants nfree-grippers ncarried nballs-A - counter d0 d1 d2 d3 - depth b0 b1 b2 b3 - bit)
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
        :effect (and (not (stack-depth d0)) (stack-depth d1) (in-stack nfree-grippers) (stack-idx nfree-grippers d1) (not (bitvalue d0 b0)) (bitvalue d1 b0) (bitvalue d1 b1) (bitvalue d1 b2) (bitvalue d1 b3))
    )
    (:action Push_nfree-grippers_d0_b1
        :precondition (and (stack-depth d0) (bitvalue d0 b1) (not (bitvalue d0 b0)))
        :effect (and (not (stack-depth d0)) (stack-depth d1) (in-stack nfree-grippers) (stack-idx nfree-grippers d1) (not (bitvalue d0 b1)) (bitvalue d0 b0) (bitvalue d1 b0) (bitvalue d1 b1) (bitvalue d1 b2) (bitvalue d1 b3))
    )
    (:action Push_nfree-grippers_d0_b2
        :precondition (and (stack-depth d0) (bitvalue d0 b2) (not (bitvalue d0 b1)) (not (bitvalue d0 b0)))
        :effect (and (not (stack-depth d0)) (stack-depth d1) (in-stack nfree-grippers) (stack-idx nfree-grippers d1) (not (bitvalue d0 b2)) (bitvalue d0 b1) (bitvalue d0 b0) (bitvalue d1 b0) (bitvalue d1 b1) (bitvalue d1 b2) (bitvalue d1 b3))
    )
    (:action Push_nfree-grippers_d0_b3
        :precondition (and (stack-depth d0) (bitvalue d0 b3) (not (bitvalue d0 b2)) (not (bitvalue d0 b1)) (not (bitvalue d0 b0)))
        :effect (and (not (stack-depth d0)) (stack-depth d1) (in-stack nfree-grippers) (stack-idx nfree-grippers d1) (not (bitvalue d0 b3)) (bitvalue d0 b2) (bitvalue d0 b1) (bitvalue d0 b0) (bitvalue d1 b0) (bitvalue d1 b1) (bitvalue d1 b2) (bitvalue d1 b3))
    )
    (:action Push_ncarried_d0_b0
        :precondition (and (stack-depth d0) (not (in-stack ncarried)) (bitvalue d0 b0))
        :effect (and (not (stack-depth d0)) (stack-depth d1) (in-stack ncarried) (stack-idx ncarried d1) (not (bitvalue d0 b0)) (bitvalue d1 b0) (bitvalue d1 b1) (bitvalue d1 b2) (bitvalue d1 b3))
    )
    (:action Push_ncarried_d0_b1
        :precondition (and (stack-depth d0) (not (in-stack ncarried)) (bitvalue d0 b1) (not (bitvalue d0 b0)))
        :effect (and (not (stack-depth d0)) (stack-depth d1) (in-stack ncarried) (stack-idx ncarried d1) (not (bitvalue d0 b1)) (bitvalue d0 b0) (bitvalue d1 b0) (bitvalue d1 b1) (bitvalue d1 b2) (bitvalue d1 b3))
    )
    (:action Push_ncarried_d0_b2
        :precondition (and (stack-depth d0) (not (in-stack ncarried)) (bitvalue d0 b2) (not (bitvalue d0 b1)) (not (bitvalue d0 b0)))
        :effect (and (not (stack-depth d0)) (stack-depth d1) (in-stack ncarried) (stack-idx ncarried d1) (not (bitvalue d0 b2)) (bitvalue d0 b1) (bitvalue d0 b0) (bitvalue d1 b0) (bitvalue d1 b1) (bitvalue d1 b2) (bitvalue d1 b3))
    )
    (:action Push_ncarried_d0_b3
        :precondition (and (stack-depth d0) (not (in-stack ncarried)) (bitvalue d0 b3) (not (bitvalue d0 b2)) (not (bitvalue d0 b1)) (not (bitvalue d0 b0)))
        :effect (and (not (stack-depth d0)) (stack-depth d1) (in-stack ncarried) (stack-idx ncarried d1) (not (bitvalue d0 b3)) (bitvalue d0 b2) (bitvalue d0 b1) (bitvalue d0 b0) (bitvalue d1 b0) (bitvalue d1 b1) (bitvalue d1 b2) (bitvalue d1 b3))
    )
    (:action Push_nballs-A_d0_b0
        :precondition (and (stack-depth d0) (not (in-stack nballs-A)) (bitvalue d0 b0))
        :effect (and (not (stack-depth d0)) (stack-depth d1) (in-stack nballs-A) (stack-idx nballs-A d1) (not (bitvalue d0 b0)) (bitvalue d1 b0) (bitvalue d1 b1) (bitvalue d1 b2) (bitvalue d1 b3))
    )
    (:action Push_nballs-A_d0_b1
        :precondition (and (stack-depth d0) (not (in-stack nballs-A)) (bitvalue d0 b1) (not (bitvalue d0 b0)))
        :effect (and (not (stack-depth d0)) (stack-depth d1) (in-stack nballs-A) (stack-idx nballs-A d1) (not (bitvalue d0 b1)) (bitvalue d0 b0) (bitvalue d1 b0) (bitvalue d1 b1) (bitvalue d1 b2) (bitvalue d1 b3))
    )
    (:action Push_nballs-A_d0_b2
        :precondition (and (stack-depth d0) (not (in-stack nballs-A)) (bitvalue d0 b2) (not (bitvalue d0 b1)) (not (bitvalue d0 b0)))
        :effect (and (not (stack-depth d0)) (stack-depth d1) (in-stack nballs-A) (stack-idx nballs-A d1) (not (bitvalue d0 b2)) (bitvalue d0 b1) (bitvalue d0 b0) (bitvalue d1 b0) (bitvalue d1 b1) (bitvalue d1 b2) (bitvalue d1 b3))
    )
    (:action Push_nballs-A_d0_b3
        :precondition (and (stack-depth d0) (not (in-stack nballs-A)) (bitvalue d0 b3) (not (bitvalue d0 b2)) (not (bitvalue d0 b1)) (not (bitvalue d0 b0)))
        :effect (and (not (stack-depth d0)) (stack-depth d1) (in-stack nballs-A) (stack-idx nballs-A d1) (not (bitvalue d0 b3)) (bitvalue d0 b2) (bitvalue d0 b1) (bitvalue d0 b0) (bitvalue d1 b0) (bitvalue d1 b1) (bitvalue d1 b2) (bitvalue d1 b3))
    )
    (:action Push_nfree-grippers_d1_b0
        :precondition (and (stack-depth d1) (bitvalue d1 b0))
        :effect (and (not (stack-depth d1)) (stack-depth d2) (in-stack nfree-grippers) (stack-idx nfree-grippers d2) (not (bitvalue d1 b0)) (bitvalue d2 b0) (bitvalue d2 b1) (bitvalue d2 b2) (bitvalue d2 b3))
    )
    (:action Push_nfree-grippers_d1_b1
        :precondition (and (stack-depth d1) (bitvalue d1 b1) (not (bitvalue d1 b0)))
        :effect (and (not (stack-depth d1)) (stack-depth d2) (in-stack nfree-grippers) (stack-idx nfree-grippers d2) (not (bitvalue d1 b1)) (bitvalue d1 b0) (bitvalue d2 b0) (bitvalue d2 b1) (bitvalue d2 b2) (bitvalue d2 b3))
    )
    (:action Push_nfree-grippers_d1_b2
        :precondition (and (stack-depth d1) (bitvalue d1 b2) (not (bitvalue d1 b1)) (not (bitvalue d1 b0)))
        :effect (and (not (stack-depth d1)) (stack-depth d2) (in-stack nfree-grippers) (stack-idx nfree-grippers d2) (not (bitvalue d1 b2)) (bitvalue d1 b1) (bitvalue d1 b0) (bitvalue d2 b0) (bitvalue d2 b1) (bitvalue d2 b2) (bitvalue d2 b3))
    )
    (:action Push_nfree-grippers_d1_b3
        :precondition (and (stack-depth d1) (bitvalue d1 b3) (not (bitvalue d1 b2)) (not (bitvalue d1 b1)) (not (bitvalue d1 b0)))
        :effect (and (not (stack-depth d1)) (stack-depth d2) (in-stack nfree-grippers) (stack-idx nfree-grippers d2) (not (bitvalue d1 b3)) (bitvalue d1 b2) (bitvalue d1 b1) (bitvalue d1 b0) (bitvalue d2 b0) (bitvalue d2 b1) (bitvalue d2 b2) (bitvalue d2 b3))
    )
    (:action Push_ncarried_d1_b0
        :precondition (and (stack-depth d1) (not (in-stack ncarried)) (bitvalue d1 b0))
        :effect (and (not (stack-depth d1)) (stack-depth d2) (in-stack ncarried) (stack-idx ncarried d2) (not (bitvalue d1 b0)) (bitvalue d2 b0) (bitvalue d2 b1) (bitvalue d2 b2) (bitvalue d2 b3))
    )
    (:action Push_ncarried_d1_b1
        :precondition (and (stack-depth d1) (not (in-stack ncarried)) (bitvalue d1 b1) (not (bitvalue d1 b0)))
        :effect (and (not (stack-depth d1)) (stack-depth d2) (in-stack ncarried) (stack-idx ncarried d2) (not (bitvalue d1 b1)) (bitvalue d1 b0) (bitvalue d2 b0) (bitvalue d2 b1) (bitvalue d2 b2) (bitvalue d2 b3))
    )
    (:action Push_ncarried_d1_b2
        :precondition (and (stack-depth d1) (not (in-stack ncarried)) (bitvalue d1 b2) (not (bitvalue d1 b1)) (not (bitvalue d1 b0)))
        :effect (and (not (stack-depth d1)) (stack-depth d2) (in-stack ncarried) (stack-idx ncarried d2) (not (bitvalue d1 b2)) (bitvalue d1 b1) (bitvalue d1 b0) (bitvalue d2 b0) (bitvalue d2 b1) (bitvalue d2 b2) (bitvalue d2 b3))
    )
    (:action Push_ncarried_d1_b3
        :precondition (and (stack-depth d1) (not (in-stack ncarried)) (bitvalue d1 b3) (not (bitvalue d1 b2)) (not (bitvalue d1 b1)) (not (bitvalue d1 b0)))
        :effect (and (not (stack-depth d1)) (stack-depth d2) (in-stack ncarried) (stack-idx ncarried d2) (not (bitvalue d1 b3)) (bitvalue d1 b2) (bitvalue d1 b1) (bitvalue d1 b0) (bitvalue d2 b0) (bitvalue d2 b1) (bitvalue d2 b2) (bitvalue d2 b3))
    )
    (:action Push_nballs-A_d1_b0
        :precondition (and (stack-depth d1) (not (in-stack nballs-A)) (bitvalue d1 b0))
        :effect (and (not (stack-depth d1)) (stack-depth d2) (in-stack nballs-A) (stack-idx nballs-A d2) (not (bitvalue d1 b0)) (bitvalue d2 b0) (bitvalue d2 b1) (bitvalue d2 b2) (bitvalue d2 b3))
    )
    (:action Push_nballs-A_d1_b1
        :precondition (and (stack-depth d1) (not (in-stack nballs-A)) (bitvalue d1 b1) (not (bitvalue d1 b0)))
        :effect (and (not (stack-depth d1)) (stack-depth d2) (in-stack nballs-A) (stack-idx nballs-A d2) (not (bitvalue d1 b1)) (bitvalue d1 b0) (bitvalue d2 b0) (bitvalue d2 b1) (bitvalue d2 b2) (bitvalue d2 b3))
    )
    (:action Push_nballs-A_d1_b2
        :precondition (and (stack-depth d1) (not (in-stack nballs-A)) (bitvalue d1 b2) (not (bitvalue d1 b1)) (not (bitvalue d1 b0)))
        :effect (and (not (stack-depth d1)) (stack-depth d2) (in-stack nballs-A) (stack-idx nballs-A d2) (not (bitvalue d1 b2)) (bitvalue d1 b1) (bitvalue d1 b0) (bitvalue d2 b0) (bitvalue d2 b1) (bitvalue d2 b2) (bitvalue d2 b3))
    )
    (:action Push_nballs-A_d1_b3
        :precondition (and (stack-depth d1) (not (in-stack nballs-A)) (bitvalue d1 b3) (not (bitvalue d1 b2)) (not (bitvalue d1 b1)) (not (bitvalue d1 b0)))
        :effect (and (not (stack-depth d1)) (stack-depth d2) (in-stack nballs-A) (stack-idx nballs-A d2) (not (bitvalue d1 b3)) (bitvalue d1 b2) (bitvalue d1 b1) (bitvalue d1 b0) (bitvalue d2 b0) (bitvalue d2 b1) (bitvalue d2 b2) (bitvalue d2 b3))
    )
    (:action Push_nfree-grippers_d2_b0
        :precondition (and (stack-depth d2) (bitvalue d2 b0))
        :effect (and (not (stack-depth d2)) (stack-depth d3) (in-stack nfree-grippers) (stack-idx nfree-grippers d3) (not (bitvalue d2 b0)) (bitvalue d3 b0) (bitvalue d3 b1) (bitvalue d3 b2) (bitvalue d3 b3))
    )
    (:action Push_nfree-grippers_d2_b1
        :precondition (and (stack-depth d2) (bitvalue d2 b1) (not (bitvalue d2 b0)))
        :effect (and (not (stack-depth d2)) (stack-depth d3) (in-stack nfree-grippers) (stack-idx nfree-grippers d3) (not (bitvalue d2 b1)) (bitvalue d2 b0) (bitvalue d3 b0) (bitvalue d3 b1) (bitvalue d3 b2) (bitvalue d3 b3))
    )
    (:action Push_nfree-grippers_d2_b2
        :precondition (and (stack-depth d2) (bitvalue d2 b2) (not (bitvalue d2 b1)) (not (bitvalue d2 b0)))
        :effect (and (not (stack-depth d2)) (stack-depth d3) (in-stack nfree-grippers) (stack-idx nfree-grippers d3) (not (bitvalue d2 b2)) (bitvalue d2 b1) (bitvalue d2 b0) (bitvalue d3 b0) (bitvalue d3 b1) (bitvalue d3 b2) (bitvalue d3 b3))
    )
    (:action Push_nfree-grippers_d2_b3
        :precondition (and (stack-depth d2) (bitvalue d2 b3) (not (bitvalue d2 b2)) (not (bitvalue d2 b1)) (not (bitvalue d2 b0)))
        :effect (and (not (stack-depth d2)) (stack-depth d3) (in-stack nfree-grippers) (stack-idx nfree-grippers d3) (not (bitvalue d2 b3)) (bitvalue d2 b2) (bitvalue d2 b1) (bitvalue d2 b0) (bitvalue d3 b0) (bitvalue d3 b1) (bitvalue d3 b2) (bitvalue d3 b3))
    )
    (:action Push_ncarried_d2_b0
        :precondition (and (stack-depth d2) (not (in-stack ncarried)) (bitvalue d2 b0))
        :effect (and (not (stack-depth d2)) (stack-depth d3) (in-stack ncarried) (stack-idx ncarried d3) (not (bitvalue d2 b0)) (bitvalue d3 b0) (bitvalue d3 b1) (bitvalue d3 b2) (bitvalue d3 b3))
    )
    (:action Push_ncarried_d2_b1
        :precondition (and (stack-depth d2) (not (in-stack ncarried)) (bitvalue d2 b1) (not (bitvalue d2 b0)))
        :effect (and (not (stack-depth d2)) (stack-depth d3) (in-stack ncarried) (stack-idx ncarried d3) (not (bitvalue d2 b1)) (bitvalue d2 b0) (bitvalue d3 b0) (bitvalue d3 b1) (bitvalue d3 b2) (bitvalue d3 b3))
    )
    (:action Push_ncarried_d2_b2
        :precondition (and (stack-depth d2) (not (in-stack ncarried)) (bitvalue d2 b2) (not (bitvalue d2 b1)) (not (bitvalue d2 b0)))
        :effect (and (not (stack-depth d2)) (stack-depth d3) (in-stack ncarried) (stack-idx ncarried d3) (not (bitvalue d2 b2)) (bitvalue d2 b1) (bitvalue d2 b0) (bitvalue d3 b0) (bitvalue d3 b1) (bitvalue d3 b2) (bitvalue d3 b3))
    )
    (:action Push_ncarried_d2_b3
        :precondition (and (stack-depth d2) (not (in-stack ncarried)) (bitvalue d2 b3) (not (bitvalue d2 b2)) (not (bitvalue d2 b1)) (not (bitvalue d2 b0)))
        :effect (and (not (stack-depth d2)) (stack-depth d3) (in-stack ncarried) (stack-idx ncarried d3) (not (bitvalue d2 b3)) (bitvalue d2 b2) (bitvalue d2 b1) (bitvalue d2 b0) (bitvalue d3 b0) (bitvalue d3 b1) (bitvalue d3 b2) (bitvalue d3 b3))
    )
    (:action Push_nballs-A_d2_b0
        :precondition (and (stack-depth d2) (not (in-stack nballs-A)) (bitvalue d2 b0))
        :effect (and (not (stack-depth d2)) (stack-depth d3) (in-stack nballs-A) (stack-idx nballs-A d3) (not (bitvalue d2 b0)) (bitvalue d3 b0) (bitvalue d3 b1) (bitvalue d3 b2) (bitvalue d3 b3))
    )
    (:action Push_nballs-A_d2_b1
        :precondition (and (stack-depth d2) (not (in-stack nballs-A)) (bitvalue d2 b1) (not (bitvalue d2 b0)))
        :effect (and (not (stack-depth d2)) (stack-depth d3) (in-stack nballs-A) (stack-idx nballs-A d3) (not (bitvalue d2 b1)) (bitvalue d2 b0) (bitvalue d3 b0) (bitvalue d3 b1) (bitvalue d3 b2) (bitvalue d3 b3))
    )
    (:action Push_nballs-A_d2_b2
        :precondition (and (stack-depth d2) (not (in-stack nballs-A)) (bitvalue d2 b2) (not (bitvalue d2 b1)) (not (bitvalue d2 b0)))
        :effect (and (not (stack-depth d2)) (stack-depth d3) (in-stack nballs-A) (stack-idx nballs-A d3) (not (bitvalue d2 b2)) (bitvalue d2 b1) (bitvalue d2 b0) (bitvalue d3 b0) (bitvalue d3 b1) (bitvalue d3 b2) (bitvalue d3 b3))
    )
    (:action Push_nballs-A_d2_b3
        :precondition (and (stack-depth d2) (not (in-stack nballs-A)) (bitvalue d2 b3) (not (bitvalue d2 b2)) (not (bitvalue d2 b1)) (not (bitvalue d2 b0)))
        :effect (and (not (stack-depth d2)) (stack-depth d3) (in-stack nballs-A) (stack-idx nballs-A d3) (not (bitvalue d2 b3)) (bitvalue d2 b2) (bitvalue d2 b1) (bitvalue d2 b0) (bitvalue d3 b0) (bitvalue d3 b1) (bitvalue d3 b2) (bitvalue d3 b3))
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
    (:action Pop_nfree-grippers_d3
        :precondition (and (stack-depth d3) (in-stack nfree-grippers) (stack-idx nfree-grippers d3))
        :effect (and (stack-depth d2) (not (stack-depth d3)) (not (in-stack nfree-grippers)) (not (stack-idx nfree-grippers d3)))
    )
    (:action Pop_ncarried_d3
        :precondition (and (stack-depth d3) (in-stack ncarried) (stack-idx ncarried d3))
        :effect (and (stack-depth d2) (not (stack-depth d3)) (not (in-stack ncarried)) (not (stack-idx ncarried d3)))
    )
    (:action Pop_nballs-A_d3
        :precondition (and (stack-depth d3) (in-stack nballs-A) (stack-idx nballs-A d3))
        :effect (and (stack-depth d2) (not (stack-depth d3)) (not (in-stack nballs-A)) (not (stack-idx nballs-A d3)))
    )
    (:action action_1
        :precondition (and (not (zero ncarried)) (zero nfree-grippers) (not (robot-at-B)))
        :effect (and (robot-at-B))
    )
    (:action action_2
        :precondition (and (zero ncarried) (not (zero nfree-grippers)) (robot-at-B))
        :effect (and (not (robot-at-B)))
    )
    (:action action_3_d1
        :precondition (and (not (zero ncarried)) (robot-at-B) (in-stack ncarried) (not (in-stack nfree-grippers)) (stack-idx ncarried d1))
        :effect (and (oneof (zero ncarried) (not (zero ncarried))) (not (zero nfree-grippers)) (bitvalue d1 b0) (bitvalue d1 b1) (bitvalue d1 b2) (bitvalue d1 b3) (bitvalue d2 b0) (bitvalue d2 b1) (bitvalue d2 b2) (bitvalue d2 b3) (bitvalue d3 b0) (bitvalue d3 b1) (bitvalue d3 b2) (bitvalue d3 b3))
    )
    (:action action_3_d2
        :precondition (and (not (zero ncarried)) (robot-at-B) (in-stack ncarried) (not (in-stack nfree-grippers)) (stack-idx ncarried d2))
        :effect (and (oneof (zero ncarried) (not (zero ncarried))) (not (zero nfree-grippers)) (bitvalue d2 b0) (bitvalue d2 b1) (bitvalue d2 b2) (bitvalue d2 b3) (bitvalue d3 b0) (bitvalue d3 b1) (bitvalue d3 b2) (bitvalue d3 b3))
    )


    (:action action_4_0_d1
        :precondition (and (not (zero nballs-A)) (not (zero nfree-grippers)) (not (robot-at-B)) (in-stack nballs-A) (not (in-stack ncarried)) (stack-idx nballs-A d1))
        :effect (and (oneof (zero nballs-A) (not (zero nballs-A))) (not (zero ncarried)) (oneof (zero nfree-grippers) (not (zero nfree-grippers))) (bitvalue d1 b0) (bitvalue d1 b1) (bitvalue d1 b2) (bitvalue d1 b3) (bitvalue d2 b0) (bitvalue d2 b1) (bitvalue d2 b2) (bitvalue d2 b3) (bitvalue d3 b0) (bitvalue d3 b1) (bitvalue d3 b2) (bitvalue d3 b3))
    )
    (:action action_4_0_d2
        :precondition (and (not (zero nballs-A)) (not (zero nfree-grippers)) (not (robot-at-B)) (in-stack nballs-A) (not (in-stack ncarried)) (stack-idx nballs-A d2))
        :effect (and (oneof (zero nballs-A) (not (zero nballs-A))) (not (zero ncarried)) (oneof (zero nfree-grippers) (not (zero nfree-grippers))) (bitvalue d2 b0) (bitvalue d2 b1) (bitvalue d2 b2) (bitvalue d2 b3) (bitvalue d3 b0) (bitvalue d3 b1) (bitvalue d3 b2) (bitvalue d3 b3))
    )
    (:action action_4_1_d1
        :precondition (and (not (zero nballs-A)) (not (zero nfree-grippers)) (not (robot-at-B)) (in-stack nfree-grippers) (not (in-stack ncarried)) (stack-idx nfree-grippers d1))
        :effect (and (oneof (zero nballs-A) (not (zero nballs-A))) (not (zero ncarried)) (oneof (zero nfree-grippers) (not (zero nfree-grippers))) (bitvalue d1 b0) (bitvalue d1 b1) (bitvalue d1 b2) (bitvalue d1 b3) (bitvalue d2 b0) (bitvalue d2 b1) (bitvalue d2 b2) (bitvalue d2 b3) (bitvalue d3 b0) (bitvalue d3 b1) (bitvalue d3 b2) (bitvalue d3 b3))
    )
    (:action action_4_1_d2
        :precondition (and (not (zero nballs-A)) (not (zero nfree-grippers)) (not (robot-at-B)) (in-stack nfree-grippers) (not (in-stack ncarried)) (stack-idx nfree-grippers d2))
        :effect (and (oneof (zero nballs-A) (not (zero nballs-A))) (not (zero ncarried)) (oneof (zero nfree-grippers) (not (zero nfree-grippers))) (bitvalue d2 b0) (bitvalue d2 b1) (bitvalue d2 b2) (bitvalue d2 b3) (bitvalue d3 b0) (bitvalue d3 b1) (bitvalue d3 b2) (bitvalue d3 b3))
    )


;   (:action action_4_X_d1
;       :precondition (and (not (zero nballs-A)) (not (zero nfree-grippers)) (not (robot-at-B)) (in-stack nballs-A) (not (in-stack ncarried)) (stack-idx nballs-A d1))
;       :effect (and (oneof (zero nballs-A) (not (zero nballs-A))) (not (zero ncarried)) (zero nfree-grippers) (bitvalue d1 b0) (bitvalue d1 b1) (bitvalue d1 b2) (bitvalue d1 b3) (bitvalue d2 b0) (bitvalue d2 b1) (bitvalue d2 b2) (bitvalue d2 b3) (bitvalue d3 b0) (bitvalue d3 b1) (bitvalue d3 b2) (bitvalue d3 b3))
;   )
;   (:action action_4_X_d2
;       :precondition (and (not (zero nballs-A)) (not (zero nfree-grippers)) (not (robot-at-B)) (in-stack nballs-A) (not (in-stack ncarried)) (stack-idx nballs-A d2))
;       :effect (and (oneof (zero nballs-A) (not (zero nballs-A))) (not (zero ncarried)) (zero nfree-grippers) (bitvalue d2 b0) (bitvalue d2 b1) (bitvalue d2 b2) (bitvalue d2 b3) (bitvalue d3 b0) (bitvalue d3 b1) (bitvalue d3 b2) (bitvalue d3 b3))
;   )


)

