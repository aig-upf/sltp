(define (domain FOND_gripper_2_2)
    (:requirements :non-deterministic)
    (:types counter depth bit)
    (:constants nfree-grippers ncarried nballs-A - counter d0 d1 d2 - depth b0 b1 b2 - bit)
    (:predicates
        (zero ?c - counter)
        (stack-depth ?d - depth)
        (stack-idx ?c - counter ?d - depth)
        (in-stack ?c - counter)
        (bitvalue ?d - depth ?b - bit)
        (robot-at-B)
    )

    (:action PUSH_d0_b0
        :parameters (?c - counter)
        :precondition (and (stack-depth d0) (bitvalue d0 b0))
        :effect (and (not (stack-depth d0)) (stack-depth d1) (in-stack ?c) (stack-idx ?c d1) (not (bitvalue d0 b0)) (bitvalue d1 b2) (bitvalue d1 b1) (bitvalue d1 b0))
    )
    (:action PUSH_d0_b1
        :parameters (?c - counter)
        :precondition (and (stack-depth d0) (bitvalue d0 b1) (not (bitvalue d0 b0)))
        :effect (and (not (stack-depth d0)) (stack-depth d1) (in-stack ?c) (stack-idx ?c d1) (not (bitvalue d0 b1)) (bitvalue d0 b0) (bitvalue d1 b2) (bitvalue d1 b1) (bitvalue d1 b0))
    )
    (:action PUSH_d0_b2
        :parameters (?c - counter)
        :precondition (and (stack-depth d0) (bitvalue d0 b2) (not (bitvalue d0 b1)) (not (bitvalue d0 b0)))
        :effect (and (not (stack-depth d0)) (stack-depth d1) (in-stack ?c) (stack-idx ?c d1) (not (bitvalue d0 b2)) (bitvalue d0 b1) (bitvalue d0 b0) (bitvalue d1 b2) (bitvalue d1 b1) (bitvalue d1 b0))
    )
    (:action PUSH_d1_b0
        :parameters (?c - counter)
        :precondition (and (stack-depth d1) (not (in-stack ?c)) (bitvalue d1 b0))
        :effect (and (not (stack-depth d1)) (stack-depth d2) (in-stack ?c) (stack-idx ?c d2) (not (bitvalue d1 b0)) (bitvalue d2 b2) (bitvalue d2 b1) (bitvalue d2 b0))
    )
    (:action PUSH_d1_b1
        :parameters (?c - counter)
        :precondition (and (stack-depth d1) (not (in-stack ?c)) (bitvalue d1 b1) (not (bitvalue d1 b0)))
        :effect (and (not (stack-depth d1)) (stack-depth d2) (in-stack ?c) (stack-idx ?c d2) (not (bitvalue d1 b1)) (bitvalue d1 b0) (bitvalue d2 b2) (bitvalue d2 b1) (bitvalue d2 b0))
    )
    (:action PUSH_d1_b2
        :parameters (?c - counter)
        :precondition (and (stack-depth d1) (not (in-stack ?c)) (bitvalue d1 b2) (not (bitvalue d1 b1)) (not (bitvalue d1 b0)))
        :effect (and (not (stack-depth d1)) (stack-depth d2) (in-stack ?c) (stack-idx ?c d2) (not (bitvalue d1 b2)) (bitvalue d1 b1) (bitvalue d1 b0) (bitvalue d2 b2) (bitvalue d2 b1) (bitvalue d2 b0))
    )
    (:action POP_d1
        :parameters (?c - counter)
        :precondition (and (stack-depth d1) (in-stack ?c) (stack-idx ?c d1))
        :effect (and (not (stack-depth d1)) (stack-depth d0) (not (in-stack ?c)) (not (stack-idx ?c d1)))
    )
    (:action POP_d2
        :parameters (?c - counter)
        :precondition (and (stack-depth d2) (in-stack ?c) (stack-idx ?c d2))
        :effect (and (not (stack-depth d2)) (stack-depth d1) (not (in-stack ?c)) (not (stack-idx ?c d2)))
    )
    (:action action_1_d1
        :precondition (and (not (zero ncarried)) (robot-at-B) (in-stack ncarried) (not (in-stack nfree-grippers)) (stack-idx ncarried d1))
        :effect (and (oneof (zero ncarried) (not (zero ncarried))) (not (zero nfree-grippers)) (bitvalue d1 b0) (bitvalue d1 b1) (bitvalue d1 b2) (bitvalue d2 b0) (bitvalue d2 b1) (bitvalue d2 b2))
    )
    (:action action_1_d2
        :precondition (and (not (zero ncarried)) (robot-at-B) (in-stack ncarried) (not (in-stack nfree-grippers)) (stack-idx ncarried d2))
        :effect (and (oneof (zero ncarried) (not (zero ncarried))) (not (zero nfree-grippers)) (bitvalue d2 b0) (bitvalue d2 b1) (bitvalue d2 b2))
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
        :effect (and (oneof (zero nballs-A) (not (zero nballs-A))) (not (zero ncarried)) (oneof (zero nfree-grippers) (not (zero nfree-grippers))) (bitvalue d1 b0) (bitvalue d1 b1) (bitvalue d1 b2) (bitvalue d2 b0) (bitvalue d2 b1) (bitvalue d2 b2))
    )
    (:action action_4_0_d2
        :precondition (and (not (zero nballs-A)) (not (zero nfree-grippers)) (not (robot-at-B)) (in-stack nballs-A) (not (in-stack ncarried)) (stack-idx nballs-A d2))
        :effect (and (oneof (zero nballs-A) (not (zero nballs-A))) (not (zero ncarried)) (oneof (zero nfree-grippers) (not (zero nfree-grippers))) (bitvalue d2 b0) (bitvalue d2 b1) (bitvalue d2 b2))
    )
    (:action action_4_1_d1
        :precondition (and (not (zero nballs-A)) (not (zero nfree-grippers)) (not (robot-at-B)) (in-stack nfree-grippers) (not (in-stack ncarried)) (stack-idx nfree-grippers d1))
        :effect (and (oneof (zero nballs-A) (not (zero nballs-A))) (not (zero ncarried)) (oneof (zero nfree-grippers) (not (zero nfree-grippers))) (bitvalue d1 b0) (bitvalue d1 b1) (bitvalue d1 b2) (bitvalue d2 b0) (bitvalue d2 b1) (bitvalue d2 b2))
    )
    (:action action_4_1_d2
        :precondition (and (not (zero nballs-A)) (not (zero nfree-grippers)) (not (robot-at-B)) (in-stack nfree-grippers) (not (in-stack ncarried)) (stack-idx nfree-grippers d2))
        :effect (and (oneof (zero nballs-A) (not (zero nballs-A))) (not (zero ncarried)) (oneof (zero nfree-grippers) (not (zero nfree-grippers))) (bitvalue d2 b0) (bitvalue d2 b1) (bitvalue d2 b2))
    )
    (:action action_5
        :precondition (and (zero ncarried) (not (zero nfree-grippers)) (robot-at-B))
        :effect (and (not (robot-at-B)))
    )
)

