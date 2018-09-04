(define (domain FOND_grid-circles_2_2)
    (:requirements :non-deterministic)
    (:types counter depth bit)
    (:constants num-rewards min-dist[at_Restrict_adjacent_Not_blocked_reward] - counter d0 d1 d2 - depth b0 b1 b2 - bit)
    (:predicates
        (zero ?c - counter)
        (stack-depth ?d - depth)
        (stack-idx ?c - counter ?d - depth)
        (in-stack ?c - counter)
        (bitvalue ?d - depth ?b - bit)
    )

    (:action Push_num-rewards_d0_b0
        :precondition (and (stack-depth d0) (bitvalue d0 b0))
        :effect (and (not (stack-depth d0)) (stack-depth d1) (in-stack num-rewards) (stack-idx num-rewards d1) (not (bitvalue d0 b0)) (bitvalue d1 b0) (bitvalue d1 b1) (bitvalue d1 b2))
    )
    (:action Push_num-rewards_d0_b1
        :precondition (and (stack-depth d0) (bitvalue d0 b1) (not (bitvalue d0 b0)))
        :effect (and (not (stack-depth d0)) (stack-depth d1) (in-stack num-rewards) (stack-idx num-rewards d1) (not (bitvalue d0 b1)) (bitvalue d0 b0) (bitvalue d1 b0) (bitvalue d1 b1) (bitvalue d1 b2))
    )
    (:action Push_num-rewards_d0_b2
        :precondition (and (stack-depth d0) (bitvalue d0 b2) (not (bitvalue d0 b1)) (not (bitvalue d0 b0)))
        :effect (and (not (stack-depth d0)) (stack-depth d1) (in-stack num-rewards) (stack-idx num-rewards d1) (not (bitvalue d0 b2)) (bitvalue d0 b1) (bitvalue d0 b0) (bitvalue d1 b0) (bitvalue d1 b1) (bitvalue d1 b2))
    )
    (:action Push_min-dist[at_Restrict_adjacent_Not_blocked_reward]_d0_b0
        :precondition (and (stack-depth d0) (not (in-stack min-dist[at_Restrict_adjacent_Not_blocked_reward])) (bitvalue d0 b0))
        :effect (and (not (stack-depth d0)) (stack-depth d1) (in-stack min-dist[at_Restrict_adjacent_Not_blocked_reward]) (stack-idx min-dist[at_Restrict_adjacent_Not_blocked_reward] d1) (not (bitvalue d0 b0)) (bitvalue d1 b0) (bitvalue d1 b1) (bitvalue d1 b2))
    )
    (:action Push_min-dist[at_Restrict_adjacent_Not_blocked_reward]_d0_b1
        :precondition (and (stack-depth d0) (not (in-stack min-dist[at_Restrict_adjacent_Not_blocked_reward])) (bitvalue d0 b1) (not (bitvalue d0 b0)))
        :effect (and (not (stack-depth d0)) (stack-depth d1) (in-stack min-dist[at_Restrict_adjacent_Not_blocked_reward]) (stack-idx min-dist[at_Restrict_adjacent_Not_blocked_reward] d1) (not (bitvalue d0 b1)) (bitvalue d0 b0) (bitvalue d1 b0) (bitvalue d1 b1) (bitvalue d1 b2))
    )
    (:action Push_min-dist[at_Restrict_adjacent_Not_blocked_reward]_d0_b2
        :precondition (and (stack-depth d0) (not (in-stack min-dist[at_Restrict_adjacent_Not_blocked_reward])) (bitvalue d0 b2) (not (bitvalue d0 b1)) (not (bitvalue d0 b0)))
        :effect (and (not (stack-depth d0)) (stack-depth d1) (in-stack min-dist[at_Restrict_adjacent_Not_blocked_reward]) (stack-idx min-dist[at_Restrict_adjacent_Not_blocked_reward] d1) (not (bitvalue d0 b2)) (bitvalue d0 b1) (bitvalue d0 b0) (bitvalue d1 b0) (bitvalue d1 b1) (bitvalue d1 b2))
    )
    (:action Push_num-rewards_d1_b0
        :precondition (and (stack-depth d1) (bitvalue d1 b0))
        :effect (and (not (stack-depth d1)) (stack-depth d2) (in-stack num-rewards) (stack-idx num-rewards d2) (not (bitvalue d1 b0)) (bitvalue d2 b0) (bitvalue d2 b1) (bitvalue d2 b2))
    )
    (:action Push_num-rewards_d1_b1
        :precondition (and (stack-depth d1) (bitvalue d1 b1) (not (bitvalue d1 b0)))
        :effect (and (not (stack-depth d1)) (stack-depth d2) (in-stack num-rewards) (stack-idx num-rewards d2) (not (bitvalue d1 b1)) (bitvalue d1 b0) (bitvalue d2 b0) (bitvalue d2 b1) (bitvalue d2 b2))
    )
    (:action Push_num-rewards_d1_b2
        :precondition (and (stack-depth d1) (bitvalue d1 b2) (not (bitvalue d1 b1)) (not (bitvalue d1 b0)))
        :effect (and (not (stack-depth d1)) (stack-depth d2) (in-stack num-rewards) (stack-idx num-rewards d2) (not (bitvalue d1 b2)) (bitvalue d1 b1) (bitvalue d1 b0) (bitvalue d2 b0) (bitvalue d2 b1) (bitvalue d2 b2))
    )
    (:action Push_min-dist[at_Restrict_adjacent_Not_blocked_reward]_d1_b0
        :precondition (and (stack-depth d1) (not (in-stack min-dist[at_Restrict_adjacent_Not_blocked_reward])) (bitvalue d1 b0))
        :effect (and (not (stack-depth d1)) (stack-depth d2) (in-stack min-dist[at_Restrict_adjacent_Not_blocked_reward]) (stack-idx min-dist[at_Restrict_adjacent_Not_blocked_reward] d2) (not (bitvalue d1 b0)) (bitvalue d2 b0) (bitvalue d2 b1) (bitvalue d2 b2))
    )
    (:action Push_min-dist[at_Restrict_adjacent_Not_blocked_reward]_d1_b1
        :precondition (and (stack-depth d1) (not (in-stack min-dist[at_Restrict_adjacent_Not_blocked_reward])) (bitvalue d1 b1) (not (bitvalue d1 b0)))
        :effect (and (not (stack-depth d1)) (stack-depth d2) (in-stack min-dist[at_Restrict_adjacent_Not_blocked_reward]) (stack-idx min-dist[at_Restrict_adjacent_Not_blocked_reward] d2) (not (bitvalue d1 b1)) (bitvalue d1 b0) (bitvalue d2 b0) (bitvalue d2 b1) (bitvalue d2 b2))
    )
    (:action Push_min-dist[at_Restrict_adjacent_Not_blocked_reward]_d1_b2
        :precondition (and (stack-depth d1) (not (in-stack min-dist[at_Restrict_adjacent_Not_blocked_reward])) (bitvalue d1 b2) (not (bitvalue d1 b1)) (not (bitvalue d1 b0)))
        :effect (and (not (stack-depth d1)) (stack-depth d2) (in-stack min-dist[at_Restrict_adjacent_Not_blocked_reward]) (stack-idx min-dist[at_Restrict_adjacent_Not_blocked_reward] d2) (not (bitvalue d1 b2)) (bitvalue d1 b1) (bitvalue d1 b0) (bitvalue d2 b0) (bitvalue d2 b1) (bitvalue d2 b2))
    )
    (:action Pop_num-rewards_d1
        :precondition (and (stack-depth d1) (in-stack num-rewards) (stack-idx num-rewards d1))
        :effect (and (stack-depth d0) (not (stack-depth d1)) (not (in-stack num-rewards)) (not (stack-idx num-rewards d1)))
    )
    (:action Pop_min-dist[at_Restrict_adjacent_Not_blocked_reward]_d1
        :precondition (and (stack-depth d1) (in-stack min-dist[at_Restrict_adjacent_Not_blocked_reward]) (stack-idx min-dist[at_Restrict_adjacent_Not_blocked_reward] d1))
        :effect (and (stack-depth d0) (not (stack-depth d1)) (not (in-stack min-dist[at_Restrict_adjacent_Not_blocked_reward])) (not (stack-idx min-dist[at_Restrict_adjacent_Not_blocked_reward] d1)))
    )
    (:action Pop_num-rewards_d2
        :precondition (and (stack-depth d2) (in-stack num-rewards) (stack-idx num-rewards d2))
        :effect (and (stack-depth d1) (not (stack-depth d2)) (not (in-stack num-rewards)) (not (stack-idx num-rewards d2)))
    )
    (:action Pop_min-dist[at_Restrict_adjacent_Not_blocked_reward]_d2
        :precondition (and (stack-depth d2) (in-stack min-dist[at_Restrict_adjacent_Not_blocked_reward]) (stack-idx min-dist[at_Restrict_adjacent_Not_blocked_reward] d2))
        :effect (and (stack-depth d1) (not (stack-depth d2)) (not (in-stack min-dist[at_Restrict_adjacent_Not_blocked_reward])) (not (stack-idx min-dist[at_Restrict_adjacent_Not_blocked_reward] d2)))
    )
    (:action action_1_d1
        :precondition (and (zero min-dist[at_Restrict_adjacent_Not_blocked_reward]) (not (zero num-rewards)) (in-stack num-rewards) (not (in-stack min-dist[at_Restrict_adjacent_Not_blocked_reward])) (stack-idx num-rewards d1))
        :effect (and (not (zero min-dist[at_Restrict_adjacent_Not_blocked_reward])) (oneof (zero num-rewards) (not (zero num-rewards))) (bitvalue d1 b0) (bitvalue d1 b1) (bitvalue d1 b2) (bitvalue d2 b0) (bitvalue d2 b1) (bitvalue d2 b2))
    )
    (:action action_1_d2
        :precondition (and (zero min-dist[at_Restrict_adjacent_Not_blocked_reward]) (not (zero num-rewards)) (in-stack num-rewards) (not (in-stack min-dist[at_Restrict_adjacent_Not_blocked_reward])) (stack-idx num-rewards d2))
        :effect (and (not (zero min-dist[at_Restrict_adjacent_Not_blocked_reward])) (oneof (zero num-rewards) (not (zero num-rewards))) (bitvalue d2 b0) (bitvalue d2 b1) (bitvalue d2 b2))
    )
    (:action action_2_d1
        :precondition (and (not (zero min-dist[at_Restrict_adjacent_Not_blocked_reward])) (not (zero num-rewards)) (in-stack min-dist[at_Restrict_adjacent_Not_blocked_reward]) (stack-idx min-dist[at_Restrict_adjacent_Not_blocked_reward] d1))
        :effect (and (oneof (zero min-dist[at_Restrict_adjacent_Not_blocked_reward]) (not (zero min-dist[at_Restrict_adjacent_Not_blocked_reward]))) (bitvalue d1 b0) (bitvalue d1 b1) (bitvalue d1 b2) (bitvalue d2 b0) (bitvalue d2 b1) (bitvalue d2 b2))
    )
    (:action action_2_d2
        :precondition (and (not (zero min-dist[at_Restrict_adjacent_Not_blocked_reward])) (not (zero num-rewards)) (in-stack min-dist[at_Restrict_adjacent_Not_blocked_reward]) (stack-idx min-dist[at_Restrict_adjacent_Not_blocked_reward] d2))
        :effect (and (oneof (zero min-dist[at_Restrict_adjacent_Not_blocked_reward]) (not (zero min-dist[at_Restrict_adjacent_Not_blocked_reward]))) (bitvalue d2 b0) (bitvalue d2 b1) (bitvalue d2 b2))
    )
)

