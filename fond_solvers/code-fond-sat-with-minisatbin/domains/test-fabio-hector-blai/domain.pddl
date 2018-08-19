(define (domain list-search)
  (:requirements :non-deterministic)
  (:predicates
    (goalNode); goal Node has been reached
    (done); fluent to express that goal has been achieved (used for disjunctive goal)
    (hasNext); current node has a successor

    ;auxiliary fluents: enable regular or auxiliary actions
    (actionEnabled); regular action enabled
    (checkGoalEnabled); auxiliary action enabled (to check disjunctive goals)
  )

  (:action next
    :precondition (and (actionEnabled) (hasNext) (not (goalNode)))
    :effect
      (oneof
        (and (not(actionEnabled)) (checkGoalEnabled) (hasNext) (goalNode))
        (and (not(actionEnabled)) (checkGoalEnabled) (not(hasNext)) (goalNode))
        (and (not(actionEnabled)) (checkGoalEnabled) (not(hasNext)) (not(goalNode)))
        (and (not(actionEnabled)) (checkGoalEnabled) (hasNext) (not(goalNode)))
      )
  )

  (:action checkGoal1
    :precondition (and (checkGoalEnabled) (goalNode))
    :effect (and (not (checkGoalEnabled)) (actionEnabled) (done))
  )

  (:action checkGoal2
    :precondition (and (checkGoalEnabled) (not (hasNext)))
    :effect (and (not (checkGoalEnabled)) (actionEnabled) (done))
  )

  (:action checkGoal3
    :precondition (and (checkGoalEnabled) (not (goalNode)) (hasNext))
    :effect (and (not (checkGoalEnabled)) (actionEnabled) (not (done)))
  )

)
