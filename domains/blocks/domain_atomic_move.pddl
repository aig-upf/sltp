;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; Atomic-move blocksworld
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
(define (domain strips-blocksworld)
  (:requirements :strips)

  (:types
    block - object)

  (:constants table - object)

  (:predicates (on ?x - block ?y - object) (clear ?x - object) (diff ?x - object ?y - object))

  (:action move
    :parameters (?b  - block ?x - object ?y - block)
    :precondition (and (clear ?b) (on ?b ?x) (clear ?y) (diff ?b ?y))
    :effect (and (not (on ?b ?x)) (clear ?x) (not (clear ?y)) (on ?b ?y)))

  (:action move-to-table
    :parameters (?b - block ?x - block)
    :precondition (and (on ?b ?x) (clear ?b))
    :effect (and (not (on ?b ?x)) (clear ?x) (on ?b table)))
)

