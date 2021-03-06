
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; Instance file automatically generated by the Tarski FSTRIPS writer
;;; 
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(define (problem gridworld-10x10)
    (:domain gridworld-strips)

    (:objects
        c1 c10 c2 c3 c4 c5 c6 c7 c8 c9 - coordinate
    )

    (:init
        (= (xpos ) c1)
        (= (ypos ) c1)
        (= (maxpos ) c10)
        (= (goal_xpos ) c10)
        (= (goal_ypos ) c10)
        (succ c8 c9)
        (succ c7 c8)
        (succ c9 c10)
        (succ c4 c5)
        (succ c6 c7)
        (succ c1 c2)
        (succ c3 c4)
        (succ c2 c3)
        (succ c5 c6)
    )

    (:goal
        (and (= (xpos ) (goal_xpos )) (= (ypos ) (goal_ypos )))
    )

    
    
    
)

