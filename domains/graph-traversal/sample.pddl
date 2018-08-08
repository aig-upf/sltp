
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; Instance file automatically generated by the Tarski FSTRIPS writer
;;; 
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(define (problem grid-circles-4x4)
    (:domain grid-circles-strips)

    (:objects
        c_0_0 c_0_1 c_0_2 c_0_3 c_1_0 c_1_1 c_1_2 c_1_3 - cell
    )

    (:init
        (reward c_1_1)
        (reward c_0_3)
        (at c_0_1)
        (adjacent c_1_2 c_1_1)
        (adjacent c_0_2 c_0_3)
        (adjacent c_0_3 c_0_2)
        (adjacent c_1_1 c_0_1)
        (adjacent c_1_2 c_1_3)
        (adjacent c_1_3 c_0_3)
        (adjacent c_1_1 c_1_0)
        (adjacent c_0_0 c_0_1)
        (adjacent c_0_1 c_0_0)
        (adjacent c_0_3 c_1_3)
        (adjacent c_0_2 c_1_2)
        (adjacent c_1_2 c_0_2)
        (adjacent c_1_0 c_0_0)
        (adjacent c_0_0 c_1_0)
        (adjacent c_1_1 c_1_2)
        (adjacent c_1_0 c_1_1)
        (adjacent c_0_2 c_0_1)
        (adjacent c_0_1 c_0_2)
        (adjacent c_0_1 c_1_1)
        (adjacent c_1_3 c_1_2)
    )

    (:goal
        (and (not (reward c_0_0)) (not (reward c_0_1)) (not (reward c_0_2)) (not (reward c_0_3)) (not (reward c_1_0)) (not (reward c_1_1)) (not (reward c_1_2)) (not (reward c_1_3)) )
    )

    
    
    
)

