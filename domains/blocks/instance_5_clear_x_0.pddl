(define (problem BLOCKS-5-CLEAR-X)
(:domain BLOCKS)
(:objects A B C D E)
(:init (CLEAR A) (CLEAR B) (CLEAR C) (CLEAR D) (CLEAR E)
(ONTABLE A) (ONTABLE B) (ONTABLE C) (ONTABLE D) (ONTABLE E) (HANDEMPTY))
(:goal (AND (CLEAR A)))
)
