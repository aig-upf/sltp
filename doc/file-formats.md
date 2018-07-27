

## `transition-matrix.dat`

Each line in the file is of the form

    s t_1 t_2 ... t_n

and corresponds to all transitions (s, t_i) that start in state s.
State IDs are not necessarily contiguous!


## `feature-matrix.dat`

Each line in the file is of the form

    s f^s_1 f^s_2 ... f^s_n

and corresponds to the denotation of all features in state s.
State IDs are not necessarily contiguous!

## `goal-states.dat`
The file contains a single line

    s_1 s_2 ... s_n

with the IDs of all goal states.
