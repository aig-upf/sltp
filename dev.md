

./run.py --instance /home/frances/projects/code/downward-benchmarks/gripper/prob01.pddl  --driver bfs --disable-static-analysis --options="max_expansions=100"




## Sanity Checks

`sort feature-matrix.txt.int | uniq -cd` should return no output, which means that there is no duplicate row in the
feature matrix.