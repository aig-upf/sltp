

# Some interesting experiments


## AAAI paper experiments

    # Blocksworld: clear a single block
    ./blocks.py aaai_clear_x_simple_hybrid  # goal: clear(x)
    
    # Blocksworld: stack one block on top of another
    ./blocks.py aaai_bw_on_x_y_completeness_opt  # goal: on(x, y)
    
    # Gripper: take all balls to a certain room
    ./gripper.py aaai_prob01
    

## Other benchmarks

    # IPC Logistics domain
    ./run.py logistics:p1
    
    # IPC Barman domain
    ./run.py barman:p1
    
    # IPC Childsnack domain
    ./run.py childsnack:p1

    # IPC Spanner domain
    ./run.py spanner:p1
    
    # IPC Grid domain
    ./run.py grid:p1

    # IPC Miconic domain
    ./run.py miconic:p1
    
    # Towers of Hanoi
    ./run.py hanoi:p1
    
    # IPC Visitall domain
    ./run.py visitall:p1
    
    # IPC depot domain
    ./run.py depot:p1
    
    # IPC Satellite domain
    ./run.py satellite:p1
    
    # Blocksworld: clear two particular blocks
    ./blocks.py clear_two_atoms  # goal: clear(x) and clear(y)
    
    # IPC Blocksworld domain: construct one arbitrary tower
    ./run.py blocks2:one_tower_inc

    # Taxi domain: Grid-like environment where a taxi needs to move to the cell where a passenger is,
    # pick her up, take her to some destination cell and drop her there
    # ./taxi.py simple  # Not working: pipeline needs to be adjusted for FSTRIPS domains

## Variations of previous benchmarks

    # Blocksworld: stack one block on top of another, one single training instance
    ./run.py blocks2:on_x_y


# Running experiments in a Slurm cluster

Assuming you have the experiment description in a file `expname.yml`, the full pipeline can be run in
a Slurm cluster in two steps: first, generate the experiment script automatically:

    ./cluster.py --exp expname

And then just run that experiment script with the standard Slurm command: 

    sbatch expname.sh
