

# Some interesting experiments


## AAAI paper experiments

    # Blocksworld: clear a single block
    ./blocks.py aaai_clear_x_simple_hybrid --all  # goal: clear(x)
    
    # Blocksworld: stack one block on top of another
    ./blocks.py aaai_bw_on_x_y_completeness_opt --all  # goal: on(x, y)
    
    # Gripper: take all balls to a certain room
    ./gripper.py aaai_prob01 --all
    

## Other benchmarks

    # IPC Logistics domain
    ./run.py logistics:p1 --all
    
    # IPC Barman domain
    ./run.py barman:p1 --all
    
    # IPC Childsnack domain
    ./run.py childsnack:p1 --all

    # IPC Spanner domain
    ./run.py spanner:p1 --all
    
    # IPC Grid domain
    ./run.py grid:p1 --all

    # IPC Miconic domain
    ./run.py miconic:p1 --all
    
    # Towers of Hanoi
    ./run.py hanoi:p1 --all
    
    # IPC Visitall domain
    ./run.py visitall:p1 --all
    
    # IPC depot domain
    ./run.py depot:p1 --all
    
    # IPC Satellite domain
    ./run.py satellite:p1 --all
    
    # Blocksworld: clear two particular blocks
    ./blocks.py clear_two_atoms --all  # goal: clear(x) and clear(y)
    
    # IPC Blocksworld domain: construct one arbitrary tower
    ./run.py blocks2:arbitrary1 --all

    # Taxi domain: Grid-like environment where a taxi needs to move to the cell where a passenger is,
    # pick her up, take her to some destination cell and drop her there
    # ./taxi.py simple --all  # Not working: pipeline needs to be adjusted for FSTRIPS domains

## Variations of previous benchmarks

    # Blocksworld: stack one block on top of another, one single training instance
    ./run.py blocks2:on_x_y --all
