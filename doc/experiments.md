

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
    ./logistics.py sample1_tg --all
    
    # IPC Barman domain
    ./barman.py prob01_tg --all
    
    # IPC Childsnack domain
    ./childsnack.py prob01_tg --all

    # IPC Spanner domain
    ./spanner.py exp1 --all
    
    # IPC Grid domain
    ./grid.py prob01 --all
    
    # IPC Miconic domain
    ./miconic.py p1 --all
    
    # Blocksworld: clear two particular blocks
    ./blocks.py clear_two_atoms --all  # goal: clear(x) and clear(y)
    
    # IPC Blocksworld domain: construct one arbitrary tower
    ./blocks2.py arbitrary1 --all 

    # Taxi domain: Grid-like environment where a taxi needs to move to the cell where a passenger is,
    # pick her up, take her to some destination cell and drop her there
    # ./taxi.py simple --all  # Not working: pipeline needs to be adjusted for FSTRIPS domains

## Variations of previous benchmarks

    # Blocksworld: stack one block on top of another, one single training instance
    ./blocks2.py on_x_y --all 
    
    # Blocksworld: stack one block on top of another, one single training instance, w./ goal-concepts
    ./blocks2.py on_x_y_gc --all
