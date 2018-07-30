CC = gcc
CXX = g++
OPT = -Wall -O3 -Wno-unused-function -Wno-unused-variable
PSVNOPT = --state_map --no_abstraction --backwards_moves --history_len=0

psvn2c_core.c:
	cp ../psvn2c_core.c ./psvn2c_core.c

psvn2c_state_map.c:
	cp ../psvn2c_state_map.c ./psvn2c_state_map.c

psvn2c_abstraction.c:
	cp ../psvn2c_abstraction.c ./psvn2c_abstraction.c

%.c: %.psvn  psvn2c_core.c psvn2c_state_map.c psvn2c_abstraction.c
	../psvn2c  $(PSVNOPT) --name=$(*F)  < $< > $@
	rm -f ./psvn2c_core.c ./psvn2c_state_map.c ./psvn2c_abstraction.c

.PRECIOUS: %.c

%.succ: %.c succ.c
	$(CC) $(OPT) succ.c -include $< -o $@

priority_queue.hpp:
	cp ../priority_queue.hpp ./priority_queue.hpp

%.dist: %.c dist.cpp priority_queue.hpp
	$(CXX) $(OPT) dist.cpp -include $< -o $@
	rm -f ./priority_queue.hpp

%.features: %.c features.cpp priority_queue.hpp
	$(CXX) $(OPT) features.cpp -include $< -o $@
	rm -f ./priority_queue.hpp

%.distSummary: %.c distSummary.cpp priority_queue.hpp
	$(CXX) $(OPT) distSummary.cpp -include $< -o $@
	rm -f ./priority_queue.hpp

.PHONY: clean
clean:
	rm -f *.succ *.dist *.distSummary *.o  priority_queue.hpp hanoi4p03d.c  hanoi4p03d_bad.c psvn2c_core.c psvn2c_state_map.c psvn2c_abstraction.c *~
