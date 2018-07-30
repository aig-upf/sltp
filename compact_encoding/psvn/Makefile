CC = gcc
CXX = g++
#OPT = -g
OPT = -g -Wall -O3


psvn2c: psvn2c.cpp psvn2c.hpp psvn.cpp psvn.hpp
	$(CXX) $(OPT) psvn2c.cpp psvn.cpp -o $@

clean:
	rm -f psvn2c *.o *~
