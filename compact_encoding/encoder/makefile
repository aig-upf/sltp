CC=clang++
CCFLAGS=-I. -Wall -std=c++11 -O3

all:			build_alt_theory

build_alt_theory:	build_alt_theory.cc sample.h
			$(CC) $(CCFLAGS) -o build_alt_theory build_alt_theory.cc
clean:			
			rm -f build_alt_theory *~
