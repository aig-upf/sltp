/* This program reads a state from stdin and prints out its successors.

Copyright (C) 2013 by the PSVN Research Group, University of Alberta
*/

#include <stdlib.h>
#include <stdio.h>
#include <inttypes.h>
#include <assert.h>
#include <sys/time.h>

#define  MAX_LINE_LENGTH 999 

int main( int argc, char **argv )
{
// VARIABLES FOR INPUT
    char str[ MAX_LINE_LENGTH +1 ] ;
    ssize_t nchars; 
    state_t state; // state_t is defined by the PSVN API. It is the type used for individual states.

// VARIABLES FOR ITERATING THROUGH state's SUCCESSORS
    state_t child;
    ruleid_iterator_t iter; // ruleid_terator_t is the type defined by the PSVN API successor/predecessor iterators.
    int ruleid ; // an iterator returns a number identifying a rule
    int childCount;

// READ A LINE OF INPUT FROM stdin
    printf("Please enter a state followed by ENTER: ");
    if ( fgets(str, sizeof str, stdin) == NULL ) {
	printf("Error: empty input line.\n");
	return 0; 
    }

// CONVERT THE STRING TO A STATE
    nchars = read_state( str, &state );
    if (nchars <= 0) {
	printf("Error: invalid state entered.\n");
	return 0; 
    }

    printf("The state you entered is: ");
    print_state( stdout, &state );
    printf("\n");

// LOOP THOUGH THE CHILDREN ONE BY ONE
    childCount = 0;
    init_fwd_iter( &iter, &state );  // initialize the child iterator 
      while( ( ruleid = next_ruleid( &iter ) ) >= 0 ) {
//    while( -1 >= 0 ) {
	apply_fwd_rule( ruleid, &state, &child );
	++childCount;
    	printf("child %d. ",childCount);
	print_state( stdout, &child );
    	//printf("  %s (cost %d)\n",get_fwd_rule_label(ruleid),get_fwd_rule_cost(ruleid));
    	printf("  %s (cost %d), goal=%d\n",get_fwd_rule_label(ruleid),get_fwd_rule_cost(ruleid),is_goal(&child));
    } // end while... no more children
    if (childCount == 0) {
	printf("Your state has no children.\n");
    }

    return 0;
} // end main
