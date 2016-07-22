#ifndef CONDENSE_H_
#define CONDENSE_H_

#include "my_structs.h"

#define INCREASING 0
#define DECREASING 1

struct condenseKey {
	int gene;		/* the lower extreme (start or end) */
	int start;		/* start of interval */
	int end;		/* end of interval */
	int diff;		/* difference of extremes of interval */
	int orientation; /* orientation of the sequence: increasing or decreasing */
};
typedef struct condenseKey CondenseKey;
typedef CondenseKey *CondenseKeyPtr;

struct setKeys {
	CondenseKeyPtr *condKeysPtrArray;
	int numKeys;
};
typedef struct setKeys SetKeys;
typedef SetKeys *SetKeysPtr;

void allocateMemoryForKeys( TreePtr phyloTreePtr, SetKeysPtr setPtr );
void freeKeys( TreePtr phyloTreePtr, SetKeysPtr setPtr ); 
void condenseLeafNodes( TreePtr phyloTreePtr, SetKeysPtr setPtr );
void reverseCondensedNodes( TreePtr phyloTreePtr, SetKeysPtr sKeysPtr, int nodeType );
void findSetKeysForCondensingLeaves( TreePtr phyloTreePtr, SetKeysPtr sKeysPtr );

#endif