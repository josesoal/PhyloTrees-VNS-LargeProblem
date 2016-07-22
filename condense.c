
#include <stdio.h>
#include <stdlib.h>

#include "condense.h"
#include "auxiliary.h"
#include "my_structs.h"
#include "tree.h"

static int findIntersection(CondenseKeyPtr sequenceAPtr, 
	CondenseKeyPtr sequenceBPtr, int *start, int *end);
static int isIncludedIn( int a, int start, int end); 
static void addUniqueKey(SetKeysPtr sKeysPtr, int *ikey, int start, int end);

void allocateMemoryForKeys(  TreePtr phyloTreePtr, SetKeysPtr setPtr ) 
{
	int i,  maxNumKeys;

	maxNumKeys = phyloTreePtr->numberGenes / 2;
	setPtr->condKeysPtrArray = malloc( maxNumKeys * sizeof( CondenseKeyPtr ) );
	if ( setPtr->condKeysPtrArray == NULL ) { nomemMessage( "setPtr->condKeysPtrArray" ); }

	for ( i = 0; i < maxNumKeys; i++ ) {
		setPtr->condKeysPtrArray[i] = malloc( sizeof( CondenseKey ) );
		if ( setPtr->condKeysPtrArray[i] == NULL ) { nomemMessage( "setPtr->condKeysPtrArray[i]" ); }
	}
}

void freeKeys( TreePtr phyloTreePtr, SetKeysPtr setPtr ) 
{
	int i,  maxNumKeys;

	maxNumKeys = phyloTreePtr->numberGenes / 2;
	for ( i = 0; i < maxNumKeys; i++ ) {
		free( setPtr->condKeysPtrArray[ i ] );
	}
	free( setPtr->condKeysPtrArray );
}

void condenseLeafNodes( TreePtr phyloTreePtr, SetKeysPtr sKeysPtr )
{
	int i, j, k, iwrite, otherExtreme;

	for ( i = 0; i < sKeysPtr->numKeys; i++ ) {
		/* determine the other extreme of the interval */
		otherExtreme = sKeysPtr->condKeysPtrArray[i]->start;
		if ( sKeysPtr->condKeysPtrArray[i]->gene == 
				sKeysPtr->condKeysPtrArray[i]->start) {
			otherExtreme = sKeysPtr->condKeysPtrArray[i]->end;
		} 
		/* condense leaves using key i */
		for ( j = 0; j < phyloTreePtr->numberLeaves; j++ ) {
			
			iwrite = 0;
			for ( k = 0; k < phyloTreePtr->numberGenes; k++ ) {

				if ( abs(phyloTreePtr->nodesPtrArray[ j ]->genome[ k ]) <= 
						sKeysPtr->condKeysPtrArray[ i ]->gene ) {

					phyloTreePtr->nodesPtrArray[ j ]->genome[ iwrite ] = 
						phyloTreePtr->nodesPtrArray[ j ]->genome[ k ];
					iwrite++;
				}
				else if ( abs(phyloTreePtr->nodesPtrArray[ j ]->genome[ k ]) > 
							otherExtreme ) {

					/* update the gen by subtracting (or adding) the difference */
					if ( phyloTreePtr->nodesPtrArray[ j ]->genome[ k ] > 0 ) {

						phyloTreePtr->nodesPtrArray[ j ]->genome[ iwrite ] =
								phyloTreePtr->nodesPtrArray[ j ]->genome[ k ] - 
								sKeysPtr->condKeysPtrArray[ i ]->diff;
					}
					else { // phyloTreePtr->nodesPtrArray[ j ]->genome[ k ] < 0

						phyloTreePtr->nodesPtrArray[ j ]->genome[ iwrite ] =
								phyloTreePtr->nodesPtrArray[ j ]->genome[ k ] + 
								sKeysPtr->condKeysPtrArray[ i ]->diff;
					}
					iwrite++;
				}
			} 
		}
		phyloTreePtr->numberGenes = iwrite;
	}	
	//show keys ...
	/*for ( i = 0; i < sKeysPtr->numKeys; i++ ) {
		printf( "*gene(%d) --> start %d, end %d\n", 
			sKeysPtr->condKeysPtrArray[ i ]->gene, 
			sKeysPtr->condKeysPtrArray[ i ]->start, 
			sKeysPtr->condKeysPtrArray[ i ]->end );
	}*/ 
}

/* reverse the condensation of leaf nodes "LEAF_NODE" ( or all nodes "EACH_NODE") */
void reverseCondensedNodes( TreePtr phyloTreePtr, SetKeysPtr sKeysPtr, int nodeType )
{
	int i, j, k, iseq, diff, gene, numberNodes;

	numberNodes = phyloTreePtr->numberNodes;
	if ( nodeType == LEAF_NODE) {
		numberNodes = phyloTreePtr->numberLeaves;
	}

	for ( i = sKeysPtr->numKeys - 1; i >= 0; i-- ) {
		diff = sKeysPtr->condKeysPtrArray[ i ]->diff;
		/* reverse condensation of nodes using key i */
		for ( j = 0; j < numberNodes; j++ ) {
			/* find gene that represents a sequence */
			for ( k = 0; k < phyloTreePtr->numberGenes; k++ ) {
				if ( abs( phyloTreePtr->nodesPtrArray[ j ]->genome[ k ] ) > 
							sKeysPtr->condKeysPtrArray[ i ]->gene ) {

					if ( phyloTreePtr->nodesPtrArray[ j ]->genome[ k ] > 0 ) {
						phyloTreePtr->nodesPtrArray[ j ]->genome[ k ] = 
							phyloTreePtr->nodesPtrArray[ j ]->genome[ k ] + diff;
					}
					else { // < 0
						phyloTreePtr->nodesPtrArray[ j ]->genome[ k ] = 
							phyloTreePtr->nodesPtrArray[ j ]->genome[ k ] - diff;
					}
				}
				else if ( abs( phyloTreePtr->nodesPtrArray[ j ]->genome[ k ] ) == 
								sKeysPtr->condKeysPtrArray[ i ]->gene ) {
					iseq = k;
					break;
				}
			}

			/* make space by shifting to the right */
			for ( k = phyloTreePtr->numberGenes-1; k > iseq; k-- ) {
				if ( abs( phyloTreePtr->nodesPtrArray[ j ]->genome[ k ] ) > 
							sKeysPtr->condKeysPtrArray[ i ]->gene ) {

					if ( phyloTreePtr->nodesPtrArray[ j ]->genome[ k ] > 0 ) {
						phyloTreePtr->nodesPtrArray[ j ]->genome[ k + diff ] = 
							phyloTreePtr->nodesPtrArray[ j ]->genome[ k ] + diff;
					}
					else{ // < 0
						phyloTreePtr->nodesPtrArray[ j ]->genome[ k + diff ] = 
							phyloTreePtr->nodesPtrArray[ j ]->genome[ k ] - diff;
					}
				}
				else { // <
					phyloTreePtr->nodesPtrArray[ j ]->genome[ k + diff ] = 
							phyloTreePtr->nodesPtrArray[ j ]->genome[ k ];
				}
			}

			/* insert the sequence  */
			if ( phyloTreePtr->nodesPtrArray[ j ]->genome[ iseq ] > 0 ) { 
				gene = sKeysPtr->condKeysPtrArray[ i ]->start; 
			}
			else { // < 0
				gene = -1 * sKeysPtr->condKeysPtrArray[ i ]->end; 	
			}
			
			if ( sKeysPtr->condKeysPtrArray[ i ]->orientation == INCREASING ) {
				for ( k = iseq; k <= iseq + diff; k++ ) {
					phyloTreePtr->nodesPtrArray[ j ]->genome[ k ] = gene;
					gene = gene + 1; 	
				}
			}
			else { // sKeysPtr->condKeysPtrArray[i]->orientation == DECREASING
				for ( k = iseq; k <= iseq + diff; k++ ) {
					phyloTreePtr->nodesPtrArray[ j ]->genome[ k ] = gene;
					gene = gene - 1; 	
				}
			}
		}//end-for
		phyloTreePtr->numberGenes = phyloTreePtr->numberGenes + diff;

	} //end-for
}

void findSetKeysForCondensingLeaves( TreePtr phyloTreePtr, SetKeysPtr sKeysPtr )
{
	int i, j, kk, hh, ikey, maxNumKeys;
	int current, next, start, end, startFound, found, sequenceComplete;

	/* allocate memory for array of setkeys */
	SetKeysPtr *arrPtr;
	maxNumKeys = phyloTreePtr->numberGenes / 2;
	arrPtr = malloc ( phyloTreePtr->numberLeaves * sizeof( SetKeysPtr ) );
	if ( arrPtr == NULL ) { nomemMessage( "arrPtr" ); }

	for ( i = 0; i < phyloTreePtr->numberLeaves; i++ ){
		arrPtr[ i ] = malloc ( sizeof( SetKeys ) );
		if ( arrPtr[ i ] == NULL ) { nomemMessage( "arrPtr[i]" ); }

		arrPtr[ i ]->numKeys = 0;
		arrPtr[ i ]->condKeysPtrArray = malloc( maxNumKeys * sizeof( CondenseKeyPtr ) );
		if ( arrPtr[ i ]->condKeysPtrArray == NULL ) { nomemMessage( "arrPtr[i]->condKeysPtrArray" ); }

		for ( j = 0; j < maxNumKeys; j++ ) {
			arrPtr[ i ]->condKeysPtrArray[j] = malloc( sizeof( CondenseKey ) );
			if ( arrPtr[ i ]->condKeysPtrArray[j] == NULL ) { nomemMessage( "arrPtr[i]->condKeysPtrArray[j]" ); }
		}
	}

	/* for each genome, enumerate all their sequences  */
	for ( i = 0; i < phyloTreePtr->numberLeaves; i++ ) {
		
		/* find sequences for genome i */
		startFound = FALSE;
		sequenceComplete = FALSE;
		ikey = 0;
		for ( j = 0; j < phyloTreePtr->numberGenes - 1; j++ ) {
			current = phyloTreePtr->nodesPtrArray[ i ]->genome[ j ];
			next 	= phyloTreePtr->nodesPtrArray[ i ]->genome[ j + 1 ];

			/* check if there exists a sequence between current gene and next */
			if ( abs( next - current ) == 1 ) {
				if ( startFound == FALSE ) {
					start = current;
					end = next;
				}
				else {
					end = next;
				}
				startFound = TRUE;
				sequenceComplete = FALSE;
			}
			else { 
				if ( startFound == TRUE ) { sequenceComplete = TRUE; }
				startFound = FALSE;
			}

			/* if a complete sequences was found, then save its interval */
			if ( ( sequenceComplete == TRUE ) || 
				(startFound == TRUE && j == phyloTreePtr->numberGenes - 2) ) {

				/* change to positive if necessary */
				if ( start < 0 ) {
					int tmp = start;
					start = end * -1;
					end = tmp * -1;
				}
				
				sequenceComplete = FALSE;
				arrPtr[ i ]->condKeysPtrArray[ ikey ]->gene 		= start < end ? start : end;
				arrPtr[ i ]->condKeysPtrArray[ ikey ]->start 		= start;
				arrPtr[ i ]->condKeysPtrArray[ ikey ]->end 			= end;
				arrPtr[ i ]->condKeysPtrArray[ ikey ]->diff 		= abs( end - start );
				arrPtr[ i ]->condKeysPtrArray[ ikey ]->orientation	= start < end ? INCREASING : DECREASING;
				ikey++; 
			} 

		}
		arrPtr[ i ]->numKeys = ikey;
	}

	/* determine if a sequence exists in the other genomes (or is included in other sequence) */
	CondenseKey sequenceA, sequenceB;
	ikey = 0;
	for ( i = 0; i < phyloTreePtr->numberLeaves; i++ ) {
		for ( j = 0; j < arrPtr[ i ]->numKeys; j++ ) {

			sequenceA = *arrPtr[ i ]->condKeysPtrArray[ j ];
			/* find sequence A in the other genomes */
			for ( kk = 0; kk < phyloTreePtr->numberLeaves; kk++ ) {
				/* skip the case when the two genomes (i and kk) are the same */
				if ( i == kk ) 
					continue;

				/* find sequence A in genome kk */
				found = FALSE;
				for ( hh = 0; hh < arrPtr[ kk ]->numKeys; hh++ ) {	
					sequenceB = *arrPtr[ kk ]->condKeysPtrArray[ hh ];
					found = findIntersection(&sequenceA, &sequenceB, &start, &end);
					if ( found == TRUE ) { break; }
				}
				/* if intersection was found update sequence A */
				if ( found == TRUE ) {
					sequenceA.gene 			= start < end ? start : end;
					sequenceA.start 		= start;
					sequenceA.end 			= end;
					sequenceA.diff 			= abs( end - start );
					sequenceA.orientation	= start < end ? INCREASING : DECREASING;	
				}
				else {
					break; /* stop the search for sequence A */
				}
			}
			/* if intersection found add to the set of keys */
			if ( found == TRUE ) {
				addUniqueKey(sKeysPtr, &ikey, sequenceA.start, sequenceA.end);
			}
			
		}
	}
	sKeysPtr->numKeys = ikey;
	
	/* update the intervals of the sequences in the order of appearance */
	for ( i = 0; i < sKeysPtr->numKeys - 1; i++ ) {
		
		for ( j = i + 1; j < sKeysPtr->numKeys; j++ ) {

			if ( sKeysPtr->condKeysPtrArray[ i ]->gene < 
					sKeysPtr->condKeysPtrArray[ j ]->gene ) {

				start = sKeysPtr->condKeysPtrArray[ j ]->start - sKeysPtr->condKeysPtrArray[ i ]->diff; 
				end = sKeysPtr->condKeysPtrArray[ j ]->end - sKeysPtr->condKeysPtrArray[ i ]->diff;
				//update interval 
				sKeysPtr->condKeysPtrArray[ j ]->gene 			= start < end ? start : end;
				sKeysPtr->condKeysPtrArray[ j ]->start 			= start;
				sKeysPtr->condKeysPtrArray[ j ]->end 			= end;
				sKeysPtr->condKeysPtrArray[ j ]->diff 			= abs( end - start );
				sKeysPtr->condKeysPtrArray[ j ]->orientation	= start < end ? INCREASING : DECREASING;
			}
		}
	} 

	/* free memory */
	for ( i = 0; i < phyloTreePtr->numberLeaves; i++ ){
		for ( j = 0; j < maxNumKeys; j++ ) {
			free( arrPtr[ i ]->condKeysPtrArray[ j ] );
		}
		free( arrPtr[ i ]->condKeysPtrArray );
		free( arrPtr[ i ] );		
	}
	free( arrPtr );
}

/* Note: elements in the input sequence A and B must be positives */
static int findIntersection(CondenseKeyPtr sequenceAPtr, 
	CondenseKeyPtr sequenceBPtr, int *start, int *end) 
{
	int found = FALSE;
	/* case: both sequences have (increasing) orientation */
	if ( sequenceAPtr->orientation == INCREASING && 
					sequenceBPtr->orientation == INCREASING) {

		if ( isIncludedIn( sequenceAPtr->start, sequenceBPtr->start, sequenceBPtr->end ) == TRUE ) {
			( *start ) = sequenceAPtr->start;
			if ( isIncludedIn( sequenceAPtr->end, sequenceBPtr->start, sequenceBPtr->end ) == TRUE ) {
				( *end ) = sequenceAPtr->end;
				found = TRUE;
			}
			else if ( isIncludedIn( sequenceBPtr->end, sequenceAPtr->start, sequenceAPtr->end ) == TRUE ) {
				( *end ) = sequenceBPtr->end;
				found = TRUE;
			}
		}
		else if ( isIncludedIn( sequenceBPtr->start, sequenceAPtr->start, sequenceAPtr->end ) == TRUE ) {
			( *start ) = sequenceBPtr->start;
			if ( isIncludedIn( sequenceBPtr->end, sequenceAPtr->start, sequenceAPtr->end ) == TRUE ) {
				( *end ) = sequenceBPtr->end;
				found = TRUE;
			}
			else if ( isIncludedIn( sequenceAPtr->end, sequenceBPtr->start, sequenceBPtr->end ) ) {
				( *end ) = sequenceAPtr->end;
				found = TRUE;
			}
		} 
	}
	/* case: both sequences have (decreasing) orientation*/
	else if ( sequenceAPtr->orientation == DECREASING && 
					sequenceBPtr->orientation == DECREASING ) {

		if ( isIncludedIn( sequenceAPtr->start, sequenceBPtr->end, sequenceBPtr->start ) == TRUE ) {
			( *start ) = sequenceAPtr->start;
			if ( isIncludedIn( sequenceAPtr->end, sequenceBPtr->end, sequenceBPtr->start ) == TRUE ) {
				( *end ) = sequenceAPtr->end;
				found = TRUE;
			}
			else if ( isIncludedIn( sequenceBPtr->end, sequenceAPtr->end, sequenceAPtr->start ) == TRUE ) {
				( *end ) = sequenceBPtr->end;
				found = TRUE;
			}
		}
		else if ( isIncludedIn( sequenceBPtr->start, sequenceAPtr->end, sequenceAPtr->start ) == TRUE ) {
			( *start ) = sequenceBPtr->start;
			if ( isIncludedIn( sequenceBPtr->end, sequenceAPtr->end, sequenceAPtr->start ) == TRUE ) {
				( *end ) = sequenceBPtr->end;
				found = TRUE;
			}
			else if ( isIncludedIn( sequenceAPtr->end, sequenceBPtr->end, sequenceBPtr->start ) ) {
				( *end ) = sequenceAPtr->end;
				found = TRUE;
			}
		} 
	}
	else {
		/* the sequences dont have the same orientation,
			in the best case they just share one element, which is
			a invalid sequence because length must be > 1 */
	}

	/* start and end can not be the same */
	if ( found == TRUE && *start == *end ) { found = FALSE; }

	return found;
}

/* Note: start must be < than end */
static int isIncludedIn( int a, int start, int end) 
{
	if ( a >= start && a <= end ) {
		return TRUE;
	}
	else {
		return FALSE;
	} 
}

static void addUniqueKey(SetKeysPtr sKeysPtr, int *ikey, int start, int end)
{
	int i, found;

	found = FALSE;
	for ( i = 0; i < *ikey; i++ ) {
		if ( sKeysPtr->condKeysPtrArray[ i ]->start == start && 
					sKeysPtr->condKeysPtrArray[ i ]->end == end ) {
			found = TRUE;
			break;
		}
	} 

	if ( found == FALSE ) {
		sKeysPtr->condKeysPtrArray[ *ikey ]->gene 			= start < end ? start : end;
		sKeysPtr->condKeysPtrArray[ *ikey ]->start 			= start;
		sKeysPtr->condKeysPtrArray[ *ikey ]->end 			= end;
		sKeysPtr->condKeysPtrArray[ *ikey ]->diff 			= abs( end - start );
		sKeysPtr->condKeysPtrArray[ *ikey ]->orientation	= start < end ? INCREASING : DECREASING;
		( *ikey )++;				
	}
}



