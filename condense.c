/*
 ============================================================================
 Authors :
 	Jose Luis Soncco-Alvarez and Mauricio Ayala-Rincon
 	Group of Theory of Computation
 	Universidade de Brasilia (UnB) - Brazil
 ============================================================================
 */

#include <stdio.h>
#include <stdlib.h>

#include "condense.h"
#include "auxiliary.h"
#include "my_structs.h"
#include "tree.h"

static int findIntersection(CondenseKeyPtr sequenceAPtr, 
	CondenseKeyPtr sequenceBPtr, int *start, int *end);
static int isIncludedIn( int a, int start, int end); 
static void addUniqueKey( SetKeysPtr setkeysPtr, int *ikey, int start, int end );

void allocateMemoryForRawData( RawDatasetPtr rdatasetPtr )
{
	int i;
	//rdatasetPtr->numberGenomes; this value must be read before entering this function
	rdatasetPtr->numberGenes = 0;
	rdatasetPtr->numberChromosomesArray = malloc( rdatasetPtr->numberGenomes * sizeof( int ) );
	if ( rdatasetPtr->numberChromosomesArray == NULL ) { nomemMessage("rdatasetPtr->numberChromosomesArray"); }

	for ( i = 0; i < rdatasetPtr->numberGenomes; i++ ) {
		rdatasetPtr->numberChromosomesArray[ i ] = 0;
	}
	rdatasetPtr->rgenomes = NULL;
}

void allocateMemoryForRawGenomes( RawDatasetPtr rdatasetPtr )
{
	int i;

	rdatasetPtr->rgenomes = malloc( rdatasetPtr->numberGenomes * sizeof( RawGenomePtr ) );
	if ( rdatasetPtr->rgenomes == NULL ) { nomemMessage("rdatasetPtr->rgenomes"); }

	for ( i = 0; i < rdatasetPtr->numberGenomes; i++ ) {
		rdatasetPtr->rgenomes[ i ] = malloc( sizeof( RawGenome ) );
		if ( rdatasetPtr->rgenomes[ i ] == NULL ) { nomemMessage("rdatasetPtr->rgenomes[ i ]"); }

        rdatasetPtr->rgenomes[ i ]->organism = NULL;
        rdatasetPtr->rgenomes[ i ]->numberChromosomes = rdatasetPtr->numberChromosomesArray[ i ];
        rdatasetPtr->rgenomes[ i ]->chromosomeType = malloc( rdatasetPtr->numberChromosomesArray[ i ] * sizeof( int ) );
        if ( rdatasetPtr->rgenomes[ i ]->chromosomeType == NULL ) { 
            nomemMessage("rdatasetPtr->rgenomes[ i ]->chromosomeType"); 
        }
		rdatasetPtr->rgenomes[ i ]->numberElements = 
						rdatasetPtr->numberGenes + rdatasetPtr->numberChromosomesArray[ i ];
		rdatasetPtr->rgenomes[ i ]->genome = 
						malloc( rdatasetPtr->rgenomes[ i ]->numberElements * sizeof( int ) );
		if ( rdatasetPtr->rgenomes[ i ]->genome == NULL ) { 
            nomemMessage("rdatasetPtr->rgenomes[ i ]->genes"); 
        }
	}
}

void allocateMemoryForKeys( int numberGenes, SetKeysPtr setkeysPtr ) 
{
	int i,  maxNumKeys;

	maxNumKeys = numberGenes / 2;
	setkeysPtr->condkeyPtrArray = malloc( maxNumKeys * sizeof( CondenseKeyPtr ) );
	if ( setkeysPtr->condkeyPtrArray == NULL ) { 
		nomemMessage( "setPtr->condkeyPtrArray" ); 
	}
	for ( i = 0; i < maxNumKeys; i++ ) {
		setkeysPtr->condkeyPtrArray[i] = malloc( sizeof( CondenseKey ) );
		if ( setkeysPtr->condkeyPtrArray[i] == NULL ) { 
			nomemMessage( "setPtr->condkeyPtrArray[i]" ); 
		}
	}
}

void freeKeys( int numberGenes, SetKeysPtr setkeysPtr ) 
{
	int i,  maxNumKeys;

	maxNumKeys = numberGenes / 2;
	for ( i = 0; i < maxNumKeys; i++ ) {
		free( setkeysPtr->condkeyPtrArray[ i ] );
	}
	free( setkeysPtr->condkeyPtrArray );
}

void freeRawDataset( RawDatasetPtr rdatasetPtr )
{
	int i;
	
	free ( rdatasetPtr->numberChromosomesArray );
	for ( i = 0; i < rdatasetPtr->numberGenomes; i++ ) {
		free( rdatasetPtr->rgenomes[ i ]->genome );
        free( rdatasetPtr->rgenomes[ i ]->chromosomeType );
        free( rdatasetPtr->rgenomes[ i ]->organism );
		free( rdatasetPtr->rgenomes[ i ] );
	}
	free( rdatasetPtr->rgenomes );
}

void findSetKeysForCondensing( RawDatasetPtr rdatasetPtr, SetKeysPtr setkeysPtr )
{
	int i, j, kk, hh, ikey, maxNumKeys;
	int current, next, start, end, startFound, found, sequenceComplete;

	/* allocate memory for array of setkeys */
	SetKeysPtr *setkeysPtrArray; 
	maxNumKeys = rdatasetPtr->numberGenes / 2;
	setkeysPtrArray = malloc ( rdatasetPtr->numberGenomes * sizeof( SetKeysPtr ) );
	if ( setkeysPtrArray == NULL ) { nomemMessage( "setkeysPtrArray" ); }

	for ( i = 0; i < rdatasetPtr->numberGenomes; i++ ) {
		setkeysPtrArray[ i ] = malloc ( sizeof( SetKeys ) );
		if ( setkeysPtrArray[ i ] == NULL ) { nomemMessage( "setkeysPtrArray[i]" ); }

		setkeysPtrArray[ i ]->numberKeys = 0;
		setkeysPtrArray[ i ]->condkeyPtrArray = malloc( maxNumKeys * sizeof( CondenseKeyPtr ) );
		if ( setkeysPtrArray[ i ]->condkeyPtrArray == NULL ) { 
			nomemMessage( "setkeysPtrArray[i]->condkeyPtrArray" ); 
		}

		for ( j = 0; j < maxNumKeys; j++ ) {
			setkeysPtrArray[ i ]->condkeyPtrArray[j] = malloc( sizeof( CondenseKey ) );
			if ( setkeysPtrArray[ i ]->condkeyPtrArray[j] == NULL ) { 
				nomemMessage( "setkeysPtrArray[i]->condkeyPtrArray[j]" ); 
			}
		}
	}

	/* for each rawGenome, enumerate all their sequences of consecutive numbers */
	for ( i = 0; i < rdatasetPtr->numberGenomes; i++ ) {
		
		/* find sequences for rawGenome i */
		startFound = FALSE;
		sequenceComplete = FALSE;
		ikey = 0;
		for ( j = 0; j < rdatasetPtr->rgenomes[ i ]->numberElements; j++ ) {
			/* AVOID case of reaching the simbol $ or @ represented by SPLIT */
			if ( rdatasetPtr->rgenomes[ i ]->genome[ j ] == SPLIT ) {
				continue;
			}

			current = rdatasetPtr->rgenomes[ i ]->genome[ j ];
			next 	= rdatasetPtr->rgenomes[ i ]->genome[ j + 1 ];

			/* check if there exists a sequence between current and next gene */
			if ( next != SPLIT && abs( next - current ) == 1 ) {
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
            if ( sequenceComplete == TRUE ) {
                //|| (startFound == TRUE && j == numberGenes - 2) ) {

				/* change to positive if necessary */
				if ( start < 0 ) {
					int tmp = start;
					start = end * -1;
					end = tmp * -1;
				}
				
				sequenceComplete = FALSE;
				setkeysPtrArray[ i ]->condkeyPtrArray[ ikey ]->gene 		= start < end ? start : end;
				setkeysPtrArray[ i ]->condkeyPtrArray[ ikey ]->start 		= start;
				setkeysPtrArray[ i ]->condkeyPtrArray[ ikey ]->end 			= end;
				setkeysPtrArray[ i ]->condkeyPtrArray[ ikey ]->diff 		= abs( end - start );
				setkeysPtrArray[ i ]->condkeyPtrArray[ ikey ]->orientation	= start < end ? INCREASING : DECREASING;
				ikey++; 
			} 

		}
		setkeysPtrArray[ i ]->numberKeys = ikey;
	}

	/* determine if a sequence exists in the other rawGenomes (or is included in other sequence) */
	CondenseKey sequenceA, sequenceB;
	ikey = 0;
	for ( i = 0; i < rdatasetPtr->numberGenomes; i++ ) {
		for ( j = 0; j < setkeysPtrArray[ i ]->numberKeys; j++ ) {

			sequenceA = *setkeysPtrArray[ i ]->condkeyPtrArray[ j ];
			/* find sequence A in the other rawGenomes */
			for ( kk = 0; kk < rdatasetPtr->numberGenomes; kk++ ) {
				/* skip the case when the two rawGenomes (i and kk) are the same */
				if ( i == kk ) 
					continue;

				/* find sequence A in rawGenome kk */
				found = FALSE;
				for ( hh = 0; hh < setkeysPtrArray[ kk ]->numberKeys; hh++ ) {	
					sequenceB = *setkeysPtrArray[ kk ]->condkeyPtrArray[ hh ];
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
				addUniqueKey(setkeysPtr, &ikey, sequenceA.start, sequenceA.end);
			}
			
		}
	}
	setkeysPtr->numberKeys = ikey;
	
	/* update the intervals of the sequences */
	for ( i = 0; i < setkeysPtr->numberKeys - 1; i++ ) {
		
		for ( j = i + 1; j < setkeysPtr->numberKeys; j++ ) {

			if ( setkeysPtr->condkeyPtrArray[ i ]->gene < 
					setkeysPtr->condkeyPtrArray[ j ]->gene ) {
                //re-calculate "start" and "end" of key j, because an
                //application of key i will modify elements of key j
                
				start = setkeysPtr->condkeyPtrArray[ j ]->start - setkeysPtr->condkeyPtrArray[ i ]->diff; 
				end = setkeysPtr->condkeyPtrArray[ j ]->end - setkeysPtr->condkeyPtrArray[ i ]->diff;
				//update interval 
				setkeysPtr->condkeyPtrArray[ j ]->gene 			= start < end ? start : end;
				setkeysPtr->condkeyPtrArray[ j ]->start 		= start;
				setkeysPtr->condkeyPtrArray[ j ]->end 			= end;
				setkeysPtr->condkeyPtrArray[ j ]->diff 			= abs( end - start );
				setkeysPtr->condkeyPtrArray[ j ]->orientation	= start < end ? INCREASING : DECREASING;
			}
		}
	} 

	/* free memory */
	for ( i = 0; i < rdatasetPtr->numberGenomes; i++ ){
		for ( j = 0; j < maxNumKeys; j++ ) {
			free( setkeysPtrArray[ i ]->condkeyPtrArray[ j ] );
		}
		free( setkeysPtrArray[ i ]->condkeyPtrArray );
		free( setkeysPtrArray[ i ] );		
	}
	free( setkeysPtrArray );
}


void condenseGenomes( RawDatasetPtr rdatasetPtr, SetKeysPtr setkeysPtr )
{
	int i, j, k, iwrite, otherExtreme;	

	for ( i = 0; i < setkeysPtr->numberKeys; i++ ) {

		// determine the other extreme of the interval 
		otherExtreme = setkeysPtr->condkeyPtrArray[i]->start;
		if ( setkeysPtr->condkeyPtrArray[i]->gene == 
				setkeysPtr->condkeyPtrArray[i]->start) {
			otherExtreme = setkeysPtr->condkeyPtrArray[i]->end;
		} 

		// condense leaves using key i 
		for ( j = 0; j < rdatasetPtr->numberGenomes; j++ ) {
			
			iwrite = 0;
			for ( k = 0; k < rdatasetPtr->rgenomes[ j ]->numberElements; k++ ) {

				if ( abs( rdatasetPtr->rgenomes[ j ]->genome[ k ] ) <= 
								setkeysPtr->condkeyPtrArray[ i ]->gene ) {

					rdatasetPtr->rgenomes[ j ]->genome[ iwrite ] = 
								rdatasetPtr->rgenomes[ j ]->genome[ k ];
					iwrite++;
				}
				else if ( abs( rdatasetPtr->rgenomes[ j ]->genome[ k ] ) > 
							otherExtreme ) {

					// update the gen by subtracting (or adding) the difference 
					if ( rdatasetPtr->rgenomes[ j ]->genome[ k ] > 0 ) {

						rdatasetPtr->rgenomes[ j ]->genome[ iwrite ] =
								rdatasetPtr->rgenomes[ j ]->genome[ k ] - 
								setkeysPtr->condkeyPtrArray[ i ]->diff;
					}
					else { // rdatasetPtr->rgenomes[ j ]->rawGenome[ k ] < 0

						rdatasetPtr->rgenomes[ j ]->genome[ iwrite ] =
								rdatasetPtr->rgenomes[ j ]->genome[ k ] + 
								setkeysPtr->condkeyPtrArray[ i ]->diff;
					}
					iwrite++;
				}

			} 
			rdatasetPtr->rgenomes[ j ]->numberElements = iwrite;
		}
	}	

	/* update number of genes */
	for ( i = 0; i < setkeysPtr->numberKeys; i++ ) {
		rdatasetPtr->numberGenes = rdatasetPtr->numberGenes - setkeysPtr->condkeyPtrArray[ i ]->diff;
	} 
}

void showKeys( SetKeysPtr setkeysPtr ) 
{
	int i;
	for ( i = 0; i < setkeysPtr->numberKeys; i++ ) {
		printf( "*gene(%d) --> start %d, end %d, diff %d\n", 
			setkeysPtr->condkeyPtrArray[ i ]->gene, 
			setkeysPtr->condkeyPtrArray[ i ]->start, 
			setkeysPtr->condkeyPtrArray[ i ]->end,
			setkeysPtr->condkeyPtrArray[ i ]->diff );
	} 
}

void unpackCondensedGenomes( RawDatasetPtr rdatasetPtr, SetKeysPtr setkeysPtr )
{
	int i, j, k, iseq, diff, gene;

	for ( i = setkeysPtr->numberKeys - 1; i >= 0; i-- ) {
		diff = setkeysPtr->condkeyPtrArray[ i ]->diff;
		// reverse condensation of nodes using key i 
		for ( j = 0; j < rdatasetPtr->numberGenomes; j++ ) {
			// find gene that represents a sequence 
			for ( k = 0; k < rdatasetPtr->rgenomes[ j ]->numberElements; k++ ) {
				if ( abs( rdatasetPtr->rgenomes[ j ]->genome[ k ] ) > 
							setkeysPtr->condkeyPtrArray[ i ]->gene ) {

					if ( rdatasetPtr->rgenomes[ j ]->genome[ k ] > 0 ) {
						rdatasetPtr->rgenomes[ j ]->genome[ k ] = 
							rdatasetPtr->rgenomes[ j ]->genome[ k ] + diff;
					}
					else { // < 0
						rdatasetPtr->rgenomes[ j ]->genome[ k ] = 
							rdatasetPtr->rgenomes[ j ]->genome[ k ] - diff;
					}
				}
				else if ( abs( rdatasetPtr->rgenomes[ j ]->genome[ k ] ) == 
								setkeysPtr->condkeyPtrArray[ i ]->gene ) {
					iseq = k;
					break;
				}
			}

			// make space by shifting to the right 
			for ( k = rdatasetPtr->rgenomes[ j ]->numberElements - 1; k > iseq; k-- ) {
				if ( abs( rdatasetPtr->rgenomes[ j ]->genome[ k ] ) > 
							setkeysPtr->condkeyPtrArray[ i ]->gene ) {

					if ( rdatasetPtr->rgenomes[ j ]->genome[ k ] > 0 ) {
						rdatasetPtr->rgenomes[ j ]->genome[ k + diff ] = 
							rdatasetPtr->rgenomes[ j ]->genome[ k ] + diff;
					}
					else{ // < 0
						rdatasetPtr->rgenomes[ j ]->genome[ k + diff ] = 
							rdatasetPtr->rgenomes[ j ]->genome[ k ] - diff;
					}
				}
				else { // < (includes case when genome[k] is zero)
					rdatasetPtr->rgenomes[ j ]->genome[ k + diff ] = 
							rdatasetPtr->rgenomes[ j ]->genome[ k ];
				}
			}

			// insert the sequence  
			if ( rdatasetPtr->rgenomes[ j ]->genome[ iseq ] > 0 ) { 
				gene = setkeysPtr->condkeyPtrArray[ i ]->start; 
			}
			else { // < 0
				gene = -1 * setkeysPtr->condkeyPtrArray[ i ]->end; 	
			}
			
			if ( setkeysPtr->condkeyPtrArray[ i ]->orientation == INCREASING ) {
				for ( k = iseq; k <= iseq + diff; k++ ) {
					rdatasetPtr->rgenomes[ j ]->genome[ k ] = gene;
					gene = gene + 1; 	
				}
			}
			else { // setkeysPtr->condkeyPtrArray[i]->orientation == DECREASING
				for ( k = iseq; k <= iseq + diff; k++ ) {
					rdatasetPtr->rgenomes[ j ]->genome[ k ] = gene;
					gene = gene - 1; 	
				}
			}
			rdatasetPtr->rgenomes[ j ]->numberElements = 
						rdatasetPtr->rgenomes[ j ]->numberElements + diff;
		}//end-for
		rdatasetPtr->numberGenes = rdatasetPtr->numberGenes + diff;

	} //end-for
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

static void addUniqueKey(SetKeysPtr setkeysPtr, int *ikey, int start, int end)
{
	int i, found;

	found = FALSE;
	for ( i = 0; i < *ikey; i++ ) {
		if ( setkeysPtr->condkeyPtrArray[ i ]->start == start && 
					setkeysPtr->condkeyPtrArray[ i ]->end == end ) {
			found = TRUE;
			break;
		}
	} 

	if ( found == FALSE ) {
		setkeysPtr->condkeyPtrArray[ *ikey ]->gene 			= start < end ? start : end;
		setkeysPtr->condkeyPtrArray[ *ikey ]->start 		= start;
		setkeysPtr->condkeyPtrArray[ *ikey ]->end 			= end;
		setkeysPtr->condkeyPtrArray[ *ikey ]->diff 			= abs( end - start );
		setkeysPtr->condkeyPtrArray[ *ikey ]->orientation	= start < end ? INCREASING : DECREASING;
		( *ikey )++;				
	}
}