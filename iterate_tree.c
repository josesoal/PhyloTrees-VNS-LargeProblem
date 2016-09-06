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

#include "random.h"
#include "auxiliary.h"
#include "iterate_tree.h"
#include "tree.h" 
#include "median_solvers.h"
#include "queue.h"
#include "invdist.h" 		/* inversion distance  */
#include "binencode.h"		/* breakpoint distance */
#include "dcjdist.h"		/* DCJ distance */

static int labelOptimizeTree_Blanchette( 
			TreePtr phyloTreePtr, int initialize, ParametersPtr paramsPtr );
static int labelOptimizeTree_KovacDCJ( 
			TreePtr phyloTreePtr, int initialize,  ParametersPtr paramsPtr );
static void generateCandidates( TreePtr phyloTreePtr, TreeNodePtr nodePtr,
						CandidatePtr **candidateMatrix, int *numCandidates, 
						int *maxCandidates, ParametersPtr paramsPtr );
static void copyCurrentNodeIntoCandidate( TreePtr phyloTreePtr
				CandidatePtr candPtr, TreeNodePtr nodePtr, int copyInverse );


static int labelOptimizeTree_GreedyCandidatesDCJ( 
			TreePtr phyloTreePtr, int initialize,  ParametersPtr paramsPtr ); 	
static void initializeTreeUsingDFS( TreePtr phyloTreePtr, 
			TreeNodePtr nodePtr, ParametersPtr paramsPtr );
static void initializeTreeWithDescendants( TreePtr phyloTreePtr, 
			TreeNodePtr nodePtr, ParametersPtr paramsPtr, int orientation );
static void initializeTreeWithDescendantsRandomly( 
			TreePtr phyloTreePtr, TreeNodePtr nodePtr );
static TreeNodePtr findNearestNeighbor( TreePtr phyloTreePtr, 
			TreeNodePtr nodePtr, int idNeighbor1, int idNeighbor2 );
static void callSolverAndLabelNode( TreePtr phyloTreePtr, TreeNodePtr nodePtr, 
	TreeNodePtr node1Ptr, TreeNodePtr node2Ptr, TreeNodePtr node3Ptr,
	struct adj_struct *adjacencyList, enum medianSolvers solver, int circular );
static void iterateTreeUsingDFS( 
	TreePtr phyloTreePtr, TreeNodePtr nodePtr, ParametersPtr paramsPtr ); 
static void improveTreebyCandidatesDCJ( TreePtr phyloTreePtr, 
	TreeNodePtr nodePtr, ParametersPtr paramsPtr ); 
static void allocateMemForCandidate( TreePtr phyloTreePtr, CandidatePtr *candPtr );
static void freeCandidate( TreePtr phyloTreePtr, CandidatePtr *candPtr ); 
static void verifyPenalization( TreePtr phyloTreePtr, 
	CandidatePtr *candidatePtrArray, int *h, ParametersPtr paramsPtr );

int labelOptimizeTree( TreePtr phyloTreePtr, ParametersPtr paramsPtr )
{
	int score;
	int initialize; /* initialize label of internal nodes ? */

	score = 0;
	if ( paramsPtr->distanceType == DCJ_DIST ) {
		if ( paramsPtr->opt == GREEDY_CANDIDATES ) {
			initialize = TRUE;
			score = labelOptimizeTree_GreedyCandidatesDCJ( 
									phyloTreePtr, initialize, paramsPtr ); 	
		}
		else if ( paramsPtr->opt == KOVAC ) {
			initialize = TRUE;
			score = labelOptimizeTree_KovacDCJ( 
									phyloTreePtr, initialize, paramsPtr);
		}
		else if ( paramsPtr->opt == HERENCSAR ) {
			//initialize = TRUE;
			//score = labelOptimizeTree_HerencsarDCJ( phyloTreePtr, initialize, paramsPtr);
			fprintf( stderr, " stderr: this function is not implemented\n" );
        	exit( EXIT_FAILURE );
		}
		else {
			fprintf( stderr, " stderr: incorrect label-optimize method for DCJ distance.\n" );
        	exit( EXIT_FAILURE );
		}
	}
	else { // paramsPtr->distanceType == INVERSION_DIST
		if ( paramsPtr->opt == BLANCHETTE ) {
			initialize = TRUE; 
			score = labelOptimizeTree_Blanchette( 
										phyloTreePtr, initialize, paramsPtr);
		}
		else {
			fprintf( stderr, " stderr: incorrect label-optimize method for reversal distance.\n" );
        	exit( EXIT_FAILURE );
		}	 
	}
	return score;
}

/* [OPTIMIZER for reversal distance] 
	algorithm for labeling and optimizing the score of a tree by 
	iteratively appliying the median of 3 genomes. */
/* NOTE: this algorithm (originally called optimize_tree) was implemented as 
	proposed by Blanchette (1997) in the paper "Breakpoint Phylogenies" */
static int labelOptimizeTree_Blanchette( 
	TreePtr phyloTreePtr, int initialize, ParametersPtr paramsPtr )
{
	int i, improved, score, newScore;
	Tree tempTree;

	/* allocate memory for temporary tree */
	tempTree.numberLeaves = phyloTreePtr->numberLeaves;
	tempTree.numberGenes = phyloTreePtr->numberGenes;
	allocateMemoryForNodes( &tempTree, paramsPtr );//--method from tree.c

	if ( initialize == TRUE ) {
		/* set leaves nodes as extremities, do the contrary for internal nodes */
		for ( i = 0; i < phyloTreePtr->numberLeaves; i++ ) {
			phyloTreePtr->nodesPtrArray[ i ]->extremity = TRUE;
		}
		for ( i = phyloTreePtr->numberLeaves; i < phyloTreePtr->numberNodes; i++ ) {
			phyloTreePtr->nodesPtrArray[ i ]->extremity = FALSE;
		}
	
		/* initialize label of internal nodes with genomes */
		initializeTreeUsingDFS( phyloTreePtr, phyloTreePtr->startingNodePtr, paramsPtr );
	}

	/* iterate the median over tree whenever the score can be improved */
	improved = TRUE;
	score = scoreTree( phyloTreePtr, phyloTreePtr->startingNodePtr->rightDescPtr, paramsPtr );

	/* improve tree by iterating */
	while ( improved ) {
		copyTreeInto( &tempTree, phyloTreePtr, FALSE, paramsPtr );/* make a copy of current tree */
		iterateTreeUsingDFS( phyloTreePtr, phyloTreePtr->startingNodePtr->rightDescPtr, paramsPtr );
		newScore = scoreTree( phyloTreePtr, phyloTreePtr->startingNodePtr->rightDescPtr, paramsPtr );
		
		if ( newScore < score ) {
			//printf("\t score: %d, new score: %d\n",score,newScore);
			score = newScore;
		}
		else { 
			improved = FALSE;	
			copyTreeInto( phyloTreePtr, &tempTree, FALSE, paramsPtr ); /* recover last tree */	
		}
	}

	freeTree( &tempTree, paramsPtr );//--method from tree.c

	return score;
}

/* iterate the tree using a depth first search (DFS) approach.
NOTE: in the first call of this procedure, the second argument 
			must be the right descendant of the starting node */
static void iterateTreeUsingDFS( TreePtr phyloTreePtr, TreeNodePtr nodePtr, ParametersPtr paramsPtr )
{
	int a, b, c, i, score, newScore, *tempGenome;
	struct adj_struct *adjacencyList; /* from structs.h */

	if ( nodePtr->type == INTERNAL_NODE ) {
		iterateTreeUsingDFS( phyloTreePtr, nodePtr->leftDescPtr, paramsPtr );
		iterateTreeUsingDFS( phyloTreePtr, nodePtr->rightDescPtr, paramsPtr );
		
		adjacencyList = 
		malloc ( ( 2*phyloTreePtr->numberGenes + 1 ) * sizeof(struct adj_struct) );
    	if ( adjacencyList == NULL ) nomemMessage( "adjacencyList" );

    	tempGenome = malloc( phyloTreePtr->numberGenes * sizeof( int ) );
    	if ( tempGenome == NULL ) nomemMessage( "tempGenome" );

    	for ( i = 0; i < phyloTreePtr->numberGenes; i++ ) {
    		tempGenome[ i ] = nodePtr->genome[ i ];
    	}

    	score = nodePtr->edgeWeight + 
				nodePtr->leftDescPtr->edgeWeight + 
				nodePtr->rightDescPtr->edgeWeight;

		/* use the 3 neighbors to nodePtr: ancestor, left desc, and right desc
		* to create a TSP instance, and call a TSP solver to label nodePtr*/
		createTspInstance( phyloTreePtr->numberGenes, nodePtr->ancestorPtr->genome, 
			nodePtr->leftDescPtr->genome, nodePtr->rightDescPtr->genome, adjacencyList );//--method from median_solvers.c
		callSolverAndLabelNode( phyloTreePtr, nodePtr, nodePtr->ancestorPtr, 
			nodePtr->leftDescPtr, nodePtr->rightDescPtr, adjacencyList, paramsPtr->solver, paramsPtr->circular );

		a = calculateDistance( nodePtr->genome, 
			nodePtr->ancestorPtr->genome, phyloTreePtr->numberGenes, paramsPtr ); 
		b = calculateDistance( nodePtr->genome, 
			nodePtr->leftDescPtr->genome, phyloTreePtr->numberGenes, paramsPtr ); 
		c = calculateDistance( nodePtr->genome, 
			nodePtr->rightDescPtr->genome, phyloTreePtr->numberGenes, paramsPtr );

		newScore = a + b + c;

		if ( newScore < score ) {
			nodePtr->edgeWeight = a; 
			nodePtr->leftDescPtr->edgeWeight = b;
			nodePtr->rightDescPtr->edgeWeight = c;
		}
		else { //recover last node
			for ( i = 0; i < phyloTreePtr->numberGenes; i++ ) {
    			nodePtr->genome[ i ] = tempGenome[ i ];
    		}
		}
			
		free( adjacencyList );	
		free( tempGenome );			
	}
}

/* [OPTIMIZER for DCJ distance] 
	algorithm for labeling and optimizing the score of a tree by 
	generating candidates (genomes) for all internal nodes using a dynamic
	programming approach */
/* NOTE: Implemented as proposed by Kovac (2011) in the paper "A Practical 
	Algorithm for Ancestral Rearrangement Reconstruction" */
static int labelOptimizeTree_KovacDCJ( 
			TreePtr phyloTreePtr, int initialize,  ParametersPtr paramsPtr )
{
	int i, j;
	int tempTree;
	int *numCandidates; /* array of num. candidates for each internal node */
	int *maxCandidates; /* max. candidates allowed for each internal node */ 
	int **scoreMatrix;
	CandidatePtr **candidateMatrix;
	int score, newScore, improved, numInternalNodes;

	/* allocate memory for temporal tree */
	tempTree.numberLeaves = phyloTreePtr->numberLeaves;
	tempTree.numberGenes = phyloTreePtr->numberGenes;
	allocateMemoryForNodes( &tempTree, paramsPtr );
	/* allocate memory (just pointers) for matrix of scores and candidates */
	numInternalNodes = phyloTreePtr->numberNodes - phyloTreePtr->numberLeaves;
	scoreMatrix = malloc( numInternalNodes * sizeof( int * ) );
	if ( scoreMatrix == NULL ) 
		nomemMessage( "scoreMatrix" );  
	candidateMatrix = malloc ( numInternalNodes * sizeof( CandidatePtr * ) );
	if ( candidateMatrix == NULL ) 
		nomemMessage( "candidateMatrix" );
	/* allocate memory for arrays that saves number of candidates */ 
	numCandidates = malloc( numInternalNodes * sizeof( int * ) );
	if ( numCandidates == NULL )
		nomemMessage( "numCandidates" );
	maxCandidates = malloc( numInternalNodes * sizeof( int * ) );
	if ( maxCandidates == NULL )
		nomemMessage( "maxCandidates" );
	/*for ( i = 0; i < numInternalNodes; i++ ) {
		scoreMatrix[ i ] = malloc( numCandidates * sizeof( int ) );
		if ( scoreMatrix[ i ] == NULL ) nomemMessage( "scoreMatrix[ i ]" ); 
	}*/

	/*for ( i = 0; i < numInternalNodes; i++ ) {
		candidateMatrix[ i ] = 
							malloc( numCandidates * sizeof ( CandidatePtr ) );
		if ( candidateMatrix[ i ] == NULL ) 
			nomemMessage( "candidateMatrix[ i ]" );
	}*/
		
	/* allocate memory for candidates genomes */
	/*for ( i = 0; i < numInternalNodes; i++) {
		for (j = 0; j < numCandidates; j++) {
			candidateMatrix[ i ][ j ] =  malloc( sizeof( struct candidate ) );
			if ( candidateMatrix[ i ][ j ] == NULL ) nomemMessage( "candidateMatrix[ i ][ j ]" );

			candidateMatrix[ i ][ j ]->genome = malloc( phyloTreePtr->numberGenes * sizeof( int ) );
			if ( candidateMatrix[ i ][ j ]->genome == NULL ) nomemMessage( "candidateMatrix[ i ][ j ]" );
		}
	}*/
	/* initialize tree */
	if ( initialize == TRUE ) {
		if ( paramsPtr->problem == LARGE_PHYLOGENY ) {
			initializeTreeWithDescendants( phyloTreePtr, 
			phyloTreePtr->startingNodePtr->rightDescPtr, paramsPtr, GO_LEFT );
		}
		else { // SMALL_PHYLOGENY
			initializeTreeWithDescendantsRandomly( 
				phyloTreePtr, phyloTreePtr->startingNodePtr->rightDescPtr );
		}
	}

	/* label and optimize tree */
	improved = TRUE;
	score = scoreTree( phyloTreePtr, 
					phyloTreePtr->startingNodePtr->rightDescPtr, paramsPtr );
	while ( improved == TRUE ) {
		copyTreeInto( &tempTree, phyloTreePtr, FALSE, paramsPtr );
		generateCandidates( phyloTreePtr, 
			phyloTreePtr->startingNodePtr->rightDescPtr, 
			candidateMatrix, numCandidates, maxCandidates, paramsPtr )
		calculateScoreCandidates( phyloTreePtr, 
			phyloTreePtr->startingNodePtr->rightDescPtr, 
			scoreMatrix, candidateMatrix, numCandidates, paramsPtr );//--TOMOD
		selectBestNeighborhood( phyloTreePtr, 
			phyloTreePtr->startingNodePtr->rightDescPtr, 
			scoreMatrix, candidateMatrix, numCandidates, paramsPtr );//--TOMOD

		/* calculate new score of the tree */
		newScore = scoreTree( phyloTreePtr, phyloTreePtr->startingNodePtr->rightDescPtr, paramsPtr );
		
		if ( newScore < score ) {
			score = newScore;
		}
		else {
			improved = FALSE;
			/* recover last state of tree */
			copyTreeInto( phyloTreePtr, &tempTree, FALSE, paramsPtr ); 
		}
	}

	/* free memory (--TOMOD) */
	for ( i = 0; i < numInternalNodes; i++ ) {
		for (j = 0; j < numCandidates; j++) {
			free( candidateMatrix[ i ][ j ]->genome );
			free( candidateMatrix[ i ][ j ] );
		}
		free( scoreMatrix[ i ] );
		free( candidateMatrix[ i ] );
		free( tempGenomes[ i ] );
	}
	free( scoreMatrix );
	free( candidateMatrix );
	free( tempGenomes ); 

	return score;
}

/* [OPTIMIZER for DCJ distance] 
	algorithm for labeling and optimizing the score of a tree by 
	generating candidates (genomes) for all internal nodes using a dynamic
	programming approach */
/* NOTE: Implemented as proposed in "An Improved Algorithm for Ancestral 
	Gene Order Reconstruction" by Herencsar */
/*static*/ int labelOptimizeTree_HerencsarDCJ( 
			TreePtr phyloTreePtr, int initialize,  ParametersPtr paramsPtr )
{
	/* 
	*
	* NOT IMMPLEMENTED */
	return 0;
}

/* NOTE: in the first call of this procedure, the second argument 
			must be the right descendant of the starting node */
static void generateCandidates( TreePtr phyloTreePtr, TreeNodePtr nodePtr,
						CandidatePtr **candidateMatrix, int *numCandidates, 
						int *maxCandidates, ParametersPtr paramsPtr ) 
{
	int i, j, k, h, index, numA, numT;
	int x, y, xpos, ypos;

	if ( nodePtr->type == INTERNAL_NODE ) { 
		generateCandidates( phyloTreePtr, nodePtr->leftDescPtr, 
			candidateMatrix, numCandidates, distanceType, circular );
		generateCandidates( phyloTreePtr, nodePtr->rightDescPtr, 
			candidateMatrix, numCandidates, distanceType, circular );

		/* re-calculate index to be zero-based */
		index = nodePtr->index - phyloTreePtr->numberLeaves; 

		/* calculate max number of candidates for the current node */
		numA = 0; //number of adjacencies
		numT = 0; //number of telomeres
		for ( i = 0; i < nodePtr->numPointsDCJ; i++ ) {
			if ( nodePtr->genomeDCJ[ i ]->type == ADJACENCY )
				numA++;
			else //nodePtr->genomeDCJ[ i ]->type == TELOMERE
				numT++;
		}
		maxCandidates[ index ] = 	
			numA + 							// for inv. of case (3)
			2 * numA * ( numA - 1 ) / 2 + 	// 2 * C(numA, 2), for case (1)
			2 * numA * numT +				// for case (2)
			numT / 2;						// for case (3)

		/* allocate memory for array of candidates for nodePtr */
		candidateMatrix[ index ] = 
			malloc( maxCandidates[ index ] * sizeof( candidatePtr ) );
		if ( candidateMatrix[ index ] = NULL )
			nomemMessage("candidateMatrix[ index ]");

		/* generate candidates for inverse of case (3) */
		h = 0;
		for ( i = 0; i < nodePtr->numPointsDCJ; i++ ) {
			if ( nodePtr->genomeDCJ[ i ]->type == ADJACENCY ) {
				/* calculate position for checking inverses genomes */
				x = nodePtr->genomeDCJ[ i ]->x;
				y = nodePtr->genomeDCJ[ i ]->y;
				xpos = x > 0 ? 2 * x - 1 : 2 * abs( x ) - 2;
				ypos = y > 0 ? 2 * y - 1 : 2 * abs( y ) - 2;
				/* if adjacency (i->x, i->y) is also a adjacency in ancestor, 
					left desc., or right desc. do nothing */
				if ( nodePtr->ancestorPtr->inverseDCJ[ xpos ] == 
							nodePtr->ancestorPtr->inverseDCJ[ ypos ] && 
					nodePtr->leftDescPtr->inverseDCJ[ xpos ] == 
							nodePtr->leftDescPtr->inverseDCJ[ ypos ] && 
					nodePtr->rightDescPtr->inverseDCJ[ xpos ] == 
							nodePtr->rightDescPtr->inverseDCJ[ ypos ] ) {
					continue;
				}

				allocateMemForCandidate( 
							phyloTreePtr, &candidateMatrix[ index ][ h ] );

				copyCurrentNodeIntoCandidate( phyloTreePtr,
							candidateMatrix[ index ][ h ], nodePtr, FALSE );
				
				/* apply DCJ on copy and generate candidate */
				applyDCJ( candidateMatrix[ index ][ h ], 
						&candidateMatrix[ index ][ h ]->numPointsDCJ, 
						i, i, TRUE );//--from dcjdist.c
			
				calculateInverseGenome( 
					candidateMatrix[ index ][ h ]->genomeDCJ, 
					candidateMatrix[ index ][ h ]->numPointsDCJ, 
					candidateMatrix[ index ][ h ]->inverseDCJ );//--from dcjdist.c

				verifyPenalization( phyloTreePtr, 
									candidateMatrix[ index ], &h, paramsPtr );
			}
		}

		/* generate candidates by choosing two different indices: case (1), (2), and (3) */ 
		for ( i = 0; i < nodePtr->numPointsDCJ - 1; i++ ) {
			for ( j = i + 1; j < nodePtr->numPointsDCJ; j++ ) {
				/* calculate position for checking inverses genomes */
				x = nodePtr->genomeDCJ[ i ]->x;
				y = nodePtr->genomeDCJ[ i ]->y;
				xpos = x > 0 ? 2 * x - 1 : 2 * abs( x ) - 2;
				ypos = y > 0 ? 2 * y - 1 : 2 * abs( y ) - 2;

				/* if adjacency (i->x, i->y) is also a adjacency in ancestor, 
					left desc., or right desc. do nothing */
				if ( nodePtr->ancestorPtr->inverseDCJ[ xpos ] == nodePtr->ancestorPtr->inverseDCJ[ ypos ] && 
					nodePtr->leftDescPtr->inverseDCJ[ xpos ] == nodePtr->leftDescPtr->inverseDCJ[ ypos ] && 
					nodePtr->rightDescPtr->inverseDCJ[ xpos ] == nodePtr->rightDescPtr->inverseDCJ[ ypos ] ) {
					continue;
				}

				/* calculate position for checking inverses genomes */
				x = nodePtr->genomeDCJ[ j ]->x;
				y = nodePtr->genomeDCJ[ j ]->y;
				xpos = x > 0 ? 2 * x - 1 : 2 * abs( x ) - 2;
				ypos = y > 0 ? 2 * y - 1 : 2 * abs( y ) - 2;

				/* if adjacency (j->x,j->y) is also a adjacency in ancestor, 
					left desc., or right desc. do nothing */
				if ( nodePtr->ancestorPtr->inverseDCJ[ xpos ] == nodePtr->ancestorPtr->inverseDCJ[ ypos ] && 
					nodePtr->leftDescPtr->inverseDCJ[ xpos ] == nodePtr->leftDescPtr->inverseDCJ[ ypos ] && 
					nodePtr->rightDescPtr->inverseDCJ[ xpos ] == nodePtr->rightDescPtr->inverseDCJ[ ypos ] ) {
					continue;
				}

				allocateMemForCandidate( 
							phyloTreePtr, &candidateMatrix[ index ][ h ] );

				copyCurrentNodeIntoCandidate( phyloTreePtr,
							candidateMatrix[ index ][ h ], nodePtr, FALSE );

				/* apply DCJ on copy and generate candidate */
				applyDCJ( candidateMatrix[ index ][ h ]->genomeDCJ, 
						&candidateMatrix[ index ][ h ]->numPointsDCJ, 
						i, j, TRUE );//--from dcjdist.c
				
				calculateInverseGenome( 
					candidateMatrix[ index ][ h ]->genomeDCJ, 
					candidateMatrix[ index ][ h ]->numPointsDCJ, 
					candidateMatrix[ index ][ h ]->inverseDCJ );//--from dcjdist.c

				verifyPenalization( phyloTreePtr, 
									candidateMatrix[ index ], &h, paramsPtr ); 

				/* do the second FORM in case one of the point is a adjacency */
				if ( nodePtr->genomeDCJ[ i ]->type == ADJACENCY || 
						nodePtr->genomeDCJ[ j ]->type == ADJACENCY ) {

					allocateMemForCandidate( 
							phyloTreePtr, &candidateMatrix[ index ][ h ] );

					copyCurrentNodeIntoCandidate( phyloTreePtr,
							candidateMatrix[ index ][ h ], nodePtr, FALSE );

					/* apply DCJ on copy and generate candidate */
					applyDCJ( candidateMatrix[ index ][ h ]->genomeDCJ, 
						&candidateMatrix[ index ][ h ]->numPointsDCJ, 
						i, j, FALSE );//--from dcjdist.c
					
					calculateInverseGenome( 
						candidateMatrix[ index ][ h ]->genomeDCJ, 
						candidateMatrix[ index ][ h ]->numPointsDCJ, 
						candidateMatrix[ index ][ h ]->inverseDCJ );//--from dcjdist.c

					verifyPenalization( phyloTreePtr, 
									candidateMatrix[ index ], &h, paramsPtr );
					
				}
			}
		}
		numCandidates[ index ] = h;
		
		//important note: free this in some function
		//for ( i = 0; i < numCandidates; i++) {
		//	freeCandidate( phyloTreePtr, &candidateMatrix[ index ][ i ] );
		//}
		//free( candidateMatrix[ index ] );
	} //end-if
}

static void copyCurrentNodeIntoCandidate( TreePtr phyloTreePtr
				CandidatePtr candPtr, TreeNodePtr nodePtr, int copyInverse ) 
{
	int k;

	candPtr->numPointsDCJ = nodePtr->numPointsDCJ;
	for ( k = 0; k < nodePtr->numPointsDCJ; k++ ) {
		candPtr->genomeDCJ[ k ]->x = nodePtr->genomeDCJ[ k ]->x;
		candPtr->genomeDCJ[ k ]->y = nodePtr->genomeDCJ[ k ]->y;
		candPtr->genomeDCJ[ k ]->type = nodePtr->genomeDCJ[ k ]->type;
	}

	if ( copyInverse == TRUE ) {
		for ( k = 0; k < 2 * phyloTreePtr->numberGenes; k++ ) {
			candPtr->inverseDCJ[ k ] = nodePtr->inverseDCJ[ k ];
		}	
	}
}

/* NOTE: in the first call of this procedure, the second argument 
			must be the right descendant of the starting node */
static void calculateScoreCandidates( TreePtr phyloTreePtr, 
	TreeNodePtr nodePtr, int **scoreMatrix, CandidatePtr **candidateMatrix, 
	int numCandidates, ParametersPtr paramsPtr ) 
{

	//CONTINUE HERE .....




	int i, j, min_left, min_right, temp, index1, index2;

	if ( nodePtr->type == INTERNAL_NODE ) {
		calculateScoreCandidates( phyloTreePtr, nodePtr->leftDescPtr, 
			scoreMatrix, candidateMatrix, numCandidates, paramsPtr );
		calculateScoreCandidates( phyloTreePtr, nodePtr->rightDescPtr,
			scoreMatrix, candidateMatrix, numCandidates, paramsPtr );

		index1 = nodePtr->index - phyloTreePtr->numberLeaves; // re-calculate index
		/* calculate min scores of candidates of current node (nodePtr)  */
		for ( i = 0; i < numCandidates; i++ ) {
			/* calculate min score of candidate "i" regarding left descendant candidates */
			if ( nodePtr->leftDescPtr->type == LEAF_NODE ) {
				min_left = calculateDistance( candidateMatrix[ index1 ][ i ]->genome, 
												nodePtr->leftDescPtr->genome, 
												phyloTreePtr->numberGenes, paramsPtr );	
			}
			else { //nodePtr->leftDescPtr->type == INTERNAL_NODE
				j = 0;
				index2 = nodePtr->leftDescPtr->index - phyloTreePtr->numberLeaves;// re-calculate index
				min_left = scoreMatrix[index2][j] + 
							calculateDistance( candidateMatrix[ index1 ][ i ]->genome, 
												candidateMatrix[ index2 ][ j ]->genome, 
												phyloTreePtr->numberGenes, paramsPtr );

				for ( j = 1; j < numCandidates; j++ ) {
					temp = scoreMatrix[index2][j] + 
							calculateDistance( candidateMatrix[ index1 ][ i ]->genome, 
												candidateMatrix[ index2 ][ j ]->genome, 
												phyloTreePtr->numberGenes, paramsPtr );
					if (temp < min_left) {
						min_left = temp;
					}
				}
			}

			/* calculate min score of candidate "i" regarding right descendant candidates */
			if ( nodePtr->rightDescPtr->type == LEAF_NODE ) {
				min_right = calculateDistance( candidateMatrix[ index1 ][ i ]->genome, 
												nodePtr->rightDescPtr->genome, 
												phyloTreePtr->numberGenes, paramsPtr );
			}
			else { //nodePtr->rightDescPtr->type == INTERNAL_NODE
				j = 0;
				index2 = nodePtr->rightDescPtr->index - phyloTreePtr->numberLeaves;// re-calculate index
				min_right = scoreMatrix[index2][j] + 
							calculateDistance( candidateMatrix[ index1 ][ i ]->genome, 
												candidateMatrix[ index2 ][ j ]->genome, 
												phyloTreePtr->numberGenes, paramsPtr);

				for ( j = 1; j < numCandidates; j++ ) {
					temp = scoreMatrix[index2][j] + 
							calculateDistance( candidateMatrix[ index1 ][ i ]->genome, 
												candidateMatrix[ index2 ][ j ]->genome, 
												phyloTreePtr->numberGenes, paramsPtr );
					if (temp < min_right) {
						min_right = temp;
					}
				}
			}
			scoreMatrix[ index1 ][ i ] = min( min_left, min_right );
		} // end for
	}
}

/* NOTE: in the first call of this procedure, the second argument 
			must be the right descendant of the starting node */
static void selectBestNeighborhood( TreePtr phyloTreePtr, TreeNodePtr nodePtr, int **scoreMatrix, 
	CandidatePtr **candidateMatrix, int numCandidates, ParametersPtr paramsPtr )
{
	int i, k, index, score, newScore, iBest;

	if ( nodePtr->type == INTERNAL_NODE ) {
		index = nodePtr->index - phyloTreePtr->numberLeaves;
		
		i = 0;
		iBest = i;
		score = scoreMatrix[ index ][ i ] + 
				calculateDistance( candidateMatrix[ index ][ i ]->genome, 
									nodePtr->ancestorPtr->genome, 
									phyloTreePtr->numberGenes, paramsPtr );

		for ( i = 1; i < numCandidates; i++ ) {
			newScore = scoreMatrix[ index ][ i ] + 
						calculateDistance( candidateMatrix[ index ][ i ]->genome, 
											nodePtr->ancestorPtr->genome, 
											phyloTreePtr->numberGenes, paramsPtr );
			if ( newScore < score ) {
				score = newScore;
				iBest = i;
			}
		}

		for ( k = 0; k < phyloTreePtr->numberGenes; k++ ) {
			nodePtr->genome[ k ] = candidateMatrix[ index ][ iBest ]->genome[ k ];
		}

		selectBestNeighborhood( phyloTreePtr, nodePtr->leftDescPtr, 
			scoreMatrix, candidateMatrix, numCandidates, paramsPtr );
		selectBestNeighborhood( phyloTreePtr, nodePtr->rightDescPtr, 
			scoreMatrix, candidateMatrix, numCandidates, paramsPtr );
	}

	//TODO...free memory allocated in function "generateCandidates"
}


/* [OPTIMIZER for DCJ distance] 
	algorithm for labeling and optimizing the score of a tree by 
	generating candidates for each internal node and greedily improve the tree */
static int labelOptimizeTree_GreedyCandidatesDCJ( 
	TreePtr phyloTreePtr, int initialize,  ParametersPtr paramsPtr ) 
{
	int improve, score, newScore;
	Tree tempTree;

	/* allocate memory for temporary tree */
	tempTree.numberLeaves = phyloTreePtr->numberLeaves;
	tempTree.numberGenes = phyloTreePtr->numberGenes;
	allocateMemoryForNodes( &tempTree, paramsPtr );//--method from tree.c

	if ( initialize == TRUE ) {
		if ( paramsPtr->problem == LARGE_PHYLOGENY ) {
			initializeTreeWithDescendants( phyloTreePtr, 
			phyloTreePtr->startingNodePtr->rightDescPtr, paramsPtr, GO_LEFT );
		}
		else { // SMALL_PHYLOGENY
			initializeTreeWithDescendantsRandomly( 
				phyloTreePtr, phyloTreePtr->startingNodePtr->rightDescPtr );
		}
	}

	improve = TRUE;
	score = scoreTree( phyloTreePtr, phyloTreePtr->startingNodePtr->rightDescPtr, paramsPtr );
	while ( improve == TRUE ) {
		copyTreeInto( &tempTree, phyloTreePtr, FALSE, paramsPtr );/* make a copy of current tree */
		improveTreebyCandidatesDCJ( phyloTreePtr, 
			phyloTreePtr->startingNodePtr->rightDescPtr, paramsPtr );
		newScore = scoreTree( phyloTreePtr, 
			phyloTreePtr->startingNodePtr->rightDescPtr, paramsPtr );

		if ( newScore < score ) {
			score = newScore;
		}
		else {
			improve = FALSE;
			copyTreeInto( phyloTreePtr, &tempTree, FALSE, paramsPtr ); /* recover last tree */
		}
	}

	freeTree( &tempTree, paramsPtr );//--method from tree.c
	return score;
}

/* Initialization proposed by Herencsar for PIVO2 see paper "An Improved
	Algorithm for Ancestral Gene Order Reconstruction" */
/* NOTE: the input (second parameter) must be the RIGHT descendant of the
	startingNode of the tree */
static void initializeTreeWithDescendants( TreePtr phyloTreePtr, 
		TreeNodePtr nodePtr, ParametersPtr paramsPtr, int orientation )
{
	int i;

	if ( nodePtr->type == INTERNAL_NODE ) {
		if ( orientation == GO_LEFT ) {
			initializeTreeWithDescendants( phyloTreePtr, nodePtr->leftDescPtr, paramsPtr, GO_RIGHT );
			initializeTreeWithDescendants( phyloTreePtr, nodePtr->rightDescPtr, paramsPtr, GO_LEFT );
		}
		else { // orientation == GO_RIGHT
			initializeTreeWithDescendants( phyloTreePtr, nodePtr->leftDescPtr, paramsPtr, GO_LEFT );
			initializeTreeWithDescendants( phyloTreePtr, nodePtr->rightDescPtr, paramsPtr, GO_RIGHT );
		}
		//p = irand(2); // 0 <= p <= 1

		if ( orientation == GO_LEFT ) {
			/* copy left descendant */
			nodePtr->numPointsDCJ = nodePtr->leftDescPtr->numPointsDCJ;
			for ( i = 0; i < nodePtr->numPointsDCJ; i++ ) {
				nodePtr->genomeDCJ[i]->x = nodePtr->leftDescPtr->genomeDCJ[i]->x;
				nodePtr->genomeDCJ[i]->y = nodePtr->leftDescPtr->genomeDCJ[i]->y;
				nodePtr->genomeDCJ[i]->type = nodePtr->leftDescPtr->genomeDCJ[i]->type;
			}

			/* copy inverse */
			for ( i = 0; i < 2 * phyloTreePtr->numberGenes; i++ ) {
				nodePtr->inverseDCJ[ i ] = nodePtr->leftDescPtr->inverseDCJ[ i ];
			}
		}
		else { // orientation == GO_RIGHT
			/* copy right descendant */
			nodePtr->numPointsDCJ = nodePtr->rightDescPtr->numPointsDCJ;
			for ( i = 0; i < nodePtr->numPointsDCJ; i++ ) {
				nodePtr->genomeDCJ[i]->x = nodePtr->rightDescPtr->genomeDCJ[i]->x;
				nodePtr->genomeDCJ[i]->y = nodePtr->rightDescPtr->genomeDCJ[i]->y;
				nodePtr->genomeDCJ[i]->type = nodePtr->rightDescPtr->genomeDCJ[i]->type;
			}

			/* copy inverse */
			for ( i = 0; i < 2 * phyloTreePtr->numberGenes; i++ ) {
				nodePtr->inverseDCJ[ i ] = nodePtr->rightDescPtr->inverseDCJ[ i ];
			}
		}

	}
}

static void initializeTreeWithDescendantsRandomly( 
									TreePtr phyloTreePtr, TreeNodePtr nodePtr )
{
	int i, p, orientation;

	if ( nodePtr->type == INTERNAL_NODE ) {
		initializeTreeWithDescendantsRandomly( phyloTreePtr, nodePtr->leftDescPtr );
		initializeTreeWithDescendantsRandomly( phyloTreePtr, nodePtr->rightDescPtr );
		
		p = irand(2); // 0 <= p <= 1
		orientation = p == 0 ? GO_LEFT : GO_RIGHT;

		if ( orientation == GO_LEFT ) {
			/* copy left descendant */
			nodePtr->numPointsDCJ = nodePtr->leftDescPtr->numPointsDCJ;
			for ( i = 0; i < nodePtr->numPointsDCJ; i++ ) {
				nodePtr->genomeDCJ[ i ]->x = 
									nodePtr->leftDescPtr->genomeDCJ[ i ]->x;
				nodePtr->genomeDCJ[ i ]->y = 
									nodePtr->leftDescPtr->genomeDCJ[ i ]->y;
				nodePtr->genomeDCJ[ i ]->type = 
									nodePtr->leftDescPtr->genomeDCJ[ i ]->type;
			}

			/* copy inverse */
			for ( i = 0; i < 2 * phyloTreePtr->numberGenes; i++ ) {
				nodePtr->inverseDCJ[ i ] = 
									nodePtr->leftDescPtr->inverseDCJ[ i ];
			}
		}
		else { // orientation == GO_RIGHT
			/* copy right descendant */
			nodePtr->numPointsDCJ = nodePtr->rightDescPtr->numPointsDCJ;
			for ( i = 0; i < nodePtr->numPointsDCJ; i++ ) {
				nodePtr->genomeDCJ[ i ]->x = 
								nodePtr->rightDescPtr->genomeDCJ[ i ]->x;
				nodePtr->genomeDCJ[ i ]->y = 
								nodePtr->rightDescPtr->genomeDCJ[ i ]->y;
				nodePtr->genomeDCJ[ i ]->type = 
								nodePtr->rightDescPtr->genomeDCJ[ i ]->type;
			}

			/* copy inverse */
			for ( i = 0; i < 2 * phyloTreePtr->numberGenes; i++ ) {
				nodePtr->inverseDCJ[ i ] = 
								nodePtr->rightDescPtr->inverseDCJ[ i ];
			}
		}

	}//end-if
}

/* NOTE: in the first call of this procedure, the second argument 
			must be the right descendant of the starting node */
static void improveTreebyCandidatesDCJ( TreePtr phyloTreePtr, 
	TreeNodePtr nodePtr, ParametersPtr paramsPtr ) 
{
	int a, b, c, newA, newB, newC, i, j, k, h, numA, numT;
	int score, newScore, maxElements, numCandidates, update;
	CandidatePtr *candidatePtrArray;
	int x, y, xpos, ypos;

	if ( nodePtr->type == INTERNAL_NODE ) {
		improveTreebyCandidatesDCJ( phyloTreePtr, nodePtr->leftDescPtr, paramsPtr );
		improveTreebyCandidatesDCJ( phyloTreePtr, nodePtr->rightDescPtr, paramsPtr );

		score = nodePtr->edgeWeight + 
				nodePtr->leftDescPtr->edgeWeight + 
				nodePtr->rightDescPtr->edgeWeight;

		/* calculate max number of candidates for the current node */
		numA = 0; //number of adjacencies
		numT = 0; //number of telomeres
		for ( i = 0; i < nodePtr->numPointsDCJ; i++ ) {
			if ( nodePtr->genomeDCJ[ i ]->type == ADJACENCY )
				numA++;
			else //nodePtr->genomeDCJ[ i ]->type == TELOMERE
				numT++;
		}
		maxElements = 	
			numA + 							// for inv. of case (3)
			2 * numA * ( numA - 1 ) / 2 + 	// 2 * C(numA, 2), for case (1)
			2 * numA * numT +				// for case (2)
			numT / 2;						// for case (3)

		/* allocate memory for array of pointers to candidate */
		candidatePtrArray = malloc( maxElements * sizeof( CandidatePtr ) );
		if ( candidatePtrArray == NULL ) nomemMessage("candidates");
		
		/* generate candidates for inverse of case (3) */
		h = 0;
		for ( i = 0; i < nodePtr->numPointsDCJ; i++ ) {
			if ( nodePtr->genomeDCJ[ i ]->type == ADJACENCY ) {
				/* calculate position for checking inverses genomes */
				x = nodePtr->genomeDCJ[ i ]->x;
				y = nodePtr->genomeDCJ[ i ]->y;
				xpos = x > 0 ? 2 * x - 1 : 2 * abs( x ) - 2;
				ypos = y > 0 ? 2 * y - 1 : 2 * abs( y ) - 2;

				/* if adjacency (i->x, i->y) is also a adjacency in ancestor, 
					left desc., or right desc. do nothing */
				if ( nodePtr->ancestorPtr->inverseDCJ[ xpos ] == nodePtr->ancestorPtr->inverseDCJ[ ypos ] && 
					nodePtr->leftDescPtr->inverseDCJ[ xpos ] == nodePtr->leftDescPtr->inverseDCJ[ ypos ] && 
					nodePtr->rightDescPtr->inverseDCJ[ xpos ] == nodePtr->rightDescPtr->inverseDCJ[ ypos ] ) {
					continue;
				}

				allocateMemForCandidate( phyloTreePtr, &candidatePtrArray[ h ] );

				/* copy genomeDCJ */
				candidatePtrArray[ h ]->numPointsDCJ = nodePtr->numPointsDCJ;
				for ( k = 0; k < nodePtr->numPointsDCJ; k++ ) {
					candidatePtrArray[ h ]->genomeDCJ[ k ]->x = nodePtr->genomeDCJ[ k ]->x;
					candidatePtrArray[ h ]->genomeDCJ[ k ]->y = nodePtr->genomeDCJ[ k ]->y;
					candidatePtrArray[ h ]->genomeDCJ[ k ]->type = nodePtr->genomeDCJ[ k ]->type;
				}

				/* apply DCJ on copy and generate candidate */
				applyDCJ( candidatePtrArray[ h ]->genomeDCJ, 
					&candidatePtrArray[ h ]->numPointsDCJ, i, i, TRUE );//--from dcjdist.c
			
				calculateInverseGenome( candidatePtrArray[ h ]->genomeDCJ, 
										candidatePtrArray[ h ]->numPointsDCJ, 
										candidatePtrArray[ h ]->inverseDCJ );//--from dcjdist.c

				verifyPenalization( phyloTreePtr, candidatePtrArray, &h, paramsPtr );
			}
		}

		/* generate candidates by choosing two different indices: case (1), (2), and (3) */ 
		for ( i = 0; i < nodePtr->numPointsDCJ - 1; i++ ) {
			for ( j = i + 1; j < nodePtr->numPointsDCJ; j++ ) {
				/* calculate position for checking inverses genomes */
				x = nodePtr->genomeDCJ[ i ]->x;
				y = nodePtr->genomeDCJ[ i ]->y;
				xpos = x > 0 ? 2 * x - 1 : 2 * abs( x ) - 2;
				ypos = y > 0 ? 2 * y - 1 : 2 * abs( y ) - 2;

				/* if adjacency (i->x, i->y) is also a adjacency in ancestor, 
					left desc., or right desc. do nothing */
				if ( nodePtr->ancestorPtr->inverseDCJ[ xpos ] == nodePtr->ancestorPtr->inverseDCJ[ ypos ] && 
					nodePtr->leftDescPtr->inverseDCJ[ xpos ] == nodePtr->leftDescPtr->inverseDCJ[ ypos ] && 
					nodePtr->rightDescPtr->inverseDCJ[ xpos ] == nodePtr->rightDescPtr->inverseDCJ[ ypos ] ) {
					continue;
				}

				/* calculate position for checking inverses genomes */
				x = nodePtr->genomeDCJ[ j ]->x;
				y = nodePtr->genomeDCJ[ j ]->y;
				xpos = x > 0 ? 2 * x - 1 : 2 * abs( x ) - 2;
				ypos = y > 0 ? 2 * y - 1 : 2 * abs( y ) - 2;

				/* if adjacency (j->x,j->y) is also a adjacency in ancestor, 
					left desc., or right desc. do nothing */
				if ( nodePtr->ancestorPtr->inverseDCJ[ xpos ] == nodePtr->ancestorPtr->inverseDCJ[ ypos ] && 
					nodePtr->leftDescPtr->inverseDCJ[ xpos ] == nodePtr->leftDescPtr->inverseDCJ[ ypos ] && 
					nodePtr->rightDescPtr->inverseDCJ[ xpos ] == nodePtr->rightDescPtr->inverseDCJ[ ypos ] ) {
					continue;
				}

				allocateMemForCandidate( phyloTreePtr, &candidatePtrArray[ h ] );

				/* copy genomeDCJ */
				candidatePtrArray[ h ]->numPointsDCJ = nodePtr->numPointsDCJ;
				for ( k = 0; k < nodePtr->numPointsDCJ; k++ ) {
					candidatePtrArray[ h ]->genomeDCJ[ k ]->x = nodePtr->genomeDCJ[ k ]->x;
					candidatePtrArray[ h ]->genomeDCJ[ k ]->y = nodePtr->genomeDCJ[ k ]->y;
					candidatePtrArray[ h ]->genomeDCJ[ k ]->type = nodePtr->genomeDCJ[ k ]->type;
				}

				/* apply DCJ on copy and generate candidate */
				applyDCJ( candidatePtrArray[ h ]->genomeDCJ, &candidatePtrArray[ h ]->numPointsDCJ, i, j, TRUE );//--from dcjdist.c
				
				calculateInverseGenome( candidatePtrArray[ h ]->genomeDCJ, 
										candidatePtrArray[ h ]->numPointsDCJ, 
										candidatePtrArray[ h ]->inverseDCJ );//--from dcjdist.c

				verifyPenalization( phyloTreePtr, candidatePtrArray, &h, paramsPtr );

				/* do the second FORM in case one of the point is a adjacency */
				if ( nodePtr->genomeDCJ[ i ]->type == ADJACENCY || 
						nodePtr->genomeDCJ[ j ]->type == ADJACENCY ) {

					allocateMemForCandidate( phyloTreePtr, &candidatePtrArray[ h ] );

					/* copy genomeDCJ */
					candidatePtrArray[ h ]->numPointsDCJ = nodePtr->numPointsDCJ;
					for ( k = 0; k < nodePtr->numPointsDCJ; k++ ) {
						candidatePtrArray[ h ]->genomeDCJ[ k ]->x = nodePtr->genomeDCJ[ k ]->x;
						candidatePtrArray[ h ]->genomeDCJ[ k ]->y = nodePtr->genomeDCJ[ k ]->y;
						candidatePtrArray[ h ]->genomeDCJ[ k ]->type = nodePtr->genomeDCJ[ k ]->type;
					}
					/* apply DCJ on copy and generate candidate */
					applyDCJ( candidatePtrArray[ h ]->genomeDCJ, &candidatePtrArray[ h ]->numPointsDCJ, i, j, FALSE );//--from dcjdist.c
					
					calculateInverseGenome( candidatePtrArray[ h ]->genomeDCJ, 
										candidatePtrArray[ h ]->numPointsDCJ, 
										candidatePtrArray[ h ]->inverseDCJ );//--from dcjdist.c

					verifyPenalization( phyloTreePtr, candidatePtrArray, &h, paramsPtr );
					
				}
			}
		}
		numCandidates = h;

		/* search for best candidate */
		update = FALSE;
		for (i = 0; i < numCandidates; i++) {
			newA = DCJdistance( candidatePtrArray[ i ]->genomeDCJ, nodePtr->ancestorPtr->genomeDCJ,
								candidatePtrArray[ i ]->inverseDCJ, nodePtr->ancestorPtr->inverseDCJ,  
								candidatePtrArray[ i ]->numPointsDCJ, nodePtr->ancestorPtr->numPointsDCJ,
								phyloTreePtr->numberGenes );
			newB = DCJdistance( candidatePtrArray[ i ]->genomeDCJ, nodePtr->leftDescPtr->genomeDCJ,
								candidatePtrArray[ i ]->inverseDCJ, nodePtr->leftDescPtr->inverseDCJ, 
								candidatePtrArray[ i ]->numPointsDCJ, nodePtr->leftDescPtr->numPointsDCJ,
								phyloTreePtr->numberGenes );
			newC = DCJdistance( candidatePtrArray[ i ]->genomeDCJ, nodePtr->rightDescPtr->genomeDCJ,
								candidatePtrArray[ i ]->inverseDCJ, nodePtr->rightDescPtr->inverseDCJ, 
								candidatePtrArray[ i ]->numPointsDCJ, nodePtr->rightDescPtr->numPointsDCJ,
								phyloTreePtr->numberGenes );
			newScore = newA + newB + newC;

			if ( newScore < score) {
				//printf("-->Improved by Kovac MOD. Optimization.\n");
				score = newScore;
				a = newA; b = newB; c = newC;
				update = TRUE;
				h = i;
			}
		}

		/* update current node with the best candidate found */
		if ( update == TRUE ) {
			nodePtr->numPointsDCJ = candidatePtrArray[ h ]->numPointsDCJ;
			for ( k = 0; k < nodePtr->numPointsDCJ; k++ ) {
				nodePtr->genomeDCJ[ k ]->x = candidatePtrArray[ h ]->genomeDCJ[ k ]->x;
				nodePtr->genomeDCJ[ k ]->y = candidatePtrArray[ h ]->genomeDCJ[ k ]->y;
				nodePtr->genomeDCJ[ k ]->type = candidatePtrArray[ h ]->genomeDCJ[ k ]->type;
			}
			nodePtr->edgeWeight = a; 
			nodePtr->leftDescPtr->edgeWeight = b; 
			nodePtr->rightDescPtr->edgeWeight = c;

			/* copy inverse of candidate */
			for ( j = 0; j < 2 * phyloTreePtr->numberGenes; j++ ) {
				nodePtr->inverseDCJ[ j ] = candidatePtrArray[ h ]->inverseDCJ[ j ];
			}
		}

		/* free memory */
		for ( i = 0; i < numCandidates; i++) {
			freeCandidate( phyloTreePtr, &candidatePtrArray[ i ] );
		}
		free( candidatePtrArray );
	} // end-if
}

static void allocateMemForCandidate( TreePtr phyloTreePtr, CandidatePtr *candPtr )
{ 
	int j;

	( *candPtr ) = malloc( sizeof( Candidate ) );
	if ( ( *candPtr ) == NULL ) nomemMessage( "(*candPtr)" );

	( *candPtr )->genomeDCJ = 
		malloc( 2 * phyloTreePtr->numberGenes * sizeof( PointDCJPtr ) );
	if ( ( *candPtr )->genomeDCJ == NULL ) 
		nomemMessage( "(*candPtr)->genomeDCJ" ); 

	for ( j = 0; j < 2 * phyloTreePtr->numberGenes; j++ ) {
		( *candPtr )->genomeDCJ[ j ] = malloc( sizeof( PointDCJ ) );
		if ( ( *candPtr )->genomeDCJ[ j ] == NULL )
			nomemMessage( "(*candPtr)->genomeDCJ[j]" );
	}

	( *candPtr )->inverseDCJ = malloc( 2 * phyloTreePtr->numberGenes * sizeof( int ) );
	if ( ( *candPtr )->inverseDCJ == NULL ) nomemMessage( "(*candPtr)->inverseDCJ" );

	( *candPtr )->numPointsDCJ = 0;
}

static void freeCandidate( TreePtr phyloTreePtr, CandidatePtr *candPtr )
{
	int j;

	free( ( *candPtr )->inverseDCJ );
	for ( j = 0; j < 2 * phyloTreePtr->numberGenes; j++ ) {
		free( ( *candPtr )->genomeDCJ[ j ] );
	}
	free( ( *candPtr )->genomeDCJ );
	free( ( *candPtr ) );
}

static void verifyPenalization( TreePtr phyloTreePtr, 
	CandidatePtr *candidatePtrArray, int *h, ParametersPtr paramsPtr ) 
{
	/* if option of penalize chromosomes was chosen and the candidate 
		was penalized then free candidate. Otherwise, increment the 
		counter of candidates */
	/* Note: function candidatePenalized is from dcjdist.c */
	if ( paramsPtr->penaltyType != -1 && 
		candidatePenalized( candidatePtrArray[ ( *h ) ], paramsPtr ) == TRUE ) 
	{
		freeCandidate( phyloTreePtr, &candidatePtrArray[ ( *h ) ] );							
	}
	else {
		( *h )++;
	}
}

/* initialize tree using a depth first search (DFS) approach.
 * the input (second parameter) must be the startingNode of the tree */
static void initializeTreeUsingDFS( TreePtr phyloTreePtr, 
		TreeNodePtr nodePtr, ParametersPtr paramsPtr )
{
	/* NOTE:
	* If a node is an extremity (==TRUE) means that it is already labeled
	* otherwise the node is not labeled */
	
	TreeNodePtr node1Ptr;
	TreeNodePtr node2Ptr;
	TreeNodePtr node3Ptr;
	struct adj_struct *adjacencyList; /* from structs.h */

	if ( nodePtr != NULL ) {
		/* if node is not an extremity, 
		* then setup a TSP instance and label the node */
		if ( nodePtr->extremity == FALSE ) {
			/* allocate memory*/
			adjacencyList = 
			malloc ( ( 2*phyloTreePtr->numberGenes + 1 ) * sizeof( struct adj_struct ) );
    		if ( adjacencyList == NULL ) nomemMessage( "adjacencyList" ); 

			/* search for the 3 nearest neighbors to nodePtr, which are labeled,
			* create a TSP instance, and call a TSP solver to label nodePtr*/
			node1Ptr = findNearestNeighbor( phyloTreePtr, nodePtr->ancestorPtr, -1, -1 );
    		node2Ptr = findNearestNeighbor( phyloTreePtr, nodePtr->leftDescPtr, node1Ptr->id, -1 );
    		node3Ptr = findNearestNeighbor( phyloTreePtr, nodePtr->rightDescPtr, node1Ptr->id, node2Ptr->id );

			createTspInstance( phyloTreePtr->numberGenes, node1Ptr->genome, 
				node2Ptr->genome, node3Ptr->genome, adjacencyList );//--method from median_solvers.c
			callSolverAndLabelNode( phyloTreePtr, nodePtr, node1Ptr, node2Ptr, 
				node3Ptr, adjacencyList, paramsPtr->solver, paramsPtr->circular );
			
			nodePtr->extremity = TRUE; /*node is now labeled*/

			/* free memory */
			free( adjacencyList );
		}
		/* visit the left and right descendant */
		initializeTreeUsingDFS( phyloTreePtr, nodePtr->leftDescPtr, paramsPtr );
		initializeTreeUsingDFS( phyloTreePtr, nodePtr->rightDescPtr, paramsPtr );	
	}
}

/* find nearest neighbor, which is labeled, using BFS */
static TreeNodePtr findNearestNeighbor( TreePtr phyloTreePtr, 
		TreeNodePtr nodePtr, int idNeighbor1, int idNeighbor2 )
{
	Queue myQueue;
	TreeNodePtr nodeAnsPtr;
	TreeNodePtr nodeAuxPtr;
	int *visited;

	/* allocate memory */
	visited = calloc( phyloTreePtr->numberNodes + 1, sizeof(int) );
	if ( visited == NULL ) nomemMessage( "visited" );

 	/* if node is labeled, and is different of the other neighbors */
	if ( nodePtr->extremity == TRUE && nodePtr->id != idNeighbor1 && 
										nodePtr->id != idNeighbor2 ) 
	{
		nodeAnsPtr = nodePtr;
	}
	else {
		visited[ nodePtr->id ] = TRUE;
		enqueue( &myQueue, nodePtr );
		
		while ( isEmpty( &myQueue ) == FALSE ) {
			nodeAuxPtr = dequeue( &myQueue );
			/* if nodeAux is labeled, and is different of the other neighbors */
			if ( nodeAuxPtr->extremity == TRUE && nodeAuxPtr->id != idNeighbor1 && 
												nodeAuxPtr->id != idNeighbor2 ) 
			{
				nodeAnsPtr = nodeAuxPtr;
				break;
			}
			else { /* otherwise visit adjacent nodes */
				if ( nodeAuxPtr->ancestorPtr != NULL && visited[ nodeAuxPtr->ancestorPtr->id ] == FALSE ) {
					visited[ nodeAuxPtr->ancestorPtr->id ] = TRUE;
					enqueue( &myQueue, nodeAuxPtr->ancestorPtr );
				}

				if ( nodeAuxPtr->leftDescPtr != NULL && visited[ nodeAuxPtr->leftDescPtr->id ] == FALSE ){
					visited[ nodeAuxPtr->leftDescPtr->id ] = TRUE;
					enqueue( &myQueue, nodeAuxPtr->leftDescPtr );
					
				}

				if ( nodeAuxPtr->rightDescPtr != NULL && visited[ nodeAuxPtr->rightDescPtr->id ] == FALSE ){					
					visited[ nodeAuxPtr->rightDescPtr->id ] = TRUE;
					enqueue( &myQueue, nodeAuxPtr->rightDescPtr );	
				}
			}
		}//end while
	}

	/* free memory */
	free( visited );

	return nodeAnsPtr;
}

/* call a median solver and label nodePtr */
/*NOTE: the median genome will returned in nodePtr*/
static void callSolverAndLabelNode( TreePtr phyloTreePtr, TreeNodePtr nodePtr, 
	TreeNodePtr node1Ptr, TreeNodePtr node2Ptr, TreeNodePtr node3Ptr,
	struct adj_struct *adjacencyList, enum medianSolvers solver, int circular )
{	
	int *tspCycle;

	/* allocate memory */
	tspCycle = malloc ( ( 2*phyloTreePtr->numberGenes + 1 ) * sizeof(int) );
    if ( tspCycle == NULL ) nomemMessage( "tspCycle" ); 

    /* call a median solver */
	if ( solver == COALESTSP ) {
		callCoalestsp( phyloTreePtr->numberGenes, node1Ptr->genome, 
	   				node2Ptr->genome, node3Ptr->genome, adjacencyList, 
	   				nodePtr->genome, tspCycle, circular );//--method from median_solvers.c
    }
    else if ( solver == BBTSP ) {
        callBbtsp( phyloTreePtr->numberGenes, node1Ptr->genome, 
        			node2Ptr->genome, node3Ptr->genome, adjacencyList, 
        			nodePtr->genome, tspCycle, circular );//--method from median_solvers.c
    }
    else if ( solver == CAPRARA_INV_MEDIAN ) {
        callCapraraInversionMedian( phyloTreePtr->numberGenes, node1Ptr->genome, 
        	node2Ptr->genome, node3Ptr->genome, nodePtr->genome, circular );
    }
    else{
        fprintf(stderr, "stderr: not a valid median solver.\n");
        exit(EXIT_FAILURE); 
    }

    /* free memory */
    free(tspCycle);
}

/* Score a tree using some metric. 
* The input (second parameter: nodePtr), 
* must be the RIGHT descendant of the STARTING node of the tree. */
int scoreTree( TreePtr phyloTreePtr, TreeNodePtr nodePtr, ParametersPtr paramsPtr )
{
	/* NOTE:
	* All nodes, except the starting node, have an EDGE with its ancestor.
	* Indeed, the STARTING node is the only node without ancestor. 
	* So, we have N-1 edges with each node representing an EDGE, where N is the number of nodes */

	int distance;

	if ( paramsPtr->distanceType == DCJ_DIST ) {
		if ( nodePtr->type == LEAF_NODE ) {
			nodePtr->edgeWeight = DCJdistance( 
								nodePtr->genomeDCJ, 
								nodePtr->ancestorPtr->genomeDCJ,
								nodePtr->inverseDCJ, 
								nodePtr->ancestorPtr->inverseDCJ,
								nodePtr->numPointsDCJ, 
								nodePtr->ancestorPtr->numPointsDCJ, 
								phyloTreePtr->numberGenes );//--from dcjdist.c
			distance = nodePtr->edgeWeight;
		}
		else { /* type is INTERNAL NODE */
			nodePtr->edgeWeight = DCJdistance( 
								nodePtr->genomeDCJ, 
								nodePtr->ancestorPtr->genomeDCJ,
								nodePtr->inverseDCJ, 
								nodePtr->ancestorPtr->inverseDCJ, 
								nodePtr->numPointsDCJ, 
								nodePtr->ancestorPtr->numPointsDCJ, 
								phyloTreePtr->numberGenes );//--from dcjdist.c 
			distance = nodePtr->edgeWeight + 
				scoreTree( phyloTreePtr, nodePtr->leftDescPtr, paramsPtr ) + 
				scoreTree( phyloTreePtr, nodePtr->rightDescPtr, paramsPtr );
		}	
	}
	else { // INVERSION_DIST or BREAKPOINT_DIST
		if ( nodePtr->type == LEAF_NODE ) {
			nodePtr->edgeWeight = calculateDistance( 
									nodePtr->genome, 
									nodePtr->ancestorPtr->genome, 
									phyloTreePtr->numberGenes, paramsPtr );
			distance = nodePtr->edgeWeight;
		}
		else { /* type is INTERNAL NODE */
			nodePtr->edgeWeight = calculateDistance( 
									nodePtr->genome, 
									nodePtr->ancestorPtr->genome, 
									phyloTreePtr->numberGenes, paramsPtr ); 
			distance = nodePtr->edgeWeight + 
				scoreTree( phyloTreePtr, nodePtr->leftDescPtr, paramsPtr ) + 
				scoreTree( phyloTreePtr, nodePtr->rightDescPtr, paramsPtr );
		}
	}
	return distance;
}

int calculateDistance( int *genome1, int *genome2, int numberGenes, ParametersPtr paramsPtr )
{
    int i;
    int distance;
    struct genome_struct *genome_list;

    /* allocate memory */
    genome_list = malloc ( 2 * sizeof ( struct genome_struct ) );
    if ( genome_list == NULL ) nomemMessage( "genome_list" ); 

    genome_list[ 0 ].genome_num = 1; //id
    genome_list[ 0 ].genes = malloc ( numberGenes * sizeof ( int ) );
    if ( genome_list[ 0 ].genes == NULL ) nomemMessage( "genome_list[ 0 ].genes" ); 

    genome_list[ 1 ].genome_num = 2; //id
    genome_list[ 1 ].genes = malloc ( numberGenes * sizeof ( int ) );
    if ( genome_list[ 1 ].genes == NULL ) nomemMessage( "genome_list[ 1 ].genes" ); 

    /* Copy genome1 and genome2 into genome_list */
    for ( i = 0; i < numberGenes; i++ ) {
        genome_list[ 0 ].genes[ i ] = genome1[ i ];
        genome_list[ 1 ].genes[ i ] = genome2[ i ];
    }

    /* calculate distance */
    distance = 0;
    if ( paramsPtr->distanceType == INVERSION_DIST ) {
        if ( paramsPtr->circular ) {
            distance = invdist_circular_nomem( 
            	genome_list, genome_list + 1, numberGenes );//from invdist.c (GRAPPA)    
        }
        else {
            distance = invdist_noncircular_nomem( 
            	genome_list, genome_list + 1, 0, numberGenes );//from invdist.c (GRAPPA)  
        } 
    }
    else {// paramsPtr->distanceType == BREAKPOINT_DIST 
        int *geneAdj = malloc ( ( numberGenes + 1 ) * 2 * sizeof ( int ) );
        if ( geneAdj ==  NULL ) nomemMessage( "geneAdj" ); 

        distance = hamming_distance( genome_list, genome_list + 1, numberGenes,
                                      paramsPtr->circular, geneAdj );//from binencode.h (GRAPPA)
        free( geneAdj );
    }

    /* free memory */
    free( genome_list[ 0 ].genes );
    free( genome_list[ 1 ].genes );
    free( genome_list );

    return distance;
}





