#include <stdio.h>
#include <stdlib.h>

#include "auxiliary.h"
#include "my_structs.h"

static void calculateInverseGenome( PointDCJPtr *genomeDCJ, int numPointsDCJ, int *inverseGenome ); 

int DCJdistance( PointDCJPtr *genome1DCJ, PointDCJPtr *genome2DCJ, 
			int numPoints1DCJ, int numPoints2DCJ, int numberGenes ) 
{
	int i, k, pos, startType, endType, pathCounter, searchOnGenome1;
	int g1coord, g2coord; /* genome 1 coordinate is x or y */
	int distance, cycles, oddPaths;
	int *inverseGenome1, *inverseGenome2;
	int *visited;
	PointDCJPtr startPoint, endPoint;

	/* allocate memory */
	inverseGenome1 = malloc( 2 * numberGenes * sizeof( int ) );
	if ( inverseGenome1 == NULL ) { nomemMessage( "inverseGenome1" ); }
	inverseGenome2 = malloc( 2 * numberGenes * sizeof( int ) );
	if ( inverseGenome2 == NULL ) { nomemMessage( "inverseGenome2" ); }
	visited = malloc( numPoints1DCJ * sizeof( int ) );
	if ( visited == NULL ) { nomemMessage( "visited" ); }

	calculateInverseGenome( genome1DCJ, numPoints1DCJ, inverseGenome1 ); 
	calculateInverseGenome( genome2DCJ, numPoints2DCJ, inverseGenome2 );

	for ( i = 0; i < numPoints1DCJ; i++ ) {
		visited[ i ] = FALSE;
	}

	/* traverse the adjacency graph to find cycles and paths */
	distance = -1;
	cycles = 0;
	oddPaths = 0;
	for ( i = 0; i < numPoints1DCJ; i++ ) {
		if ( visited[ i ] == FALSE ) {
			startPoint 		= genome1DCJ[ i ];
			startType 		= genome1DCJ[ i ]->type; 
			g1coord 		= genome1DCJ[ i ]->x;
			searchOnGenome1 = FALSE; /* if FALSE search on genome 2*/
			pathCounter 	= 0;
			visited[ i ] 	= TRUE;

			while ( TRUE ) {
				/* go one step forward */
				if ( searchOnGenome1 == TRUE ) {
					/* calculate index based on g2coord */
					pos = g2coord > 0 ? 2 * g2coord - 1 : 2 * abs( g2coord ) - 2;
					/* find the point in genome1 */
					k = inverseGenome1[ pos ];
					if ( g2coord == genome1DCJ[ k ]->x ) {
						g1coord = genome1DCJ[ k ]->y;
					}
					else { // g2coord == node1Ptr->genomeDCJ[ k ]->y;
						g1coord = genome1DCJ[ k ]->x;
					}
					endPoint 		= genome1DCJ[ k ];
					endType 		= genome1DCJ[ k ]->type;
					searchOnGenome1 = FALSE; // search on genome2
					visited[ k ] 	= TRUE;
					pathCounter++;
				}
				else { //search on genome2
					/* calculate index based on g2coord */
					pos = g1coord > 0 ? 2 * g1coord - 1 : 2 * abs( g1coord ) - 2;
					/* find the point in genome1 */
					k = inverseGenome2[ pos ];
					if ( g1coord == genome2DCJ[ k ]->x ) {
						g2coord = genome2DCJ[ k ]->y;
					}
					else { // g1coord == node2Ptr->genomeDCJ[ k ]->y;
						g2coord = genome2DCJ[ k ]->x;
					}
					endPoint 		= genome2DCJ[ k ];
					endType 		= genome2DCJ[ k ]->type;
					searchOnGenome1 = TRUE;
					pathCounter++;
				}

				/* count the cycles and odd paths if found */
				if ( startType == TELOMERE && endType == TELOMERE ) {
					if ( pathCounter % 2 != 0 ) {
						oddPaths++;
					}
					break;	
				}

				if ( startType == ADJACENCY && 
						endType == ADJACENCY && startPoint == endPoint ) {
					cycles++;
					break;	
				}

				if ( startType == ADJACENCY && endType == TELOMERE ) {
					startType 		= TELOMERE;
					g1coord 		= genome1DCJ[ i ]->y;
					searchOnGenome1 = FALSE; // search on genome2 
				}
			}//end-while
		}//end-if
	}

	distance = numberGenes - ( cycles + oddPaths / 2 );

	/* free memory */
	free( inverseGenome1 );
	free( inverseGenome2 );
	free( visited );

	return distance;
}

/* NOTE: the inverseGenome has a lenght of 2 * numberGenes, each element
of this array represent the position the extreme of a gene (-a or +a) 
in the array genomeDCJ. Suposse we have 3 genes 1, 2, and 3, the
array inverseGenome represents the following extremes:
-1, +1, -2, +2, -3, +3
So in order to obtaing the position of -3, we have to do the following
calculation: 2 * (-3) - 2 = 4. 
In general, if the gene is (+) we do: 2 * gene - 1
otherwise, if the gene is (-) we do: 2 * abs(gene) - 2
*/
static void calculateInverseGenome( PointDCJPtr *genomeDCJ, int numPointsDCJ, int *inverseGenome ) 
{
	int i, k;

	for ( i = 0; i < numPointsDCJ; i++ ) {
		/*  */
		if ( genomeDCJ[ i ]->x > 0 )
			k = 2 * genomeDCJ[ i ]->x - 1;
		else // genomeDCJ[ i ]->x < 0
			k = 2 * abs( genomeDCJ[ i ]->x ) - 2;

		inverseGenome[ k ] = i;
		
		if ( genomeDCJ[ i ]->type == ADJACENCY ) {
			if ( genomeDCJ[ i ]->y > 0 )
				k = 2 * genomeDCJ[ i ]->y - 1;
			else // genomeDCJ[ i ]->y < 0
				k = 2 * abs( genomeDCJ[ i ]->y ) - 2;

			inverseGenome[ k ] = i;
		} 
	}
}	

/* NOTE: the application of a DCJ operation is according to the
	definition 1 that appears in the Bergeron's paper: 
	"A unifying view of genome rearrangements" */
/* Parameters:
	TreeNodePtr nodePtr : contains the genome  
	int i, int j 		: the two positions for applying DCJ
	int firstForm 		: in case TRUE is applied the first form,
		and in case FALSE is applied the second form. 
		The case (a) and case (b) have two forms for applying DCJ.
*/
void applyDCJ( PointDCJPtr *genomeDCJ, int *numPointsDCJ, int i, int j, int firstForm ) 
{
	int k;
	PointDCJ tempU, tempV;

	tempU.x = genomeDCJ[i]->x;
	tempU.y = genomeDCJ[i]->y;
	tempV.x = genomeDCJ[j]->x;
	tempV.y = genomeDCJ[j]->y;

	if ( i != j ) {
		if ( genomeDCJ[ i ]->type == ADJACENCY && 
				genomeDCJ[ j ]->type == ADJACENCY ) {
			/* CASE (a) */
			if ( firstForm == TRUE ) {
				genomeDCJ[i]->x = tempU.x;
				genomeDCJ[i]->y = tempV.x;
				genomeDCJ[j]->x = tempV.y;
				genomeDCJ[j]->y = tempU.y;
			}
			else { // secondForm
				genomeDCJ[i]->x = tempU.x;
				genomeDCJ[i]->y = tempV.y;
				genomeDCJ[j]->x = tempU.y;
				genomeDCJ[j]->y = tempV.x;
			}
		}
		else if ( genomeDCJ[i]->type == ADJACENCY && 
					genomeDCJ[j]->type == TELOMERE ) {
			/* CASE (b) */
			if ( firstForm == TRUE ) {
				genomeDCJ[i]->x = tempU.x;
				genomeDCJ[i]->y = tempV.x;
				genomeDCJ[j]->x = tempU.y;
				genomeDCJ[j]->y = tempU.y;			
			}
			else { // secondForm
				genomeDCJ[i]->x = tempU.y;
				genomeDCJ[i]->y = tempV.x;
				genomeDCJ[j]->x = tempU.x;
				genomeDCJ[j]->y = tempU.x;
			}
		}
		else if ( genomeDCJ[i]->type == TELOMERE && 
					genomeDCJ[j]->type == ADJACENCY ) {
			/* CASE (b) */
			if ( firstForm == TRUE ) {
				genomeDCJ[i]->x = tempV.y;
				genomeDCJ[i]->y = tempV.y;
				genomeDCJ[j]->x = tempV.x;
				genomeDCJ[j]->y = tempU.x;
			}
			else { // secondForm
				genomeDCJ[i]->x = tempV.x;
				genomeDCJ[i]->y = tempV.x;
				genomeDCJ[j]->x = tempV.y;
				genomeDCJ[j]->y = tempU.x;
			}
		}
		else { // genomeDCJ[i]->type == TELOMERE && genomeDCJ[j] == TELOMERE
			/* CASE (c) */
			genomeDCJ[i]->x = tempU.x;
			genomeDCJ[i]->y = tempV.x;
			genomeDCJ[i]->type = ADJACENCY;

			/* eliminate genomeDCJ[j] by shifting elements to the left from position j + 1 to the end */
			for ( k = j + 1; k < (*numPointsDCJ); k++ ) {
				genomeDCJ[ k - 1 ]->x = genomeDCJ[ k ]->x;
				genomeDCJ[ k - 1 ]->y = genomeDCJ[ k ]->y;
				genomeDCJ[ k - 1 ]->type = genomeDCJ[ k ]->type; 
			}
			(*numPointsDCJ)--;
		}	
	}
	else { // i == j
		/* INVERSE of CASE (c) */
		if ( genomeDCJ[i]->type == ADJACENCY ) {
			genomeDCJ[i]->x = tempU.x;
			genomeDCJ[i]->y = tempU.x;
			genomeDCJ[i]->type = TELOMERE;
			
			genomeDCJ[ (*numPointsDCJ) ]->x = tempU.y;
			genomeDCJ[ (*numPointsDCJ) ]->y = tempU.y;
			genomeDCJ[ (*numPointsDCJ) ]->type = TELOMERE;
			(*numPointsDCJ)++;
		}
		else {
			fprintf( stderr, " stderr: can't apply DCJ on just one TELOMERE.\n" );
        	exit( EXIT_FAILURE );
		}
	}

}















