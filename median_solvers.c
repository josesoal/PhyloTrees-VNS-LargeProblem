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

#include "median_solvers.h"
#include "convert.h"
#include "cheaptsp.h"
#include "bbtsp.h"
#include "inversion_median_alberto.h"

/* Call CONVERT2_TO_TSP method 
* inputs: numberGenes, genome1, genome2, genome3
* outputs: adjacencyList
*/
void createTspInstance(int numberGenes, int *genome1, int *genome2, int *genome3, 
    struct adj_struct *adjacencyList)
{
    int CIRCULAR = 1; // true
    struct genome_struct g1;
    struct genome_struct g2;
    struct genome_struct g3;
    struct adj_struct *adj_pool;

    /* allocate memory */ 
    adj_pool = malloc ( ( 14 * numberGenes ) * sizeof ( struct adj_struct ) );
    if ( adj_pool == NULL )
        fprintf ( stderr, "ERROR: adj_pool NULL\n" );

    g1.genes = malloc (numberGenes * sizeof(int));
    if ( g1.genes == NULL )
        fprintf ( stderr, "stderr: var g1.genes, no memory avaliable.\n" );

    g2.genes = malloc (numberGenes * sizeof(int));
    if ( g2.genes == NULL )
        fprintf ( stderr, "stderr: var g2.genes, no memory avaliable.\n" );

    g3.genes = malloc (numberGenes * sizeof(int));
    if ( g3.genes == NULL )
        fprintf ( stderr, "stderr: var g3.genes, no memory avaliable.\n" );

    /* copy genomes */
    int i;
    for (i=0; i<numberGenes; i++){
        g1.genes[i] = genome1[i];
        g2.genes[i] = genome2[i];
        g3.genes[i] = genome3[i];
    }
    /* convert2_to_tsp method */
    convert2_to_tsp (&g1, &g2, &g3, adjacencyList, 
        adj_pool, numberGenes, CIRCULAR); //--method from  convert.c (GRAPPA)

    /* free memory */
    free(adj_pool);
    free(g1.genes);
    free(g2.genes);
    free(g3.genes);
}

/* Call COALESTSP method
* inputs: numberGenes, genome1, genome2, genome3, adjList
* ouputs: medianGenome, tspCycle
* return score of tspCycle 
*/
int callCoalestsp(int numberGenes, int *genome1, int *genome2, int *genome3, 
    struct adj_struct *adjList, int *medianGenome, int *tspCycle, int circular)
{
	/* int ncount: (input) NUMBER of GENES (must be multiplied by 2) 
	* int *tour: (output) MEDIAN GENOME related with outcycle
	* int use_median: NOT USED---
	* int *g1, int *g2, int *g3: (input) INPUT GENOMES
    * struct adj_struct *adj_list: 
	*	(input) ARRAY of adj_struct with 2*genes + 1
	*	firt must be called convert2_to_tsp
    * intpair_t * neighbors: (empty, *not sure*) ARRAY of intpair_t with 2*genes + 1
    * int *stack: (empty) ARRAY with 2*genes + 1 elements
    * int *outcycle: (output) SOLUTION, a cycle for TSP, ARRAY with 2*genes + 1 elements
    * int *degree : (empty) ARRAY with 2*genes + 1 elements
    * int *otherEnd: (empty) ARRAY with 2*genes + 1 elements
    * edge_t * edges : (empty) ARRAY of edge_t with 7*genes
    * int CIRCULAR: NOT USED---
	*/

    int use_median = 0; //not used
    intpair_t *neighbors; 
    int *stack; 
    int *outcycle; 
    int *degree; 
    int *otherEnd; 
    edge_t *edges;

    /* allocate memory */
    neighbors = malloc ( ( 2*numberGenes + 1 ) * sizeof ( intpair_t ) );
    if ( neighbors == NULL )
        fprintf ( stderr, "ERROR: neighbors NULL\n" );

    stack = malloc ( ( 2*numberGenes + 1 ) * sizeof ( int ) );
    if ( stack == NULL )
        fprintf ( stderr, "ERROR: stack NULL\n" );

    outcycle = malloc ( ( 2*numberGenes + 1 ) * sizeof ( int ) );
    if ( outcycle == NULL )
        fprintf ( stderr, "ERROR: outcycle NULL\n" );

    degree = malloc ( ( 2*numberGenes + 1 ) * sizeof ( int ) );
    if ( degree == NULL )
        fprintf ( stderr, "ERROR: degree NULL\n" );

    otherEnd = malloc ( ( 2*numberGenes + 1 ) * sizeof ( int ) );
    if ( otherEnd == NULL )
        fprintf ( stderr, "ERROR: otherEnd NULL\n" );

    edges = malloc ( ( 7*numberGenes ) * sizeof ( edge_t ) );
    if ( edges == NULL )
        fprintf ( stderr, "ERROR: edges NULL\n" );

    /* coalestsp method */
	int score = coalestsp (numberGenes*2, medianGenome, use_median, genome1, genome2, 
        genome3, adjList, neighbors, stack, tspCycle, degree, otherEnd,
        edges, circular); //--method from (GRAPPA)

    /* free memory */
    free(neighbors);
    free(stack);
    free(outcycle);
    free(degree);
    free(otherEnd);
    free(edges);

    return score;
}

/* Call BBTSP method
* inputs: numberGenes, genome1, genome2, genome3, adjList
* ouputs: medianGenome, tspCycle
* return score of tspCycle 
*/

int callBbtsp(int numberGenes, int *genome1, int *genome2, int *genome3, 
    struct adj_struct *adjList, int *medianGenome, int *tspCycle, int circular){
    /* int ncount: (input) NUMBER of GENES (must be multiplied by 2) 
    * int *tour: (output) MEDIAN GENOME related with outcycle
    * int use_median: (input) true (1)
    * int *g1, int *g2, int *g3: (input) INPUT GENOMES
    * struct adj_struct *adj_list: 
    *   (input) ARRAY of adj_struct with 2*genes + 1
    *   firt must be called convert2_to_tsp
    * intpair_t * neighbors: (empty, *not sure*) ARRAY of intpair_t with 2*genes + 1
    * int *stack: (empty) ARRAY with 2*genes + 1 elements
    * int *outcycle: (output) SOLUTION, a cycle for TSP, ARRAY with 2*genes + 1 elements
    * int *degree : (empty) ARRAY with 2*genes + 1 elements
    * int *otherEnd: (empty) ARRAY with 2*genes + 1 elements
    * edge_t * edges : (empty) ARRAY of edge_t with 7*genes
    * int CIRCULAR: (input)
    */

    int use_median = TRUE; 
    intpair_t *neighbors; 
    int *stack; 
    int *outcycle; 
    int *degree; 
    int *otherEnd; 
    edge_t *edges;

    /* allocate memory */
    neighbors = malloc ( ( 2*numberGenes + 1 ) * sizeof ( intpair_t ) );
    if ( neighbors == NULL )
        fprintf ( stderr, "ERROR: neighbors NULL\n" );

    stack = malloc ( ( 2*numberGenes + 1 ) * sizeof ( int ) );
    if ( stack == NULL )
        fprintf ( stderr, "ERROR: stack NULL\n" );

    outcycle = malloc ( ( 2*numberGenes + 1 ) * sizeof ( int ) );
    if ( outcycle == NULL )
        fprintf ( stderr, "ERROR: outcycle NULL\n" );

    degree = malloc ( ( 2*numberGenes + 1 ) * sizeof ( int ) );
    if ( degree == NULL )
        fprintf ( stderr, "ERROR: degree NULL\n" );

    otherEnd = malloc ( ( 2*numberGenes + 1 ) * sizeof ( int ) );
    if ( otherEnd == NULL )
        fprintf ( stderr, "ERROR: otherEnd NULL\n" );

    edges = malloc ( ( 7*numberGenes ) * sizeof ( edge_t ) );
    if ( edges == NULL )
        fprintf ( stderr, "ERROR: edges NULL\n" );

    /* bbtsp method */
    int score =bbtsp ( numberGenes*2, medianGenome, use_median, 
        genome1, genome2, genome3, adjList, neighbors, stack, 
        outcycle, degree, otherEnd, edges, circular ); //--method from (GRAPPA)

    /* free memory */
    free(neighbors);
    free(stack);
    free(outcycle);
    free(degree);
    free(otherEnd);
    free(edges);

    return score;
}

void callCapraraInversionMedian(int numberGenes,
    int *genome1, int *genome2, int *genome3, int *medianGenome, int circular)
{
    /*Begin allocate distmem*/
    int NUM_GENES = numberGenes;
    distmem_t distmem;
    distmem.hammingArr =
    ( int * ) malloc ( ( NUM_GENES + 1 ) * 2 * sizeof ( int ) );
    if ( distmem.hammingArr == ( int * ) NULL )
        fprintf ( stderr, "ERROR: hammingArr NULL\n" );
    distmem.perm1 =
    ( int * ) malloc ( ( 2 * NUM_GENES + 2 ) * sizeof ( int ) );
    if ( distmem.perm1 == ( int * ) NULL )
        fprintf ( stderr, "ERROR: perm1 NULL\n" );
    distmem.perm2 =
    ( int * ) malloc ( ( 2 * NUM_GENES + 2 ) * sizeof ( int ) );
    if ( distmem.perm2 == ( int * ) NULL )
        fprintf ( stderr, "ERROR: perm2 NULL\n" );
    distmem.perm =
    ( int * ) malloc ( ( 2 * NUM_GENES + 2 ) * sizeof ( int ) );
    if ( distmem.perm == ( int * ) NULL )
        fprintf ( stderr, "ERROR: perm NULL\n" );
    distmem.done =
    ( int * ) malloc ( ( 2 * NUM_GENES + 2 ) * sizeof ( int ) );
    if ( distmem.done == ( int * ) NULL )
        fprintf ( stderr, "ERROR: done NULL\n" );
    distmem.greyEdges =
    ( int * ) malloc ( ( 2 * NUM_GENES + 2 ) * sizeof ( int ) );
    if ( distmem.greyEdges == ( int * ) NULL )
        fprintf ( stderr, "ERROR: greyEdges NULL\n" );
    distmem.stack =
    ( int * ) malloc ( ( 2 * NUM_GENES + 2 ) * sizeof ( int ) );
    if ( distmem.stack == ( int * ) NULL )
        fprintf ( stderr, "ERROR: stack NULL\n" );
    distmem.oriented =
    ( int * ) malloc ( ( 2 * NUM_GENES + 2 ) * sizeof ( int ) );
    if ( distmem.oriented == ( int * ) NULL )
        fprintf ( stderr, "ERROR: oriented NULL\n" );
    distmem.cc = ( int * ) malloc ( ( 2 * NUM_GENES + 2 ) * sizeof ( int ) );
    if ( distmem.cc == ( int * ) NULL )
        fprintf ( stderr, "ERROR: cc NULL\n" );
    distmem.labeled =
    ( int * ) malloc ( ( 2 * NUM_GENES + 2 ) * sizeof ( int ) );
    if ( distmem.labeled == ( int * ) NULL )
        fprintf ( stderr, "ERROR: labeled NULL\n" );
    distmem.components = ( component_t * )
    malloc ( ( 2 * NUM_GENES + 2 ) * sizeof ( component_t ) );
    if ( distmem.components == ( component_t * ) NULL )
        fprintf ( stderr, "ERROR: components NULL\n" );
    distmem.uf = UFalloc ( 2 * NUM_GENES + 2 );
    /*End allocate distmem*/
    
    struct genome_struct g1,g2,g3;
    struct genome_struct *genomes[3];
    g1.genes = genome1;
    g2.genes = genome2;
    g3.genes = genome3;
    genomes[0] = &g1;
    genomes[1] = &g2;
    genomes[2] = &g3;
    
    init_global_variables (numberGenes, &distmem );
    if ( circular ){
        albert_inversion_median_circular ( genomes, numberGenes, medianGenome );
    }
    else{
        albert_inversion_median_noncircular ( gen, numberGenes, medianGenome );
    }
    
    /*Begin free distmem*/
    UFfree ( distmem.uf );
    free ( distmem.components );
    free ( distmem.labeled );
    free ( distmem.cc );
    free ( distmem.oriented );
    free ( distmem.stack );
    free ( distmem.greyEdges );
    free ( distmem.done );
    free ( distmem.perm1 );
    free ( distmem.perm2 );
    free ( distmem.perm );
    free ( distmem.hammingArr );
    /*End free distmem*/
}