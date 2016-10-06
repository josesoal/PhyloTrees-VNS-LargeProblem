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
#include <string.h>

#include "auxiliary.h"
#include "tree.h"
#include "random.h"
#include "iterate_tree.h"
#include "dcjdist.h"
#include "stack_array.h" 
#include "int_queue.h"

static int countNumDifferentLeaves( RawDatasetPtr rdatasetPtr ); 
static void readGenomes( 
    TreePtr phyloTreePtr, ParametersPtr paramsPtr, RawDatasetPtr rdatasetPtr );
static void readGenomesDCJ( TreePtr phyloTreePtr, RawDatasetPtr rdatasetPtr );
static void countGenomesPerLeaf( TreePtr phyloTreePtr, 
                RawDatasetPtr rdatasetPtr, MultipleLeafPtr **multiple );
static void readGenomesDCJ_MultiLeaf( 
    TreePtr phyloTreePtr, RawDatasetPtr rdatasetPtr, MultipleLeafPtr *multiple );
static void convertRawGenomeIntoDCJ_LeafNode( 
            TreePtr phyloTreePtr, int i, RawDatasetPtr rdatasetPtr, int r );
static void convertRawGenomeIntoDCJ_LeafCandidate( MultipleLeafPtr *multiple, 
                    int i, int m, RawDatasetPtr rdatasetPtr, int r );

static void readNewickFormat( char *newickTree, char *filename );
static int recoverNamesFromNewickFormat( char *newickTree, GNode *nodes );
static void generateGraphOfNodes( 
                        char *newickTree, GNode *nodes, 
                        int numLeaves, int graph[ MAX_NODES ][ MAX_NODES ] );
static void createTopologyFromGraphByBFS( 
                        TreePtr phyloTreePtr, GNode *nodes, 
                        int graph[ MAX_NODES ][ MAX_NODES ], int numLeaves );
static void createTreeTopologyRandomly( TreePtr phyloTreePtr );
static void createTreeRandomLeaf_FirstBestEdge( 
    TreePtr phyloTreePtr, ParametersPtr paramsPtr, 
    MultipleLeafPtr *multiple, TreePtr previousTreePtr, int iteration );
static void createTreeWith3LeavesRandomly( TreePtr phyloTreePtr );
static void selectLeafNodeByName( 
    TreePtr phyloTreePtr, TreeNodePtr *nodePtr, char * name );
static void selectLeafNodeRandomly( TreePtr phyloTreePtr, TreeNodePtr *nodePtr );
static void selectInternalNode( TreePtr phyloTreePtr, TreeNodePtr *nodePtr );
static void selectInternalNodeByIndex( 
                    TreePtr phyloTreePtr, TreeNodePtr *nodePtr, int index );
static TreeNodePtr getLeafNodePointerByName( 
                                            TreePtr phyloTreePtr, char *name );
static TreeNodePtr getInternalNodePointerByIndex(
                                            TreePtr phyloTreePtr, int index );
static void selectEdgeRandomly( TreePtr phyloTreePtr, 
  TreeNodePtr *startNodePtr, TreeNodePtr *endNodePtr );
static void addRandomLeafNodeIntoEdge( TreePtr phyloTreePtr, 
  TreeNodePtr startNodePtr, TreeNodePtr endNodePtr );
static void selectAvaliableLeaf( TreePtr phyloTreePtr, TreeNodePtr *nodePtr );
static void linkNodeIntoEdge( TreePtr phyloTreePtr, TreeNodePtr nodePtr, 
    TreeNodePtr internalNodePtr, TreeNodePtr startNodePtr, TreeNodePtr endNodePtr );
static void unlinkNodeFromEdge( TreePtr phyloTreePtr, TreeNodePtr nodePtr, 
    TreeNodePtr internalNodePtr, TreeNodePtr startNodePtr, TreeNodePtr endNodePtr );
static void writeRecursiveNewickFormat( FILE *filePtr, 
                char *filename, TreeNodePtr currentNodePtr, int option );

void allocateMemoryForNodes( TreePtr phyloTreePtr, ParametersPtr paramsPtr )
{
	int i, j;	

    phyloTreePtr->startingNodePtr = NULL;

	/* allocate memory for array of pointers to tree nodes */
	phyloTreePtr->numberNodes = (2 * phyloTreePtr->numberLeaves) - 2; /*  2N - 2 */
 	phyloTreePtr->nodesPtrArray = malloc(phyloTreePtr->numberNodes * sizeof(TreeNodePtr));
 	if ( phyloTreePtr->nodesPtrArray == NULL ) nomemMessage( "nodesPtrArray" );

    /* calculate avaliable nodes */
    phyloTreePtr->avaliableLeavesNodes = phyloTreePtr->numberLeaves;
    phyloTreePtr->avaliableInternalNodes = 
                phyloTreePtr->numberNodes - phyloTreePtr->numberLeaves;

    /* allocate memory for leaf nodes(N) and internal nodes(N-2) */
    for( i = 0; i < phyloTreePtr->numberNodes; i++ ) {
    	phyloTreePtr->nodesPtrArray[ i ] = malloc( sizeof( TreeNode ) );
		if ( phyloTreePtr->nodesPtrArray[ i ] == NULL ) nomemMessage( "nodesPtrArray[ i ]" );  
        
        if ( i < phyloTreePtr->numberLeaves ) 
            phyloTreePtr->nodesPtrArray[ i ]->type = LEAF_NODE;
        else
            phyloTreePtr->nodesPtrArray[ i ]->type = INTERNAL_NODE; 

    	phyloTreePtr->nodesPtrArray[ i ]->id           = i + 1;
    	phyloTreePtr->nodesPtrArray[ i ]->index        = i; 
    	phyloTreePtr->nodesPtrArray[ i ]->organism     = NULL; 
    	phyloTreePtr->nodesPtrArray[ i ]->genome       = NULL;
        phyloTreePtr->nodesPtrArray[ i ]->genomeDCJ    = NULL;
        phyloTreePtr->nodesPtrArray[ i ]->inverseDCJ   = NULL;
        phyloTreePtr->nodesPtrArray[ i ]->numPointsDCJ = 0;
    	phyloTreePtr->nodesPtrArray[ i ]->edgeWeight   = -1;
    	phyloTreePtr->nodesPtrArray[ i ]->ancestorPtr  = NULL;
    	phyloTreePtr->nodesPtrArray[ i ]->leftDescPtr  = NULL;
    	phyloTreePtr->nodesPtrArray[ i ]->rightDescPtr = NULL;
    	phyloTreePtr->nodesPtrArray[ i ]->avaliable    = TRUE;
    }
    
    /* allocate memory for genomes */
    if ( paramsPtr->distanceType == INVERSION_DIST ) {
        for( i = 0; i < phyloTreePtr->numberNodes; i++ ) {
            phyloTreePtr->nodesPtrArray[ i ]->genome = 
                calloc( phyloTreePtr->numberGenes, sizeof( int ) );
            if ( phyloTreePtr->nodesPtrArray[ i ]->genome == NULL ) 
                nomemMessage( "genome" );
        }
    }
    else if ( paramsPtr->distanceType == DCJ_DIST ) {
        for( i = 0; i < phyloTreePtr->numberNodes; i++ ) {
            phyloTreePtr->nodesPtrArray[ i ]->genomeDCJ = 
                malloc( 2 * phyloTreePtr->numberGenes * sizeof( PointDCJPtr ) );
            if ( phyloTreePtr->nodesPtrArray[ i ]->genomeDCJ == NULL ) 
                nomemMessage( "genomeDCJ" );
            
            for ( j = 0; j < 2 * phyloTreePtr->numberGenes; j++ ) {
                phyloTreePtr->nodesPtrArray[ i ]->genomeDCJ[ j ] = 
                    malloc( sizeof( PointDCJ ) );
                if ( phyloTreePtr->nodesPtrArray[ i ]->genomeDCJ[ j ] == NULL ) 
                    nomemMessage( "genomeDCJ[ j ]" ); 
            }

            phyloTreePtr->nodesPtrArray[ i ]->inverseDCJ = 
                malloc( 2 * phyloTreePtr->numberGenes * sizeof( int ) );
            if ( phyloTreePtr->nodesPtrArray[ i ]->inverseDCJ == NULL )
                nomemMessage( "inverseDCJ" );

            /* this is defined in funct. readGenomesDCJ(...) */
            // phyloTreePtr->nodesPtrArray[ i ]->numPointsDCJ; 
        }
    }
    else {
        fprintf( stderr, " stderr: incorrect distance.\n" );
        exit( EXIT_FAILURE );
    }
}

void freeTree( TreePtr phyloTreePtr, ParametersPtr paramsPtr )
{
	int i, j;
    
    if ( paramsPtr->distanceType == INVERSION_DIST ) {
        for ( i = 0; i < phyloTreePtr->numberNodes; i++ ) {
            free( phyloTreePtr->nodesPtrArray[ i ]->genome );
        }
    }
    else if ( paramsPtr->distanceType == DCJ_DIST ) {
        for ( i = 0; i < phyloTreePtr->numberNodes; i++ ) {
            for ( j = 0; j < 2 * phyloTreePtr->numberGenes; j++ ) {
                free( phyloTreePtr->nodesPtrArray[ i ]->genomeDCJ[ j ] );
            }   
            free( phyloTreePtr->nodesPtrArray[ i ]->genomeDCJ );
            free( phyloTreePtr->nodesPtrArray[ i ]->inverseDCJ );
        }
    }
    else {
        fprintf( stderr, " stderr: incorrect distance [freeTree(...)].\n" );
        exit( EXIT_FAILURE );
    }

    for ( i = 0; i < phyloTreePtr->numberNodes; i++ ) {
        if ( phyloTreePtr->nodesPtrArray[ i ]->type == LEAF_NODE ) {
            free( phyloTreePtr->nodesPtrArray[ i ]->organism );
        }
        free( phyloTreePtr->nodesPtrArray[ i ] );
    }
    free( phyloTreePtr->nodesPtrArray );
}

/* Copy phyloTree2 into phyloTree1 
NOTE:
Before using this procedure, the following MUST be done:
- set the field phylotreePtr->numberLeaves
- set the field phylotreePtr->numberGenes
- call the procedure allocateMemoryForNodes(phylotreePtr)
*/
void copyTreeInto( TreePtr phyloTree1Ptr, 
        TreePtr phyloTree2Ptr, int copyStructure, ParametersPtr paramsPtr )
{
	int i, j;
	int len = 0;

	/* copy information of leaf nodes(N) and internal nodes (N-2) */
    for( i = 0; i < phyloTree1Ptr->numberNodes; i++ ) {
    	phyloTree1Ptr->nodesPtrArray[ i ]->id = phyloTree2Ptr->nodesPtrArray[ i ]->id; 
    	phyloTree1Ptr->nodesPtrArray[ i ]->index = phyloTree2Ptr->nodesPtrArray[ i ]->index; 

    	if ( phyloTree1Ptr->nodesPtrArray[ i ]->type == LEAF_NODE ) {
    		len = strlen( phyloTree2Ptr->nodesPtrArray[ i ]->organism );
    		phyloTree1Ptr->nodesPtrArray[ i ]->organism = malloc( (len+1) * sizeof(char) );
        	strcpy( phyloTree1Ptr->nodesPtrArray[ i ]->organism, 
                        phyloTree2Ptr->nodesPtrArray[ i ]->organism );
    	}
        //Note:type is not copied, it is already setted in nodes mem alloc.
    	phyloTree1Ptr->nodesPtrArray[ i ]->edgeWeight    = phyloTree2Ptr->nodesPtrArray[ i ]->edgeWeight;
    	phyloTree1Ptr->nodesPtrArray[ i ]->avaliable     = phyloTree2Ptr->nodesPtrArray[ i ]->avaliable;
    	phyloTree1Ptr->nodesPtrArray[ i ]->extremity     = phyloTree2Ptr->nodesPtrArray[ i ]->extremity; 
    }

    /* copy information of genomes */
    if ( paramsPtr->distanceType == INVERSION_DIST ) {
        for( i = 0; i < phyloTree1Ptr->numberNodes; i++ ) {
            phyloTree1Ptr->nodesPtrArray[ i ]->genomeDCJ = NULL;
            phyloTree1Ptr->nodesPtrArray[ i ]->numPointsDCJ = -1;

            for( j = 0; j < phyloTree1Ptr->numberGenes; j++ ) {
                phyloTree1Ptr->nodesPtrArray[ i ]->genome[ j ] = 
                            phyloTree2Ptr->nodesPtrArray[ i ]->genome[ j ];
            } 
        }
    }
    else if ( paramsPtr->distanceType == DCJ_DIST ) {
        for( i = 0; i < phyloTree1Ptr->numberNodes; i++ ) {
            phyloTree1Ptr->nodesPtrArray[ i ]->genome = NULL;

            phyloTree1Ptr->nodesPtrArray[ i ]->numPointsDCJ = 
                            phyloTree2Ptr->nodesPtrArray[ i ]->numPointsDCJ;

            for ( j = 0; j < phyloTree1Ptr->nodesPtrArray[ i ]->numPointsDCJ; j++ ) {
                phyloTree1Ptr->nodesPtrArray[ i ]->genomeDCJ[ j ]->x = 
                    phyloTree2Ptr->nodesPtrArray[ i ]->genomeDCJ[ j ]->x; 
                phyloTree1Ptr->nodesPtrArray[ i ]->genomeDCJ[ j ]->y = 
                    phyloTree2Ptr->nodesPtrArray[ i ]->genomeDCJ[ j ]->y;
                phyloTree1Ptr->nodesPtrArray[ i ]->genomeDCJ[ j ]->type = 
                    phyloTree2Ptr->nodesPtrArray[ i ]->genomeDCJ[ j ]->type; 
            }

            for ( j = 0; j < 2 * phyloTree1Ptr->numberGenes; j++ ) {
                phyloTree1Ptr->nodesPtrArray[ i ]->inverseDCJ[ j ] = 
                    phyloTree2Ptr->nodesPtrArray[ i ]->inverseDCJ[ j ]; 
            }
        }
    }
    else {
        fprintf( stderr, " stderr: incorrect distance.\n" );
        exit( EXIT_FAILURE );
    }

    if ( copyStructure == FALSE ) return; 

    /* copy avaliable nodes */
    phyloTree1Ptr->avaliableLeavesNodes     = phyloTree2Ptr->avaliableLeavesNodes;
    phyloTree1Ptr->avaliableInternalNodes   = phyloTree2Ptr->avaliableInternalNodes;

    /* copy the structure of the tree */
    if ( phyloTree2Ptr->startingNodePtr == NULL ) {
        phyloTree1Ptr->startingNodePtr = NULL; 
    }
    else {
        for ( i = 0; i < phyloTree2Ptr->numberLeaves; i++ ) {
            if ( phyloTree2Ptr->nodesPtrArray[ i ]->id == phyloTree2Ptr->startingNodePtr->id )
                break;
        }
        phyloTree1Ptr->startingNodePtr = phyloTree1Ptr->nodesPtrArray[ i ]; 
    }

    for( i = 0; i < phyloTree1Ptr->numberNodes; i++ ) {
    	if ( phyloTree2Ptr->nodesPtrArray[ i ]->ancestorPtr == NULL ) {
    		phyloTree1Ptr->nodesPtrArray[ i ]->ancestorPtr = NULL;
    	}
    	else {
    		j = phyloTree2Ptr->nodesPtrArray[ i ]->ancestorPtr->index;
    		phyloTree1Ptr->nodesPtrArray[ i ]->ancestorPtr = phyloTree1Ptr->nodesPtrArray[ j ];
    	}

    	if ( phyloTree2Ptr->nodesPtrArray[ i ]->leftDescPtr == NULL ) {
    		phyloTree1Ptr->nodesPtrArray[ i ]->leftDescPtr = NULL;
    	}
    	else {
    		j = phyloTree2Ptr->nodesPtrArray[ i ]->leftDescPtr->index;
    		phyloTree1Ptr->nodesPtrArray[ i ]->leftDescPtr = phyloTree1Ptr->nodesPtrArray[ j ];
    	}

    	if ( phyloTree2Ptr->nodesPtrArray[ i ]->rightDescPtr == NULL ){
    		phyloTree1Ptr->nodesPtrArray[ i ]->rightDescPtr = NULL;
    	}
    	else{
    		j = phyloTree2Ptr->nodesPtrArray[ i ]->rightDescPtr->index;
    		phyloTree1Ptr->nodesPtrArray[ i ]->rightDescPtr = phyloTree1Ptr->nodesPtrArray[ j ];
    	}
    }
}

void readNumberLeavesAndGenes( 
    TreePtr phyloTreePtr, ParametersPtr paramsPtr, RawDatasetPtr rdatasetPtr ) 
{
    int numDifferentLeaves;

    numDifferentLeaves = countNumDifferentLeaves( rdatasetPtr );

    if ( numDifferentLeaves < rdatasetPtr->numberGenomes ) {
        paramsPtr->useMultipleGenomesOneLeaf = TRUE;

        phyloTreePtr->numberLeaves = numDifferentLeaves;
        phyloTreePtr->numberGenes = rdatasetPtr->numberGenes;
    }
    else { //one leaf per genome
        paramsPtr->useMultipleGenomesOneLeaf = FALSE;

        phyloTreePtr->numberLeaves = rdatasetPtr->numberGenomes;
        phyloTreePtr->numberGenes = rdatasetPtr->numberGenes;
    }

    /* discard some undesired cases */
    if ( paramsPtr->distanceType == INVERSION_DIST && 
            paramsPtr->useMultipleGenomesOneLeaf == TRUE ) {
        fprintf( stderr, 
            " stderr: the program does not support using the reversal" );
        fprintf( stderr, " distance and alternative multiple genomes for a leaf.\n" );
        exit( EXIT_FAILURE );
    }
}

static int countNumDifferentLeaves( RawDatasetPtr rdatasetPtr ) 
{
    int i, j, numberLeaves, count; 
    int *visited;

    visited = malloc( rdatasetPtr->numberGenomes * sizeof( int ) );
    if ( visited == NULL )
        nomemMessage( "visited" );
    for ( i = 0; i < rdatasetPtr->numberGenomes; i++ ) { 
        visited[ i ] = FALSE; 
    }

    /* count number of different leaves */
    numberLeaves = 0;
    for ( i = 0; i < rdatasetPtr->numberGenomes; i++ ) {
        count = TRUE;
        for ( j = 0; j < rdatasetPtr->numberGenomes; j++ ) {
            if ( i != j && 
                strcmp( rdatasetPtr->rgenomes[ i ]->organism,
                rdatasetPtr->rgenomes[ j ]->organism ) == 0 &&
                visited[ j ] == TRUE ) 
            {
                count = FALSE;
                break;
            }
        }
        visited[ i ] = TRUE;
        if ( count == TRUE ) {
            numberLeaves++;
        }              
    }

    free( visited );

    return numberLeaves;
}

/* NOTE: before calling this function, you have to call the functions
    "readNumberLeavesAndGenes" and "allocateMemoryForNodes" */
void readGenomesFromRawData( TreePtr phyloTreePtr, ParametersPtr paramsPtr, 
            RawDatasetPtr rdatasetPtr, MultipleLeafPtr **multiple ) 
{
    if ( paramsPtr->distanceType == INVERSION_DIST ) {
        readGenomes( phyloTreePtr, paramsPtr, rdatasetPtr );
    }
    else if ( paramsPtr->distanceType == DCJ_DIST ) {
        if ( paramsPtr->useMultipleGenomesOneLeaf == TRUE ) {
            countGenomesPerLeaf( phyloTreePtr, rdatasetPtr, multiple );
            allocateMemoryForLeafCandidates( phyloTreePtr, *multiple );
            readGenomesDCJ_MultiLeaf( phyloTreePtr, rdatasetPtr, *multiple );
        }
        else { //use one genome per leaf
            readGenomesDCJ( phyloTreePtr, rdatasetPtr );  
        }
    }
    else {
        fprintf( stderr, 
            " stderr: incorrect distance [readGenomesFromRawData(...)].\n" );
        exit( EXIT_FAILURE );
    }
}

static void readGenomes( 
    TreePtr phyloTreePtr, ParametersPtr paramsPtr, RawDatasetPtr rdatasetPtr )
{
    int i, j, k;
    int endSymbol, previousEndSymbol;

    endSymbol = 0; /* init var */
    previousEndSymbol = 0; /* init var */

    for( i = 0; i < phyloTreePtr->numberLeaves; i++ ) {
        phyloTreePtr->nodesPtrArray[ i ]->organism = rdatasetPtr->rgenomes[ i ]->organism;

        if ( rdatasetPtr->numberChromosomesArray[ i ] > 1 ) {
            fprintf( stderr, " stderr: reversal distance just support single-chromosome genomes.\n" );
            exit( EXIT_FAILURE );
        }
        else {
            k = 0;
            /* read a genome i */
            for ( j = 0; j < rdatasetPtr->rgenomes[ i ]->numberElements; j++ ) {
                if ( rdatasetPtr->rgenomes[ i ]->genome[ j ] == SPLIT ) {
                    endSymbol = rdatasetPtr->rgenomes[ i ]->chromosomeType[ 0 ];
                }
                else {
                    phyloTreePtr->nodesPtrArray[ i ]->genome[ k ] = 
                            rdatasetPtr->rgenomes[ i ]->genome[ j ];
                    k++;
                }    
            }

            if ( i > 0 && ( endSymbol != previousEndSymbol ) ) {
                fprintf( stderr, " stderr: genomes must be all circular (@)" ); 
                fprintf( stderr, " or all linear ($) for reversal distance.\n" ); 
                exit( EXIT_FAILURE );
            }
            previousEndSymbol = endSymbol;
        }
    }

    /* determine if all element of dataset are circular or linear */
    if ( endSymbol == LINEAR_SYM ) paramsPtr->circular = FALSE;
    if ( endSymbol == CIRCULAR_SYM ) paramsPtr->circular = TRUE;
}

static void readGenomesDCJ( TreePtr phyloTreePtr, RawDatasetPtr rdatasetPtr )
{
    int i;

    for( i = 0; i < phyloTreePtr->numberLeaves; i++ ) {
        phyloTreePtr->nodesPtrArray[ i ]->organism = 
                            rdatasetPtr->rgenomes[ i ]->organism;

        convertRawGenomeIntoDCJ_LeafNode( 
                                    phyloTreePtr, i, rdatasetPtr, i );
    }
}

static void countGenomesPerLeaf( TreePtr phyloTreePtr, 
                RawDatasetPtr rdatasetPtr, MultipleLeafPtr **multiple )
{
    int i , j, count, index, advanceIndex;
    int *visited;

    visited = malloc( rdatasetPtr->numberGenomes * sizeof( int ) );
    if ( visited == NULL )
            nomemMessage( "visited" );
    for ( i = 0; i < rdatasetPtr->numberGenomes; i++ ) { 
        visited[ i ] = FALSE; 
    }

    /* allocate memory for array multiple */
    ( *multiple ) = malloc( 
        phyloTreePtr->numberLeaves * sizeof( MultipleLeafPtr ) );
    if ( ( *multiple ) == NULL )
        nomemMessage( "( *multiple )" );
    for ( i = 0; i < phyloTreePtr->numberLeaves; i++ ) {
        ( *multiple )[ i ] = malloc( sizeof( MultipleLeaf ) );
        if ( ( *multiple )[ i ] == NULL ) 
            nomemMessage( "( *multiple )[ i ]" );

        ( *multiple )[ i ]->numLeafCandidates = 0;
    }

    /* count just multiple genomes per leaf */
    index = 0;
    for ( i = 0; i < rdatasetPtr->numberGenomes; i++ ) {
        count = 0;
        advanceIndex = TRUE;
        /* verify that genome "i" is part of a set of multiple genomes */
        for ( j = 0; j < rdatasetPtr->numberGenomes; j++ ) {
            if ( i != j && 
                strcmp( rdatasetPtr->rgenomes[ i ]->organism,
                rdatasetPtr->rgenomes[ j ]->organism ) == 0 )    
            {
                if ( visited[ j ] == FALSE ) {
                    count++;
                    visited[ j ] = TRUE;
                }
                else {
                    advanceIndex = FALSE;
                    break;
                }         
            }
        }

        visited[ i ] = TRUE;

        if ( count > 0 ) {
            /* number of candidates are "count" + 1 ( we use +1, because 
                the genome at "i" index was not count ) */
            ( *multiple )[ index ]->numLeafCandidates = count + 1;
        }

        if ( advanceIndex == TRUE )
            index++;
    }

    free( visited );
    /* Note: multiple is freed in function "freeMultipleLeafs()" */
}

void allocateMemoryForLeafCandidates( 
                TreePtr phyloTreePtr, MultipleLeafPtr *multiple )
{
    int i, j;

    for ( i = 0; i < phyloTreePtr->numberLeaves; i++ ) {

        if ( multiple[ i ]->numLeafCandidates > 0 ) {
            multiple[ i ]->candidates = 
                        malloc( multiple[ i ]->numLeafCandidates * 
                                                sizeof( CandidatePtr ) );
            if ( multiple[ i ]->candidates == NULL )
                nomemMessage( "multiple[ i ]->candidates" );

            for ( j = 0; j < multiple[ i ]->numLeafCandidates; j++) {
                allocateMemForCandidate( 
                    phyloTreePtr, &multiple[ i ]->candidates[ j ] );//from iterate_tree.c
            }
        }
        else {
            multiple[ i ]->candidates = NULL;
        }
    }
}

void freeMultipleLeafs( 
    TreePtr phyloTreePtr, MultipleLeafPtr **multiple, ParametersPtr paramsPtr )
{
    int i, j;
    CandidatePtr someCandidatePtr;

    if ( paramsPtr->useMultipleGenomesOneLeaf == TRUE ){    
        
        for ( i = 0; i < phyloTreePtr->numberLeaves; i++ ) {
            if ( ( *multiple )[ i ]->numLeafCandidates > 0 ) {
                for ( j = 0; j < ( *multiple )[ i ]->numLeafCandidates; j++) {
                    someCandidatePtr = ( *multiple )[ i ]->candidates[ j ];
                    freeCandidate( phyloTreePtr, &someCandidatePtr );//from iterate_tree.c
                }
                free( ( *multiple )[ i ]->candidates );                
            }
            free( ( *multiple )[ i ] ); 
        }
        free( ( *multiple ) );
    } 
}

static void readGenomesDCJ_MultiLeaf( 
    TreePtr phyloTreePtr, RawDatasetPtr rdatasetPtr, MultipleLeafPtr *multiple )
{
    int i, j, k, h, unique, index, advanceIndex;
    int *visited;

    visited = malloc( rdatasetPtr->numberGenomes * sizeof( int ) );
    if ( visited == NULL )
            nomemMessage( "visited" );
    for ( i = 0; i < rdatasetPtr->numberGenomes; i++ ) { 
        visited[ i ] = FALSE; 
    }

    index = 0;
    k = 0; //index of phyloTreePtr
    for ( i = 0; i < rdatasetPtr->numberGenomes; i++ ) {
        unique = TRUE;
        advanceIndex = TRUE;
        h = 0;

        for ( j = 0; j < rdatasetPtr->numberGenomes; j++ ) {

            if ( i != j && 
                strcmp( rdatasetPtr->rgenomes[ i ]->organism,
                rdatasetPtr->rgenomes[ j ]->organism ) == 0 )    
            {
                if ( visited[ j ] == FALSE ) {
                    if ( h == 0 ) {
                        convertRawGenomeIntoDCJ_LeafCandidate(
                                        multiple, index, h, rdatasetPtr, i );
                        h++;
                    }

                    convertRawGenomeIntoDCJ_LeafCandidate(
                                        multiple, index, h, rdatasetPtr, j );
                    h++;

                    visited[ j ] = TRUE;
                }
                else {
                    unique = FALSE;
                    advanceIndex = FALSE;
                    break;
                }
            }
        }

        visited[ i ] = TRUE;

        if ( unique == TRUE ) {
            convertRawGenomeIntoDCJ_LeafNode( 
                                    phyloTreePtr, k, rdatasetPtr, i );
            k++;
        }

        if ( advanceIndex == TRUE )
            index++;
    }

    free( visited );
}

static void convertRawGenomeIntoDCJ_LeafNode( 
            TreePtr phyloTreePtr, int i, RawDatasetPtr rdatasetPtr, int r )
{
    int j, k, h, init;
    int temp, a, b, num, newChromosome;

    phyloTreePtr->nodesPtrArray[ i ]->organism = 
                            rdatasetPtr->rgenomes[ r ]->organism;

    init = 0; k = 1; h = 0;
    temp = 0; a = 0; b = 0;
    newChromosome = TRUE;
    /* convert raw genome into dcj format */
    for ( j = 0; j < rdatasetPtr->rgenomes[ r ]->numberElements; j++ ) {
        /* NOTE: how to encode a gene either positive or negative
        gen a: tail (-a), head (+a)
        gen -a: head (+a), tail (-a)*/ 

        num = rdatasetPtr->rgenomes[ r ]->genome[ j ];

        if ( newChromosome == TRUE ) {
            temp = -1 * num;
            a = num;
            newChromosome = FALSE;
        }
        else if ( num == SPLIT ) { 
            if ( rdatasetPtr->rgenomes[ r ]->chromosomeType[ h ] == LINEAR_SYM ) {
                phyloTreePtr->nodesPtrArray[ i ]->genomeDCJ[ init ]->x = temp;                         
                phyloTreePtr->nodesPtrArray[ i ]->genomeDCJ[ init ]->y = temp;
                phyloTreePtr->nodesPtrArray[ i ]->genomeDCJ[ init ]->type = TELOMERE;
                        
                phyloTreePtr->nodesPtrArray[ i ]->genomeDCJ[ k ]->x = a;                         
                phyloTreePtr->nodesPtrArray[ i ]->genomeDCJ[ k ]->y = a;
                phyloTreePtr->nodesPtrArray[ i ]->genomeDCJ[ k ]->type = TELOMERE;
                k++; init = k;
                k++;
                newChromosome = TRUE;
            } 
            else if ( rdatasetPtr->rgenomes[ r ]->chromosomeType[ h ] == CIRCULAR_SYM ) {
                phyloTreePtr->nodesPtrArray[ i ]->genomeDCJ[ init ]->x = a;                         
                phyloTreePtr->nodesPtrArray[ i ]->genomeDCJ[ init ]->y = temp;
                phyloTreePtr->nodesPtrArray[ i ]->genomeDCJ[ init ]->type = ADJACENCY;
                init = k;
                k++;
                newChromosome = TRUE;
            }
            h++;
        }
        else { // num is a number
            b = -1 * num;
            phyloTreePtr->nodesPtrArray[ i ]->genomeDCJ[ k ]->x = a;                         
            phyloTreePtr->nodesPtrArray[ i ]->genomeDCJ[ k ]->y = b;
            phyloTreePtr->nodesPtrArray[ i ]->genomeDCJ[ k ]->type = ADJACENCY;
            a = num;
            k++;
        }

    }//end-for
    phyloTreePtr->nodesPtrArray[ i ]->numPointsDCJ = k - 1;

    /* generate inverse of dcj genome */
    calculateInverseGenome( 
        phyloTreePtr->nodesPtrArray[ i ]->genomeDCJ, 
        phyloTreePtr->nodesPtrArray[ i ]->numPointsDCJ, 
        phyloTreePtr->nodesPtrArray[ i ]->inverseDCJ );//--from dcjdist.c
}

static void convertRawGenomeIntoDCJ_LeafCandidate( MultipleLeafPtr *multiple, 
                    int i, int m, RawDatasetPtr rdatasetPtr, int r )
{
    int j, k, h, init;
    int temp, a, b, num, newChromosome;

    init = 0; k = 1; h = 0;
    temp = 0; a = 0; b = 0;
    newChromosome = TRUE;
    /* convert raw genome into dcj format */
    for ( j = 0; j < rdatasetPtr->rgenomes[ r ]->numberElements; j++ ) {
        /* NOTE: how to encode a gene either positive or negative
        gen a: tail (-a), head (+a)
        gen -a: head (+a), tail (-a)*/ 

        num = rdatasetPtr->rgenomes[ r ]->genome[ j ];

        if ( newChromosome == TRUE ) {
            temp = -1 * num;
            a = num;
            newChromosome = FALSE;
        }
        else if ( num == SPLIT ) { 
            if ( rdatasetPtr->rgenomes[ r ]->chromosomeType[ h ] == LINEAR_SYM ) {
                multiple[ i ]->candidates[ m ]->genomeDCJ[ init ]->x = temp;                         
                multiple[ i ]->candidates[ m ]->genomeDCJ[ init ]->y = temp;
                multiple[ i ]->candidates[ m ]->genomeDCJ[ init ]->type = TELOMERE;
                        
                multiple[ i ]->candidates[ m ]->genomeDCJ[ k ]->x = a;                         
                multiple[ i ]->candidates[ m ]->genomeDCJ[ k ]->y = a;
                multiple[ i ]->candidates[ m ]->genomeDCJ[ k ]->type = TELOMERE;
                k++; init = k;
                k++;
                newChromosome = TRUE;
            } 
            else if ( rdatasetPtr->rgenomes[ r ]->chromosomeType[ h ] == CIRCULAR_SYM ) {
                multiple[ i ]->candidates[ m ]->genomeDCJ[ init ]->x = a;                         
                multiple[ i ]->candidates[ m ]->genomeDCJ[ init ]->y = temp;
                multiple[ i ]->candidates[ m ]->genomeDCJ[ init ]->type = ADJACENCY;
                init = k;
                k++;
                newChromosome = TRUE;
            }
            h++;

        }
        else { // num is a number
            b = -1 * num;
            multiple[ i ]->candidates[ m ]->genomeDCJ[ k ]->x = a;                         
            multiple[ i ]->candidates[ m ]->genomeDCJ[ k ]->y = b;
            multiple[ i ]->candidates[ m ]->genomeDCJ[ k ]->type = ADJACENCY;
            a = num;
            k++;
        }

    }//end-for
    multiple[ i ]->candidates[ m ]->numPointsDCJ = k - 1;

    /* generate inverse of dcj genome */
    calculateInverseGenome( 
        multiple[ i ]->candidates[ m ]->genomeDCJ, 
        multiple[ i ]->candidates[ m ]->numPointsDCJ, 
        multiple[ i ]->candidates[ m ]->inverseDCJ );//--from dcjdist.c
}

/* IMPORTANT NOTE: This function is used just for the Small-Phylogny case */
void createTopologyFromNewickFormat( TreePtr phyloTreePtr, ParametersPtr paramsPtr )
{
    /* Note: constants MAX_NEWICK_LEN  and MAX_NODES are in stack_array.h*/
    char newickTree[ MAX_NEWICK_LEN ];
    int numNodes, numLeaves, i, j;
    GNode nodes[ MAX_NODES ]; 

    readNewickFormat( newickTree,  paramsPtr->newickFile );
    numLeaves = recoverNamesFromNewickFormat( newickTree, nodes );
    numNodes = numLeaves + ( numLeaves - 2 ); /* leaves + internal nodes */

    int graph[ MAX_NODES ][ MAX_NODES ];
    for ( i = 0; i < numNodes; i++ ) {
        for ( j = 0; j < numNodes; j++ ) {
            graph[ i ][ j ] = 0;
        }
    }

    generateGraphOfNodes( newickTree, nodes, numLeaves, graph );
    createTopologyFromGraphByBFS( phyloTreePtr, nodes, graph, numLeaves );
}

/* IMPORTANT NOTE: The input is an unrooted phylogenetic tree, where 
    each internal node has exactly 3 adjacent nodes. Given the above
    constraints, the newick format of the tree must be written in a  way 
    that is always considered a pair between "(" and ")". For example 
    ((A,B),(B,C)) and (A,((B,C),D)) are valid.   */
static void readNewickFormat( char *newickTree, char *filename )
{
    FILE *filePtr;
    int c; /* use int (not char) for the EOF */
    int k;
    char buffer[ MAX_NEWICK_LEN ];

    if ( ( filePtr = fopen( filename, "r" ) ) == NULL ) {
        fprintf( stderr, " stderr: %s file could not be opened.\n", filename );
        exit( EXIT_FAILURE );
    }
    else {
            /* read newick format into buffer */
            k = 0;
            while ( ( c = fgetc( filePtr ) ) != EOF ) {
                if ( c == '\n' ){
                    buffer[ k ] = '\0';
                    break;
                }
                else{
                    buffer[ k ] = c;
                    k++;
                    if ( k + 1 > MAX_NEWICK_LEN ) {
                        fprintf( stderr, " stderr: increment MAX_NEWICK_LEN value!\n" );
                        exit( EXIT_FAILURE );
                    }
                }
            }
            strcpy( newickTree, buffer ); 
    }
    fclose( filePtr );
}

/* return number of leaves in the newick format */
static int recoverNamesFromNewickFormat( char *newickTree, GNode *nodes ) 
{
    int i, j , k, len;
    char buffer[ MAX_STRING_LEN ];

    i = 0; 
    j = 0; 
    k = 0;
    len = strlen( newickTree );
    while ( i < len ) {
        if ( newickTree[ i ] == '(' || 
                    newickTree[ i ] == ')' || 
                    newickTree[ i ] == ','  ) {
            if ( k > 0 ) {
                buffer[ k ] = '\0';
                strcpy( nodes[ j ].name, buffer );
                k = 0;
                j++;
            }
        }
        else {
            buffer[ k ] = newickTree[ i ];
            k++;
            if ( k + 1 > MAX_STRING_LEN) {
                fprintf( stderr, " stderr: increment MAX_STRING_LEN value!\n" );
                exit( EXIT_FAILURE );
            }
        }
        i++;
    }

    if ( j > MAX_LEAVES ) {
        fprintf( stderr, " stderr: increment MAX_LEAVES value!\n" );
        exit( EXIT_FAILURE );
    }

    return j;
}

static void generateGraphOfNodes( 
    char *newickTree, GNode *nodes, int numLeaves, int graph[ MAX_NODES ][ MAX_NODES ] )
{
    Stack parenthesisSt;
    Stack namesSt;
    char *c1, *c2;
    int i, j, k, len;
    int numNodes, indexInternal, index_c1, index_c2;
    char buffer[ MAX_STRING_LEN ];

    numNodes = numLeaves + ( numLeaves - 2 ); /* leaves + internal nodes */
    indexInternal = numLeaves; /* internal node starts at numLeaves */
    len = strlen( newickTree );
    parenthesisSt.top = -1; /* init empty stack */
    namesSt.top = -1; /* init empty stack */
 
    i = 0; 
    k = 0;
    buffer[ 0 ] = newickTree[ i ];
    buffer[ 1 ] = '\0';
    push( &parenthesisSt, buffer );
    i++;
    while ( ! isStackEmpty( &parenthesisSt ) && i < len ) {
        /* put parenthesis and names into stacks */
        if ( newickTree[ i ] == '(' || newickTree[ i ] == ')' ) {
            if ( k > 0 ) {
                buffer[ k ] = '\0';
                push( &namesSt, buffer );
                k = 0; 
            }
            buffer[ 0 ] = newickTree[ i ];
            buffer[ 1 ] = '\0';
            push( &parenthesisSt, buffer );
        }
        else if ( newickTree[ i ] == ',' ) {
            if ( k > 0 ) {
                buffer[ k ] = '\0';
                push( &namesSt, buffer );
                k = 0;
            }
        }
        else {
            buffer[ k ] = newickTree[ i ];
            k++;
        }
       
        /* verify if parenthesis stack have "(" and ")" */
        char str[ 3 ];
        if ( strcmp( parenthesisSt.s[ parenthesisSt.top - 1 ], "(" ) == 0 
            && strcmp( parenthesisSt.s[ parenthesisSt.top ], ")" ) == 0 ) 
        {
            pop( &parenthesisSt );
            pop( &parenthesisSt );
            c1 = pop( &namesSt );
            c2 = pop( &namesSt );

            /* find the index of c1 and c2 */
            index_c1 = -1;
            index_c2 = -1;
            for ( j = 0; j< numNodes ;j++ ) {
                if ( strcmp( c1, nodes[ j ].name ) == 0 ) {
                    index_c1 = j;
                }
                if ( strcmp( c2, nodes[ j ].name ) == 0 ) {
                    index_c2 = j;         
                }
                if ( index_c1 != -1 && index_c2 != -1 ) {
                    break;
                }
            }

            /* create a relation in the graph */
            if ( isStackEmpty( &parenthesisSt ) ) {
                /* stack is empty create a relation between c1 and c2 */
                graph[ index_c1 ][ index_c2 ] = 1;
                graph[ index_c2 ][ index_c1 ] = 1; 
            }
            else {
                /* stack is not empty, create a relation with an internal node */
                graph[ index_c1 ][ indexInternal ] = 1;
                graph[ indexInternal ][ index_c1 ] = 1;
                graph[ index_c2 ][ indexInternal ] = 1;
                graph[ indexInternal ][ index_c2 ] = 1;
                /* create a name for the internal node and push into names stack */
                str[ 0 ] = 'i';
                str[ 1 ] = indexInternal + '0';
                str[ 2 ] = '\0';
                strcpy( nodes[ indexInternal ].name, str );
                push( &namesSt, str );
                indexInternal++;
            }
        }
        i++;
    }
}

static void createTopologyFromGraphByBFS( 
                        TreePtr phyloTreePtr, GNode *nodes, 
                        int graph[ MAX_NODES ][ MAX_NODES ], int numLeaves )
{   
    int i, j, index, numNodes;
    IntQueue iqueue;
    TreeNodePtr nodePtr, node1Ptr, node2Ptr, internalPtr;
    
    iqueue.headPtr = NULL;
    iqueue.tailPtr = NULL;
    numNodes = numLeaves + ( numLeaves - 2 );
    for ( i = 0; i < numNodes; i++ ){
        nodes[ i ].visited = FALSE;
    }
    
    /*  first leaf (nodes[0]) is the root node */
    nodes[ 0 ].visited = TRUE;
    selectLeafNodeByName( phyloTreePtr, &nodePtr, nodes[ 0 ].name );
    phyloTreePtr->startingNodePtr = nodePtr;
    nodePtr->ancestorPtr = NULL;
    nodePtr->leftDescPtr = NULL;
    /* create right descendant of root */
    for ( i = numLeaves; i < numNodes; i++ ) {
        if ( graph[ 0 ][ i ] == 1 ) {
            selectInternalNodeByIndex( phyloTreePtr, &internalPtr, i );
            nodePtr->rightDescPtr = internalPtr;
            internalPtr->ancestorPtr = nodePtr;
            /* since nodes[0] is a leaf it only 
                has one internal adjacent node */
            break;
        } 
    }

    /* the right descendant of root is the first element of queue */
    nodes[ i ].visited = TRUE;
    enqueue_i( &iqueue, i );

    while ( ! isEmpty_i( &iqueue ) ) {
        index = dequeue_i( &iqueue );

        /* get the pointer of the node */
        if ( index < numLeaves ) {
            nodePtr = getLeafNodePointerByName( 
                                phyloTreePtr, nodes[ index ].name );
        }
        else {
            nodePtr = getInternalNodePointerByIndex( phyloTreePtr, index );
        }
        /* make left and right descendant NULL */
        nodePtr->leftDescPtr = NULL;
        nodePtr->rightDescPtr = NULL;
        /* create left descendant if exists */
        for ( i = 0; i < numNodes; i++ ) {
            if ( graph[ index ][ i ] == 1 && nodes[ i ].visited == FALSE ) 
            {
                if ( i < numLeaves ) {
                    selectLeafNodeByName( 
                        phyloTreePtr, &node1Ptr, nodes[ i ].name );
                }
                else {
                    selectInternalNodeByIndex( phyloTreePtr, &node1Ptr, i );
                }
                nodePtr->leftDescPtr = node1Ptr;
                node1Ptr->ancestorPtr = nodePtr;

                nodes[ i ].visited = TRUE;
                enqueue_i( &iqueue, i );
                break; 
            }            
        }
        /* create right descendant if exists */
        for ( j = i + 1; j < numNodes; j++ ) {
            if ( graph[ index ][ j ] == 1 && nodes[ j ].visited == FALSE ) 
            {
                if ( j < numLeaves ) {
                    selectLeafNodeByName( 
                        phyloTreePtr, &node2Ptr, nodes[ j ].name );
                }
                else {
                    selectInternalNodeByIndex( phyloTreePtr, &node2Ptr, j );
                }
                nodePtr->rightDescPtr = node2Ptr;
                node2Ptr->ancestorPtr = nodePtr;

                nodes[ j ].visited = TRUE;
                enqueue_i( &iqueue, j );
                break; 
            }            
        }
    }//end-while
}

int createInitialTreeTopology( 
    TreePtr phyloTreePtr, ParametersPtr paramsPtr, 
    MultipleLeafPtr *multiple, TreePtr previousTreePtr, int iteration )
{
    int i, maxIterations, score, newScore;
    Tree temporalTree1, temporalTree2;
    /*NOTE: execute many times the initialization method 
    and choose the best tree generated */
    maxIterations = 3;//10; //at least must be 1

    /* execute one init method */
    if ( paramsPtr->initMethod == R_LEAF_R_EDGE ) {
        createTreeTopologyRandomly( phyloTreePtr );
        score = labelOptimizeTree( phyloTreePtr, paramsPtr, multiple, previousTreePtr, iteration );
    }
    else if ( paramsPtr->initMethod == R_LEAF_1BEST_EDGE ){
        temporalTree1.numberLeaves = phyloTreePtr->numberLeaves;
        temporalTree1.numberGenes = phyloTreePtr->numberGenes;
        allocateMemoryForNodes( &temporalTree1, paramsPtr ); 
        temporalTree2.numberLeaves = phyloTreePtr->numberLeaves;
        temporalTree2.numberGenes = phyloTreePtr->numberGenes;
        allocateMemoryForNodes( &temporalTree2, paramsPtr );

        copyTreeInto( &temporalTree1, phyloTreePtr, TRUE, paramsPtr );//make a fresh copy
        createTreeRandomLeaf_FirstBestEdge( 
            &temporalTree1, paramsPtr, multiple, previousTreePtr, iteration );

        score = labelOptimizeTree( 
            &temporalTree1, paramsPtr, multiple, previousTreePtr, iteration );//from iterate_tree.c

        if ( DEBUG ){ printf( "[Initial tree score: %d]\n",score ); }

        for ( i = 1; i < maxIterations; i++ ) {
            copyTreeInto( &temporalTree2, phyloTreePtr, TRUE, paramsPtr );//make a fresh copy
            createTreeRandomLeaf_FirstBestEdge( 
                &temporalTree2, paramsPtr, multiple, previousTreePtr, iteration );
            newScore = labelOptimizeTree( 
                &temporalTree2, paramsPtr, multiple, previousTreePtr, iteration );//from iterate_tree.c
            if ( DEBUG ){ printf( "[Initial tree score: %d]\n", newScore );}

            if ( newScore < score ) {
                score = newScore;
                copyTreeInto( &temporalTree1, &temporalTree2, TRUE, paramsPtr );
            }
        }

        copyTreeInto( phyloTreePtr, &temporalTree1, TRUE, paramsPtr );

        freeTree( &temporalTree1, paramsPtr );//--method from tree.c
        freeTree( &temporalTree2, paramsPtr );//--method from tree.c

    }
    else {
        fprintf( stderr, " stderr: incorrect initialization method\n" );
        exit( EXIT_FAILURE );
    }

    return score;
}


/* 1st initialization method:
* Create a tree topology by inserting a "random leaf" into a "random edge" */
static void createTreeTopologyRandomly( TreePtr phyloTreePtr )
{
	int i;
	TreeNodePtr startNodePtr;	/* start node of an edge */
	TreeNodePtr endNodePtr;		/* end node of an edge */

    startNodePtr = NULL; /* init var */
    endNodePtr = NULL; /* init var */

	/* create initial tree with 3 leaves  */
	createTreeWith3LeavesRandomly(phyloTreePtr); 
	/* add remaining leaves randomly */
	for(i=3; i < phyloTreePtr->numberLeaves; i++){
		selectEdgeRandomly(phyloTreePtr, &startNodePtr, &endNodePtr);
		addRandomLeafNodeIntoEdge(phyloTreePtr, startNodePtr, endNodePtr);
	}	
}

/* 2nd initialization method:
* Create a tree by inserting a "random leaf" into 
* the edge that leads to the "first best" score */
static void createTreeRandomLeaf_FirstBestEdge( 
        TreePtr phyloTreePtr, ParametersPtr paramsPtr, 
        MultipleLeafPtr *multiple, TreePtr previousTreePtr, int iteration ) 
{
	int i, j, score, newScore;
	TreeNodePtr node1Ptr, node2Ptr, node3Ptr, nodeTmpPtr, internalNodePtr;

    node1Ptr        = NULL; /* init var */ 
    node2Ptr        = NULL; /* init var */  
    node3Ptr        = NULL; /* init var */  
    nodeTmpPtr      = NULL; /* init var */  
    internalNodePtr = NULL; /* init var */ 

    /*[CREATE INITIAL TREE with 3 LEAVES]*/
	/* first select a random leaf as  1st node */
    if ( paramsPtr->useOutgroup == TRUE ) {
        selectLeafNodeByName( phyloTreePtr, &node1Ptr, paramsPtr->outgroup );
    }
    else {
	   selectLeafNodeRandomly( phyloTreePtr, &node1Ptr );//1st random leaf
    }

	/* select the 2nd leaf that makes the first best score*/
	selectAvaliableLeaf( phyloTreePtr, &node2Ptr );
    if ( paramsPtr->distanceType == DCJ_DIST ) {
        score = DCJdistance( node1Ptr->genomeDCJ, node2Ptr->genomeDCJ,
                            node1Ptr->inverseDCJ, node2Ptr->inverseDCJ,     
                            node1Ptr->numPointsDCJ, node2Ptr->numPointsDCJ, 
                            phyloTreePtr->numberGenes );//--from dcjdist.c
    }
    else { // INVERSION_DIST or BREAKPOINT_DIST
	   score = calculateDistance( node1Ptr->genome, node2Ptr->genome, 
                                phyloTreePtr->numberGenes, paramsPtr );//--from iterate_tree.c
    }
	for ( i = 0; i < phyloTreePtr->numberLeaves; i++ ) {
		if ( phyloTreePtr->nodesPtrArray[i]->avaliable == TRUE ) {
			nodeTmpPtr = phyloTreePtr->nodesPtrArray[ i ];
			nodeTmpPtr->avaliable = FALSE; 
			phyloTreePtr->avaliableLeavesNodes--;

            if ( paramsPtr->distanceType == DCJ_DIST ) {
                newScore = DCJdistance( node1Ptr->genomeDCJ, nodeTmpPtr->genomeDCJ,
                            node1Ptr->inverseDCJ, nodeTmpPtr->inverseDCJ, 
                            node1Ptr->numPointsDCJ, nodeTmpPtr->numPointsDCJ, 
                            phyloTreePtr->numberGenes );//--from dcjdist.c
            }
            else { // INVERSION_DIST or BREAKPOINT_DIST
                newScore = calculateDistance( node1Ptr->genome, nodeTmpPtr->genome, 
                            phyloTreePtr->numberGenes, paramsPtr );//--from iterate_tree.c
            }

			if ( newScore < score ) {
				score = newScore;
				node2Ptr->avaliable = TRUE; //let node2Ptr be avaliable
				phyloTreePtr->avaliableLeavesNodes++;
				node2Ptr = nodeTmpPtr; //update pointer
				break;
			}
			else {
				nodeTmpPtr->avaliable = TRUE; //let nodeTmpPtr be avaliable
				phyloTreePtr->avaliableLeavesNodes++;
			}
		}
	}

	/* select the 3rd leaf that leads to the first best score and create initial tree */
	selectInternalNode( phyloTreePtr, &internalNodePtr );
	selectAvaliableLeaf( phyloTreePtr, &node3Ptr );
	phyloTreePtr->startingNodePtr = node1Ptr;
	
	node1Ptr->ancestorPtr = NULL;
	node1Ptr->leftDescPtr = NULL; 
	node1Ptr->rightDescPtr = internalNodePtr;
	internalNodePtr->ancestorPtr = node1Ptr;

	internalNodePtr->leftDescPtr = node2Ptr;
	internalNodePtr->rightDescPtr = node3Ptr;
	node2Ptr->ancestorPtr = internalNodePtr;
	node3Ptr->ancestorPtr = internalNodePtr;

    score = labelOptimizeTree( 
        phyloTreePtr, paramsPtr, multiple, previousTreePtr, iteration );//--from iterate_tree.c 

	for ( i = 0; i < phyloTreePtr->numberLeaves; i++ ) {
		if ( phyloTreePtr->nodesPtrArray[i]->avaliable == TRUE ) {
			nodeTmpPtr = phyloTreePtr->nodesPtrArray[ i ];
			nodeTmpPtr->avaliable = FALSE;
			phyloTreePtr->avaliableLeavesNodes--;

			/* link to the temporal node instead of node 3 */
			internalNodePtr->rightDescPtr = nodeTmpPtr;
			nodeTmpPtr->ancestorPtr = internalNodePtr;
			newScore = labelOptimizeTree( 
                phyloTreePtr, paramsPtr, multiple, previousTreePtr, iteration );//--from iterate_tree.c

			if ( newScore < score ) {
				score = newScore;
				/* undo select node 3 */
				node3Ptr->avaliable = TRUE; 
				node3Ptr->ancestorPtr = NULL;
				phyloTreePtr->avaliableLeavesNodes++;
				/* update pointer */
				node3Ptr = nodeTmpPtr; 
				break;
			}
			else {
				/* recover last state */
				internalNodePtr->rightDescPtr = node3Ptr;
				node3Ptr->ancestorPtr = internalNodePtr;
				/* undo select temporal node */
				nodeTmpPtr->avaliable = TRUE; 
				nodeTmpPtr->ancestorPtr = NULL;
				phyloTreePtr->avaliableLeavesNodes++;
			}
		}
	}

    /* [ADD REMAINING LEAVES INTO TREE] */

	TreeNodePtr startNodePtr, tmpStartNodePtr;	/* start node of an edge */
	TreeNodePtr endNodePtr, tmpEndNodePtr;		/* end node of an edge */
	TreeNodePtr nodePtr;

    startNodePtr    = NULL; /* init var */ 
    tmpStartNodePtr = NULL; /* init var */
    endNodePtr      = NULL; /* init var */ 
    tmpEndNodePtr   = NULL; /* init var */
    nodePtr         = NULL; /* init var */

	for( i = 3; i < phyloTreePtr->numberLeaves; i++ ) {
		selectLeafNodeRandomly( phyloTreePtr, &nodePtr );
		selectInternalNode( phyloTreePtr, &internalNodePtr );

		/* select and initial edge from the tree */
        j = 0;
		while ( j < phyloTreePtr->numberNodes ) {
			/*if node is not the root */
			if ( phyloTreePtr->nodesPtrArray[ j ]->id != phyloTreePtr->startingNodePtr->id &&
                phyloTreePtr->nodesPtrArray[ j ]->id != nodePtr->id && 
                phyloTreePtr->nodesPtrArray[ j ]->id != internalNodePtr->id) { 

                if ( phyloTreePtr->nodesPtrArray[ j ]->avaliable == FALSE ) {
                    startNodePtr = 	phyloTreePtr->nodesPtrArray[ j ];
                    endNodePtr = startNodePtr->ancestorPtr;
                    j++;
                    break;
                }
            }
            j++;
		}

		/* link selected leaf into the initial edge by using the internal node */
        linkNodeIntoEdge( phyloTreePtr, nodePtr, internalNodePtr, startNodePtr, endNodePtr );

		/* select the edge that leads to the first best score */
        score = labelOptimizeTree( 
            phyloTreePtr, paramsPtr, multiple, previousTreePtr, iteration );//--from iterate_tree.c

        while ( j < phyloTreePtr->numberNodes ) {
            /*if node is not the root */
            if ( phyloTreePtr->nodesPtrArray[ j ]->id != phyloTreePtr->startingNodePtr->id &&
                phyloTreePtr->nodesPtrArray[ j ]->id != nodePtr->id && 
                phyloTreePtr->nodesPtrArray[ j ]->id != internalNodePtr->id ){

                /* select a temporal edge */
                if ( phyloTreePtr->nodesPtrArray[j]->avaliable == FALSE ){
                    tmpStartNodePtr =  phyloTreePtr->nodesPtrArray[ j ];
                    tmpEndNodePtr = tmpStartNodePtr->ancestorPtr;

                    /* unlink the selected node from previous initial edge */
                    unlinkNodeFromEdge( phyloTreePtr, nodePtr, internalNodePtr, startNodePtr, endNodePtr );
                    /* link selected edge into temporal edge of the tree */
                    linkNodeIntoEdge( phyloTreePtr, nodePtr, internalNodePtr, tmpStartNodePtr, tmpEndNodePtr );
                    /* test if this new tree is better than the former */

                    newScore = labelOptimizeTree( 
                        phyloTreePtr, paramsPtr, multiple, previousTreePtr, iteration );//--from iterate_tree.c                    
                    if ( newScore < score ) {
                        if ( DEBUG ) {
                            printf( "[initial tree improved]" );
                            showTreeNewickFormat( phyloTreePtr->startingNodePtr, SHOW_BY_ID );
                        }
                        score = newScore;
                        break;
                    }
                    else{
                        unlinkNodeFromEdge( phyloTreePtr, nodePtr, internalNodePtr, tmpStartNodePtr, tmpEndNodePtr );
                        linkNodeIntoEdge( phyloTreePtr, nodePtr, internalNodePtr, startNodePtr, endNodePtr );
                    }
                }
            }    
            j++;
        }//end while

	}//end for
}

/* 3rd initialization method:
* [new] Create a tree by inserting a "random leaf" into 
* the edge that leads to the "best" tree score */
/*static*/ void createTreeRandomLeaf_BestEdge( TreePtr phyloTreePtr )
{

}

/* 4th initialization method:
* Create a tree by selecting randomly the pair (leaf-edge) from a list 
* that has the best tree scores */
/*static*/ void createTreeBestLeafEdgePair( TreePtr phyloTreePtr )
{

}


static void createTreeWith3LeavesRandomly( TreePtr phyloTreePtr )
{
	TreeNodePtr node1Ptr, node2Ptr, node3Ptr, internalNodePtr;

    node1Ptr = NULL; /* init var */
    node2Ptr = NULL; /* init var */
    node3Ptr = NULL; /* init var */
    internalNodePtr = NULL;	/* init var */

	/* select 3 tree leaf nodes randomly and one internal node */
	selectLeafNodeRandomly(phyloTreePtr, &node1Ptr);
	selectLeafNodeRandomly(phyloTreePtr, &node2Ptr);
	selectLeafNodeRandomly(phyloTreePtr, &node3Ptr);
	selectInternalNode(phyloTreePtr, &internalNodePtr);

	/* startingNode of the tree points to node 1 */ 
	phyloTreePtr->startingNodePtr = node1Ptr;
	
	/* let the node 1 be a leaf with just a right descendant: internalNode*/
	node1Ptr->ancestorPtr = NULL;
	node1Ptr->leftDescPtr = NULL; 
	node1Ptr->rightDescPtr = internalNodePtr;
	internalNodePtr->ancestorPtr = node1Ptr;
	/* let internalNode be the ancestor of node 2 and 3 */
	internalNodePtr->leftDescPtr = node2Ptr;
	internalNodePtr->rightDescPtr = node3Ptr;
	node2Ptr->ancestorPtr = internalNodePtr;
	node3Ptr->ancestorPtr = internalNodePtr;
	/* node2 and node 3 are leaves so descendants are NULL*/
	node2Ptr->leftDescPtr = NULL;
	node2Ptr->rightDescPtr = NULL;
	node3Ptr->leftDescPtr = NULL;
	node3Ptr->rightDescPtr = NULL;
}


static void selectLeafNodeByName( 
    TreePtr phyloTreePtr, TreeNodePtr *nodePtr, char * name )
{
    int i;
    for ( i = 0; i < phyloTreePtr->numberLeaves; i++ ) {
        if ( strcmp( phyloTreePtr->nodesPtrArray[ i ]->organism, name ) == 0 ) {
            phyloTreePtr->nodesPtrArray[ i ]->avaliable = FALSE;
            phyloTreePtr->avaliableLeavesNodes--;
            ( *nodePtr ) = phyloTreePtr->nodesPtrArray[ i ];
            return;
        }
    }

    fprintf(stderr, " stderr: %s leaf not found \n", name); 
    exit(EXIT_FAILURE);    
}

static void selectLeafNodeRandomly( TreePtr phyloTreePtr, TreeNodePtr *nodePtr )
{
	int i, count, low, hi;

	/* random number in the interval  1<=i<=avaliableLeaves */
	low = 1;
	hi = phyloTreePtr->avaliableLeavesNodes;
	count = low + irand( hi ); /* irand from random.c */
	//count = 1 + (rand() % phyloTreePtr->avaliableLeavesNodes); 

	for ( i = 0; i < phyloTreePtr->numberLeaves; i++ ) {
		if ( phyloTreePtr->nodesPtrArray[ i ]->avaliable == TRUE ) {
			count--;
			if ( count == 0 ) {
				phyloTreePtr->nodesPtrArray[ i ]->avaliable = FALSE;
				phyloTreePtr->avaliableLeavesNodes--;
				( *nodePtr ) = phyloTreePtr->nodesPtrArray[ i ];
				break;
			}
		}
	}	
}
 
static void selectInternalNode( TreePtr phyloTreePtr, TreeNodePtr *nodePtr )
{
	int i;
	for(i=phyloTreePtr->numberLeaves; i<phyloTreePtr->numberNodes; i++){
		if (phyloTreePtr->nodesPtrArray[ i ]->avaliable == TRUE){
			phyloTreePtr->nodesPtrArray[ i ]->avaliable = FALSE;
			phyloTreePtr->avaliableInternalNodes--;
			(*nodePtr) = phyloTreePtr->nodesPtrArray[ i ];
			break;
		}
	}
}

/* IMPORTANT NOTE: This function is used just for the Small-Phylogny case */
static void selectInternalNodeByIndex( 
                        TreePtr phyloTreePtr, TreeNodePtr *nodePtr, int index )
{
    if (phyloTreePtr->nodesPtrArray[ index ]->avaliable == TRUE) {
        phyloTreePtr->avaliableInternalNodes--;
        ( *nodePtr ) = phyloTreePtr->nodesPtrArray[ index ];
    }
    else {
        fprintf( stderr, " stderr: not avaliable internal node by index\n" ); 
        exit( EXIT_FAILURE ); 
    }
}

/* IMPORTANT NOTE: This function is used just for the Small-Phylogny case */
static TreeNodePtr getLeafNodePointerByName( TreePtr phyloTreePtr, char *name )
{
    int i;
    TreeNodePtr nodePtr = NULL;

    for ( i = 0; i < phyloTreePtr->numberLeaves; i++ ) {
        if ( strcmp( phyloTreePtr->nodesPtrArray[ i ]->organism, name ) == 0 ) {
            nodePtr = phyloTreePtr->nodesPtrArray[ i ];
            break;
        }
    }

    if ( nodePtr == NULL ) {
        fprintf( stderr, "stderr: pointer to leaf %s not found\n", name ); 
        exit( EXIT_FAILURE );
    }

    return nodePtr;
}

/* IMPORTANT NOTE: This function is used just for the Small-Phylogny case */
static TreeNodePtr getInternalNodePointerByIndex(
                                            TreePtr phyloTreePtr, int index ) 
{
    return phyloTreePtr->nodesPtrArray[ index ];
}

static void selectEdgeRandomly( TreePtr phyloTreePtr, 
	TreeNodePtr *startNodePtr, TreeNodePtr *endNodePtr )
{
	/* NOTE: 
	* We assumed that all nodes except the root have and egde
	* with its ancestor, so when we want to choose an arbitrary edge
	* we simply choose an arbitrary node except the root.
	* The edge is formed by the chosen node and its ancestor */
	
	int i, counter, nodesInTheTree;

	/* choose randomly an edge (from nodes used in the tree) to be modified */
	nodesInTheTree = phyloTreePtr->numberNodes - 
		(phyloTreePtr->avaliableLeavesNodes + phyloTreePtr->avaliableInternalNodes);
	nodesInTheTree--; /* decrement one because of the root */
	counter = 1 + irand(nodesInTheTree);/* 1<=i<=nodesInTheTree */
	//counter = 1 + (rand() % nodesInTheTree); /* 1<=i<=nodesInTheTree */

	for (i=0;i<phyloTreePtr->numberNodes;i++){
		/*if node is the root, then discard this case */
		if (phyloTreePtr->nodesPtrArray[i]->id == phyloTreePtr->startingNodePtr->id)
			continue;
		/*if node is in the tree, decrement the counter*/
		if (phyloTreePtr->nodesPtrArray[i]->avaliable == FALSE){
			counter--;
			/*if counter is zero, then we reached the desired node*/
			if (counter == 0){ 
				(*startNodePtr) = phyloTreePtr->nodesPtrArray[i];
				(*endNodePtr) = (*startNodePtr)->ancestorPtr; /* endNodePtr is ancestor of startNodePtr */
				break;
			}
		}
	}
}

/* add a new leaf node into the edge (startNode, endNode) previously selected by creating 
* a new edge (new leaf node, new internal node) */
static void addRandomLeafNodeIntoEdge( TreePtr phyloTreePtr, 
	TreeNodePtr startNodePtr, TreeNodePtr endNodePtr )
{
	int direction;					/* direction of descendant (left or right) */
	//TreeNodePtr startNodePtr;		/* (param) start node of an edge */
	//TreeNodePtr endNodePtr;		/* (param) end node of an edge */
	TreeNodePtr leafNodePtr; 		/* new leaf to be added */
	TreeNodePtr internalNodePtr; 	/* new internal node to be added */
	TreeNodePtr testNodePtr;

	/* NOTE: 
	* the endNode is the ancestor of the startNode */

    leafNodePtr     = NULL; /* init var */
    internalNodePtr = NULL; /* init var */
    testNodePtr     = NULL; /* init var */ 

	/* select a new random leaf node and a new internal node  */
	selectLeafNodeRandomly(phyloTreePtr, &leafNodePtr);
	selectInternalNode(phyloTreePtr, &internalNodePtr);

	/* determine if startNode is the right or left descendant of endNode */
	testNodePtr = endNodePtr->leftDescPtr;
	if (testNodePtr == NULL || testNodePtr->id != startNodePtr->id){
		direction = RIGHT_DESC; /* startNode is right descendant of endNode */
	}
	else{
		direction = LEFT_DESC; /* startNode is left descendant of endNode */
	}

	/* set internalNode as descendant (left or right) of endNode */
	if (direction == RIGHT_DESC){
		endNodePtr->rightDescPtr = internalNodePtr;
	}
	else{ /* direction is LEFT */
		endNodePtr->leftDescPtr = internalNodePtr;
	}
	/* let endNode be the ancestor of internalNode */
	internalNodePtr->ancestorPtr = endNodePtr; 

	/* set leafNode and startNode as descendants of internalNode */
	if (direction == RIGHT_DESC){
		internalNodePtr->leftDescPtr = leafNodePtr; 
		internalNodePtr->rightDescPtr = startNodePtr;
	}
	else{ /* direction is LEFT */
		internalNodePtr->leftDescPtr = startNodePtr; 
		internalNodePtr->rightDescPtr = leafNodePtr; 
	}
	/* let internalNode be the ancestor of leafNode */
	leafNodePtr->ancestorPtr = internalNodePtr;
	leafNodePtr->leftDescPtr = NULL;
	leafNodePtr->rightDescPtr = NULL;
	/* let internalNode be the ancestor of startNode */
	startNodePtr->ancestorPtr = internalNodePtr;
}

static void selectAvaliableLeaf( TreePtr phyloTreePtr, TreeNodePtr *nodePtr )
{
    int i;
    for (i=0; i<phyloTreePtr->numberLeaves; i++){
        /* select an available leaf */
        if (phyloTreePtr->nodesPtrArray[i]->avaliable == TRUE){
            (*nodePtr) = phyloTreePtr->nodesPtrArray[i];
            (*nodePtr)->avaliable = FALSE;
            phyloTreePtr->avaliableLeavesNodes--;
            break;
        }
    }
}

static void linkNodeIntoEdge( TreePtr phyloTreePtr, TreeNodePtr nodePtr, 
    TreeNodePtr internalNodePtr, TreeNodePtr startNodePtr, TreeNodePtr endNodePtr )
{
    /* link selected into edge of the tree */
    if (endNodePtr->leftDescPtr != NULL && endNodePtr->leftDescPtr->id == startNodePtr->id){
        endNodePtr->leftDescPtr = internalNodePtr;
        internalNodePtr->ancestorPtr = endNodePtr;

        internalNodePtr->leftDescPtr = startNodePtr;
        startNodePtr->ancestorPtr = internalNodePtr;

        internalNodePtr->rightDescPtr = nodePtr;
        nodePtr->ancestorPtr = internalNodePtr;
    }
    else{//endNodePtr->rightDescPtr->id == startNodePtr->id
        endNodePtr->rightDescPtr = internalNodePtr;
        internalNodePtr->ancestorPtr = endNodePtr;

        internalNodePtr->rightDescPtr = startNodePtr;
        startNodePtr->ancestorPtr = internalNodePtr;

        internalNodePtr->leftDescPtr = nodePtr;
        nodePtr->ancestorPtr = internalNodePtr;
    }
}

static void unlinkNodeFromEdge( TreePtr phyloTreePtr, TreeNodePtr nodePtr, 
    TreeNodePtr internalNodePtr, TreeNodePtr startNodePtr, TreeNodePtr endNodePtr )
{
    if (endNodePtr->leftDescPtr != NULL && endNodePtr->leftDescPtr->id == internalNodePtr->id){
        endNodePtr->leftDescPtr = startNodePtr;
        startNodePtr->ancestorPtr = endNodePtr;
    }
    else{//endNodePtr->rightDescPtr->id == internalNodePtr->id
        endNodePtr->rightDescPtr = startNodePtr;
        startNodePtr->ancestorPtr = endNodePtr;
    }

    internalNodePtr->ancestorPtr = NULL;
    internalNodePtr->leftDescPtr = NULL;
    internalNodePtr->rightDescPtr = NULL;    
    nodePtr->ancestorPtr = NULL;
}

/*
NOTE: nodePtr can not be the root node */
void unlinkNode( TreeNodePtr nodePtr, TreeNodePtr *internalNodePtr )
{
    (*internalNodePtr) = nodePtr->ancestorPtr;
 
    if ( (*internalNodePtr)->leftDescPtr != NULL && 
            (*internalNodePtr)->leftDescPtr->id == nodePtr->id ){
        (*internalNodePtr)->leftDescPtr = NULL;
        nodePtr->ancestorPtr = NULL;
    }
    else if ( (*internalNodePtr)->rightDescPtr != NULL && 
                (*internalNodePtr)->rightDescPtr->id == nodePtr->id ){
        (*internalNodePtr)->rightDescPtr = NULL;
        nodePtr->ancestorPtr = NULL;
    }
    else{
        fprintf(stderr, " stderr: internal node is not ancestor of leaf node\n");
        exit(EXIT_FAILURE);
    }    
}

void relinkNode( TreeNodePtr nodePtr, TreeNodePtr internalNodePtr )
{
    if (internalNodePtr->leftDescPtr == NULL){
        internalNodePtr->leftDescPtr = nodePtr;
        nodePtr->ancestorPtr = internalNodePtr;
    }
    else if (internalNodePtr->rightDescPtr == NULL){
        internalNodePtr->rightDescPtr = nodePtr;
        nodePtr->ancestorPtr = internalNodePtr;
    }
    else{
        fprintf(stderr, " stderr: internal node has no NULL descendant\n");
        exit(EXIT_FAILURE);
    }
}

/* NOTE: in the first call of this function, the first parameter must be 
    the starting node of the phylogenetic tree */
void showTreeNewickFormat( TreeNodePtr currentNodePtr, int option )
{
	if (currentNodePtr->ancestorPtr == NULL) /* if currentNode is starting node */
	{ 
		
		printf("-->Newick Format\t: ");
		if (option == SHOW_BY_ID){ 
			printf("(%d,", currentNodePtr->id);
		}
		else if (option == SHOW_BY_NAME){
			printf("(%s,", currentNodePtr->organism);
		}
		else{
        	fprintf(stderr, " stderr: incorrect option for showing tree.\n");
        	exit(EXIT_FAILURE);
    	}
		showTreeNewickFormat(currentNodePtr->rightDescPtr, option);
		printf(")\n");
	}
	else if (currentNodePtr->leftDescPtr == NULL && 
		currentNodePtr->rightDescPtr == NULL) /* if currentNode is a LEAF node*/
	{
		if (option == SHOW_BY_ID){ 
			printf("%d", currentNodePtr->id);
		}
		else if (option == SHOW_BY_NAME){
			printf("%s", currentNodePtr->organism);
		}
		else{
        	fprintf(stderr, " stderr: incorrect option for showing tree.\n");
        	exit(EXIT_FAILURE);
    	}
	}
	else{ /* currentNode is INTERNAL node */
		printf("(");
		showTreeNewickFormat(currentNodePtr->leftDescPtr, option);
		printf(",");
		showTreeNewickFormat(currentNodePtr->rightDescPtr, option);
		printf(")");
	}

}

/* NOTE: in the first call of this function, the second parameter must be 
    the starting node of the phylogenetic tree */
void writeNewickFormatToFile( char *filename, TreeNodePtr currentNodePtr, int option )
{
    FILE *filePtr;  

    if ( ( filePtr = fopen( filename, "w" ) ) == NULL ) {
        fprintf( stderr, " stderr: %s file could not be opened\n", filename );
        exit( EXIT_FAILURE );
    }
    else {
        writeRecursiveNewickFormat( filePtr, filename, currentNodePtr, option );    
    }

    fclose(filePtr);
}

static void writeRecursiveNewickFormat( FILE *filePtr, 
                char *filename, TreeNodePtr currentNodePtr, int option )
{
    if ( currentNodePtr->ancestorPtr == NULL ) /* if currentNode is starting node */
    { 
        if ( option == SHOW_BY_ID) { 
            fprintf( filePtr, "(%d,", currentNodePtr->id );
        }
        else if (option == SHOW_BY_NAME) {
            fprintf( filePtr, "(%s,", currentNodePtr->organism );
        }
        else{
            fprintf( stderr, " stderr: incorrect option for showing tree.\n" );
            exit( EXIT_FAILURE );
        }
        writeRecursiveNewickFormat( filePtr, filename, currentNodePtr->rightDescPtr, option );
        fprintf( filePtr, ")\n" );
    }
    else if ( currentNodePtr->leftDescPtr == NULL && 
        currentNodePtr->rightDescPtr == NULL ) /* if currentNode is a LEAF node*/
    {
        if ( option == SHOW_BY_ID ) { 
            fprintf( filePtr, "%d", currentNodePtr->id );
        }
        else if ( option == SHOW_BY_NAME ) {
            fprintf( filePtr, "%s", currentNodePtr->organism );
        }
        else{
            fprintf( stderr, " stderr: incorrect option for showing tree.\n" );
            exit( EXIT_FAILURE );
        }
    }
    else{ /* currentNode is INTERNAL node */
        fprintf( filePtr, "(" );
        writeRecursiveNewickFormat( filePtr, filename, currentNodePtr->leftDescPtr, option );
        fprintf( filePtr, "," );
        writeRecursiveNewickFormat( filePtr, filename, currentNodePtr->rightDescPtr, option );
        fprintf( filePtr, ")" );
    }
}

void showNode( TreePtr phyloTreePtr, TreeNodePtr nodePtr)
{
	int i = nodePtr->index;

    printf( "----------------\n" );
	printf( "Pointer value: %p\n", phyloTreePtr->nodesPtrArray[ i ] );
	printf( "Node id\t\t: %d\n", phyloTreePtr->nodesPtrArray[ i ]->id );
			
	if ( phyloTreePtr->nodesPtrArray[ i ]->type == LEAF_NODE ) {
	    printf( "Organism\t: %s\n", phyloTreePtr->nodesPtrArray[ i ]->organism );
	    printf( "Type\t\t: %s\n","LEAF NODE" );
	}
	else{
	   printf( "Type\t\t: %s\n","INTERNAL NODE" );
	}
			
	if ( phyloTreePtr->nodesPtrArray[ i ]->ancestorPtr == NULL )
	    printf("Ancestor Id\t: %s\n", "NULL");
	else
        printf("Ancestor Id\t: %d (Pointer: %p)\n", 
				phyloTreePtr->nodesPtrArray[ i ]->ancestorPtr->id,
				phyloTreePtr->nodesPtrArray[ i ]->ancestorPtr);

    if ( phyloTreePtr->nodesPtrArray[ i ]->leftDescPtr == NULL )
		printf( "Left Desc. Id\t: %s\n", "NULL");
	else
		printf( "Left Desc. Id\t: %d (Pointer: %p)\n", 
				phyloTreePtr->nodesPtrArray[ i ]->leftDescPtr->id,
				phyloTreePtr->nodesPtrArray[ i ]->leftDescPtr );
			
	if ( phyloTreePtr->nodesPtrArray[ i ]->rightDescPtr == NULL )
		printf( "Right Desc. Id\t: %s\n", "NULL" );
	else
		printf( "Right Desc. Id\t: %d (Pointer: %p)\n", 
				phyloTreePtr->nodesPtrArray[ i ]->rightDescPtr->id,
				phyloTreePtr->nodesPtrArray[ i ]->rightDescPtr );
	printf("----------------\n");		
}

void showNodesArray( TreePtr phyloTreePtr )
{
	int i;
	for( i = 0; i < phyloTreePtr->numberNodes; i++ ) {
        //if (phyloTreePtr->nodesPtrArray[i]->avaliable == TRUE)
        //    continue;
		showNode( phyloTreePtr, phyloTreePtr->nodesPtrArray[ i ] );
	}
}

void showGenomes( TreePtr phyloTreePtr, ParametersPtr paramsPtr ) 
{
    int i, j;
    printf( "Number of genomes: %d\n", phyloTreePtr->numberLeaves );
    printf( "Number of genes: %d\n", phyloTreePtr->numberGenes );

    if ( paramsPtr->distanceType == INVERSION_DIST ) {
        for( i = 0; i < phyloTreePtr->numberNodes; i++ ) {
            for( j = 0; j < phyloTreePtr->numberGenes; j++ ) {
                printf( "%d, ", phyloTreePtr->nodesPtrArray[ i ]->genome[ j ] );
            }  
            printf("\n"); 
        }
        printf("----------------\n");
    }
    else if ( paramsPtr->distanceType == DCJ_DIST ) {
        for( i = 0; i < phyloTreePtr->numberNodes; i++ ) {
            printf( "%s: ",  phyloTreePtr->nodesPtrArray[ i ]->organism );
            for( j = 0; j < phyloTreePtr->nodesPtrArray[ i ]->numPointsDCJ; j++ ) {
                if ( phyloTreePtr->nodesPtrArray[ i ]->genomeDCJ[ j ]->type == TELOMERE ) {
                    printf( "{%d}, ", phyloTreePtr->nodesPtrArray[ i ]->genomeDCJ[ j ]->x );
                }
                else { // phyloTreePtr->nodesPtrArray[ i ]->genomeDCJ[ j ]->type == ADJACENCY
                    printf( "{%d,%d}, ", phyloTreePtr->nodesPtrArray[ i ]->genomeDCJ[ j ]->x,
                                    phyloTreePtr->nodesPtrArray[ i ]->genomeDCJ[ j ]->y );
                }
            }  
            printf("\n\n"); 
        }
        printf("----------------\n");
    } 
}

/* show the multiple candidates of the leaves that hold: ( num cand > 0) */
void showMultipleLeaves( 
    TreePtr phyloTreePtr, MultipleLeafPtr *multiple, ParametersPtr paramsPtr )
{
    int i, j, k;

    if ( paramsPtr->distanceType == DCJ_DIST ) {
        for( i = 0; i < phyloTreePtr->numberLeaves; i++ ) {
            if ( multiple[ i ]->numLeafCandidates > 0 ) {

                printf( "%s: \n",  phyloTreePtr->nodesPtrArray[ i ]->organism );

                for( j = 0; j < multiple[ i ]->numLeafCandidates; j++ ) {
                    printf("Candidate %d: ", j );
                    for( k = 0; k < multiple[ i ]->candidates[ j ]->numPointsDCJ; k++ ) {

                        if ( multiple[ i ]->candidates[ j ]->genomeDCJ[ k ]->type == TELOMERE ) {
                            printf( "{%d}, ", multiple[ i ]->candidates[ j ]->genomeDCJ[ k ]->x );
                        }
                        else { // genomeDCJ[ j ]->type == ADJACENCY
                            printf( "{%d,%d}, ", multiple[ i ]->candidates[ j ]->genomeDCJ[ k ]->x,
                                            multiple[ i ]->candidates[ j ]->genomeDCJ[ k ]->y );
                        }
                    }
                    printf("\n");
                }  
                printf("\n\n");
            }
        }//end-for   
    }
    else {
        fprintf( stderr, "Multiple leaves are just allowed for DCJ distance.\n" );
    }
}

/* ---------------------------------------------------------------------- */
/* IMPORTANT NOTE: This function is used just for the Small-Phylogny case */
void showResultsSmallPhylogeny( TreePtr phyloTreePtr, 
                        enum distances dist, int score, double timediff )
{
    printf("-----------------------------\n");
    printf("Program for Small Phylogeny\n");
    printf("-----------------------------\n");
    printf("Number of genomes\t: %d\n", phyloTreePtr->numberLeaves);
    printf("Number of genes\t\t: %d\n", phyloTreePtr->numberGenes);
    if ( dist == INVERSION_DIST )
        printf("Reversal Score\t\t: %d\n", score);
    else if ( dist == DCJ_DIST )
        printf("DCJ Score\t\t: %d\n", score);
    else
        printf("Breakpoint Score\t: %d\n", score);
    printf("Total time\t\t: %.2f seconds\n", timediff);
    showTreeNewickFormat(phyloTreePtr->startingNodePtr, SHOW_BY_ID);//--from tree.c
    showTreeNewickFormat(phyloTreePtr->startingNodePtr, SHOW_BY_NAME);//--from tree.c
    //showNodesArray(phyloTreePtr);//--method from tree.c
}