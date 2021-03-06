/*
 ============================================================================
 Project    : Greedy Algorithm for Small Phylogeny Problem
 File       : main.c
 Authors    : Jose Luis Soncco-Alvarez and Mauricio Ayala-Rincon
                Group of Theory of Computation
                Universidade de Brasilia (UnB) - Brazil
 ============================================================================
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <limits.h>

#include "my_structs.h"
#include "condense.h"
#include "tree.h"
#include "iterate_tree.h"
#include "measure_time.h"
#include "random.h"

#include "dcjdist.h"

static void initParameters( ParametersPtr paramsPtr );
static void readCommandLine( int argc, char *argv[], ParametersPtr paramsPtr );
static int readNumberGenomes( char *filename );
static void readNumberGenesAndChromosomes( char *filename, RawDatasetPtr rdatasetPtr );
static void readRawData( char *filename, RawDatasetPtr rdatasetPtr );

int main( int argc, char **argv )
{
    struct timeval t_ini, t_fin;
    gettimeofday( &t_ini, NULL );//---------------------------take start time--
    
    Tree            phyloTree, minTree;
    RawDataset      rdataset; 
    SetKeys         setkeys;  /* set of keys for condensing genomes */
    int             minScore, score, i;
    Parameters      params;

    MultipleLeafPtr *multiple; /*in case of multiple genomes per leaf*/

    initParameters( &params );
    readCommandLine( argc, argv, &params );
    srand( params.seed ); if ( DEBUG ) { printf( "\nseed: %d\n", params.seed ); }

    /* read raw data */
    rdataset.numberGenomes = readNumberGenomes( params.testsetName );
    allocateMemoryForRawData( &rdataset );//--from condense.c
    readNumberGenesAndChromosomes( params.testsetName, &rdataset );
    allocateMemoryForRawGenomes( &rdataset );//--from condense.c
    readRawData( params.testsetName, &rdataset );

    /* condense raw data */
    allocateMemoryForKeys( rdataset.numberGenes, &setkeys );//--from condense.c
    findSetKeysForCondensing( &rdataset, &setkeys );//--from condense.c
    condenseGenomes( &rdataset, &setkeys );//--from condense.c

    minScore = INT_MAX;
    for ( i = 0; i < params.iterations; i++ ) {
        /* read genomes from raw data into phylogenetic tree */
        readNumberLeavesAndGenes( &phyloTree, &params, &rdataset );//--from tree.c     
        allocateMemoryForNodes( &phyloTree, &params );//--from tree.c    
        readGenomesFromRawData( 
            &phyloTree, &params, &rdataset, &multiple );//--from tree.c

        /* create topology and optimize tree */
        createTopologyFromNewickFormat( &phyloTree, &params );//from tree.c 
        score = labelOptimizeTree( 
                    &phyloTree, &params, multiple, &minTree, i );//--iterate tree.c
        if ( DEBUG ) { printf( "iteration %d result : %d\n", i, score ); }
        if ( score < minScore ) { 
            minScore = score;
            /* copy phyloTree into a minTree */
            minTree.numberLeaves = phyloTree.numberLeaves;
            minTree.numberGenes = phyloTree.numberGenes;
            allocateMemoryForNodes( &minTree, &params );//--from tree.c
            copyTreeInto( &minTree, &phyloTree, TRUE, &params);//--from tree.c
        }

        freeMultipleLeafs( &phyloTree, &multiple, &params );//--from tree.c
        freeTree( &phyloTree, &params );//--from tree.c
    }

    gettimeofday( &t_fin, NULL );//---------------------------take final time--
    double timediff = timeval_diff( &t_fin, &t_ini );//--from measure_time.h   

    /* show results */
    if ( SHOW_JUST_SCORE == TRUE ) 
        printf( "%d %.2f\n", minScore, timediff );
    else 
        showResultsSmallPhylogeny( &minTree, 
            params.distanceType, minScore, timediff );//--from tree.c

    /* free memory */
    freeTree( &minTree, &params );//--from tree.c
    freeKeys( rdataset.numberGenes, &setkeys );//--from condense.c
    freeRawDataset( &rdataset );//--from condense.c

	return 0;
}

static void initParameters( ParametersPtr paramsPtr )
{
    paramsPtr->problem          = SMALL_PHYLOGENY;
    paramsPtr->seed             = time( NULL);
    paramsPtr->testsetName      = "";
    paramsPtr->newickFile       = "";
    paramsPtr->circular         = TRUE; //This value is set in function readGenomes() at tree.c
    paramsPtr->distanceType     = INVERSION_DIST;
    paramsPtr->solver           = CAPRARA_INV_MEDIAN;
    paramsPtr->useOutgroup      = FALSE;
    paramsPtr->outgroup         = "";
    paramsPtr->preferredType    = ANY;
    paramsPtr->initMethod       = R_LEAF_1BEST_EDGE;
    paramsPtr->opt              = BLANCHETTE; // (*)
    paramsPtr->useMultipleGenomesOneLeaf = FALSE;
    paramsPtr->iterations       = 1;
}

static void readCommandLine( int argc, char *argv[], ParametersPtr paramsPtr )
{   
    int i;
    char option;

    /* show parameter options */
    if ( argc == 1 ) {
        fprintf( stdout, "\nParameter Options:\n" );
        fprintf( stdout, "\t-d : evolutionary distance \n");
        fprintf( stdout, "\t\t -d rev : reversal\n\t\t -d dcj : double-cut-join\n");
        fprintf( stdout, "\t-f : dataset filename\n" );
        fprintf( stdout, "\t\t -f filename\n" );
        fprintf( stdout, "\t-k : topology in Newick format\n" );
        fprintf( stdout, "\t\t -k filename\n" );
        fprintf( stdout, "\t-s : [optional] seed \n" );
        fprintf( stdout, "\t\t -s some_seed\n" );
        fprintf( stdout, "\t\t(seed taken from system time by default if option is omitted)\n\n" );
        fprintf( stdout, "\t*** The following options work just for the DCJ distance (-d dcj) : \n\n");
        fprintf( stdout, "\t-i : [optional] number of iterations\n" );
        fprintf( stdout, "\t\t -i number\n" );
        fprintf( stdout, "\t\t(one iteration is used by default if option is omitted)\n" );
        fprintf( stdout, "\t-o : [optional] optimization method for DCJ\n" );
        fprintf( stdout, "\t\t -o gre : Greedy Candidates opt.\n" );
        fprintf( stdout, "\t\t -o kov : Kovac opt.\n" );
        fprintf( stdout, "\t\t(Kovac opt is used by default if option is omitted)\n" );
        fprintf( stdout, "\t-r : [optional] preferred dcj genome structure \n" );
        fprintf( stdout, "\t\t -r 0 : any genome structure\n" );
        fprintf( stdout, "\t\t -r 1 : one circular chromosome\n" );
        fprintf( stdout, "\t\t -r 2 : one or more linear chromosomes\n" );    
        fprintf( stdout, "\t\t -r 3 : one circular, or one or more linear chromosomes\n" );
        fprintf( stdout, "\t\t( -r 0 is used by default if option is omitted)\n\n" );
        exit( EXIT_FAILURE );
    }

    /* read parameters from command line */
    i = 1;
    while ( i < argc) {
        if ( argv[ i ][ 0 ] == '-' && (i + 1) < argc ) { 
            option = argv[ i ][ 1 ]; 
            switch ( option ) {
                case 'd':
                    if ( strcmp( argv[ i + 1 ], "rev" ) == 0 ) {
                        paramsPtr->distanceType = INVERSION_DIST;
                        paramsPtr->opt = BLANCHETTE;
                    }
                    else if ( strcmp( argv[ i + 1 ], "dcj" ) == 0 ) {
                        paramsPtr->distanceType = DCJ_DIST;
                        paramsPtr->opt = KOVAC;
                    }
                    else {
                        fprintf( stderr, 
                            " stderr: incorrect distance (-d).\n" );
                        exit( EXIT_FAILURE );
                    }
                    break;
                case 'f':
                    paramsPtr->testsetName = argv[ i + 1 ]; 
                    break;
                case 'k':
                    paramsPtr->newickFile = argv[ i + 1 ]; 
                    break;
                case 'i':
                    paramsPtr->iterations = atoi( argv[ i + 1 ] );
                    break;
                case 'o':
                    if ( strcmp( argv[ i + 1 ], "gre" ) == 0 ) {
                        paramsPtr->opt = GREEDY_CANDIDATES;
                    }
                    else if ( strcmp( argv[ i + 1 ], "kov" ) == 0 ) {
                        paramsPtr->opt = KOVAC;
                    }
                    else {
                        fprintf( stderr, 
                            " stderr: incorrect opt method (-o).\n" );
                        exit( EXIT_FAILURE ); 
                    }
                    break;
                case 's':
                    paramsPtr->seed = atoi( argv[ i + 1 ] );
                    break;
                case 'r':
                    if ( strcmp( argv[ i + 1 ], "0" ) == 0 ) {
                        paramsPtr->preferredType = ANY;
                    }
                    else if ( strcmp( argv[ i + 1 ], "1" ) == 0 ) {
                        paramsPtr->preferredType = ONE_CIRCULAR_CH;
                    }
                    else if ( strcmp( argv[ i + 1 ], "2" ) == 0 ) {
                        paramsPtr->preferredType = ONE_OR_MORE_LINEAR_CH;
                    }
                    else if ( strcmp( argv[ i + 1 ], "3" ) == 0 ) {
                        paramsPtr->preferredType = 
                            ONE_CIRCULAR_or_ONE_OR_MORE_LINEAR_CH;
                    }
                    else {
                        fprintf( stderr, 
                            " stderr: incorrect penalty type (-r).\n" );
                        exit( EXIT_FAILURE ); 
                    }
                    break;
                default:
                    fprintf( stderr, 
                        " stderr: incorrect option: %c.\n", option );
                    exit( EXIT_FAILURE );
            }
            i = i + 2;  
        }
        else{
            fprintf( stderr, " stderr: incorrect options or parameters.\n" );
            exit( EXIT_FAILURE );
        }
    }

    /* discard some undesired cases */
    if ( paramsPtr->distanceType == INVERSION_DIST && 
            paramsPtr->preferredType > 0 ) {
        fprintf( stderr, 
            " stderr: the program does not support using the reversal" );
        fprintf( stderr, " distance and a preferred genome structure.\n" );
        exit( EXIT_FAILURE );
    }
}

static int readNumberGenomes( char *filename )
{   
    FILE *filePtr;
    int c; /* use int (not char) for the EOF */
    int charCounter = 0;
    int newlineCounter = 0;

    if ( ( filePtr = fopen( filename , "r" ) ) == NULL ) {
        fprintf( stderr, " stderr: %s file could not be opened\n", filename  );
        exit( EXIT_FAILURE );
    }
    else {
        /* count the newline characters */
        while ( ( c = fgetc( filePtr ) ) != EOF ) {
            if ( c == '\n' ) { 
                /* if there exist at least one char before '\n'
                * increment the line counter. */
                if (charCounter > 0) { 
                    newlineCounter++;
                    charCounter=0;
                }
            }
            else { /* c is any char different from '\n' */
                /* dont count white spaces*/
                if ( c != ' ' ) {
                    charCounter++;
                }
            }
        }
        /* if last char before EOF was not '\n' and
        chars different from ' ' were found */
        if (charCounter > 0) { 
            newlineCounter++;
            charCounter = 0;
        }
    }

    fclose( filePtr );

    //printf("newlineCounter: %d\n", newlineCounter); 
    if ( newlineCounter == 0 || newlineCounter % 2 != 0 ) {
        fprintf( stderr, " stderr: data of %s file has incorrect format\n", filename  );
        exit( EXIT_FAILURE );
    }

    if ( newlineCounter / 2 <= 3) {
        fprintf( stderr, " stderr: the number of genomes of %s file must be  greater than 3.\n", filename );
        exit( EXIT_FAILURE );
    }
    return newlineCounter / 2;
}

/* NOTE: before calling this function initialize with zero: numberGenes, 
        and elements of numberChromosomesArray */
static void readNumberGenesAndChromosomes( char *filename, RawDatasetPtr rdatasetPtr )  
{
    FILE *filePtr;
    int c, lastc; /* use int (not char) for the EOF */
    int count, i, k;

    if ( ( filePtr = fopen( filename, "r" ) ) == NULL ) {
        fprintf( stderr, 
            " stderr: %s file could not be opened\n", filename );
        exit( EXIT_FAILURE );
    }
    else {              
        /* read first line and discard */
        while ( ( c = fgetc( filePtr ) ) != EOF ) {
            if ( c == '\n' )
                break;
        }
        /* read second line and count number of genes */
        i = 0; /* digit counter*/
        while ( ( c = fgetc( filePtr ) ) != EOF ) {
            if ( c == ' ' || c == '@' || c == '$' || c == '\n' ) {
                /* if the digit counter "i" has at least one digit 
                * increment the genes counter */
                if ( i > 0 ) {
                    rdatasetPtr->numberGenes++; 
                    //( *numberGenes )++;
                }
                i = 0; /* re-start digit counter */

                if ( c == '@' || c == '$' ) {
                    rdatasetPtr->numberChromosomesArray[ 0 ]++;
                }

                if ( c == '\n' )
                    break;
            }
            else { /* c should be a digit */
                i++;
            }
        }

        /* verify if the other genomes have the same number of genes */
        for ( k = 1; k < rdatasetPtr->numberGenomes; k++ ) {
            count = 0;
            /* read name_of_genome line and discard */
            while ( ( c = fgetc( filePtr ) ) != EOF ) {
                if ( c == '\n' )
                    break;
            }
            /* read next line and count number of genes */
            i = 0;
            c = fgetc( filePtr );
            lastc = EOF;

            while ( c != EOF ) {
                if ( c != ' '  && c != '\n' ) {  
                    lastc = c;
                }
                if ( c == ' ' || c == '@' || c == '$' || c == '\n' ) {
                    if ( i > 0 ) {
                        count++;
                    }
                    i = 0;

                    if ( c == '@' || c == '$' ) {
                        rdatasetPtr->numberChromosomesArray[ k ]++;
                    }

                    if ( c == '\n' )
                        break;
                }
                else { /* c should be a number */
                    i++;
                }
                c = fgetc( filePtr );
            }

            if ( lastc != '@' && lastc != '$' ) {
                fprintf( stderr, 
                    " stderr: there is a genome whose last char is not @ or $ in %s\n", filename ); 
                exit( EXIT_FAILURE );
            }

            if ( rdatasetPtr->numberGenes != count ) {
                fprintf( stderr, 
                    " stderr: number of genes are not the same in %s\n", filename ); 
                exit( EXIT_FAILURE );
            } 
        }
    }

    fclose( filePtr );
}

/* Read "Raw" genomes, that is, genomes before condensation */
static void readRawData( char *filename, RawDatasetPtr rdatasetPtr )
{
    FILE *filePtr;
    int c; /* use int (not char) for the EOF */
    int i, j, k, h;
    char buffer[ MAX_STRING_LEN ];

    if ( ( filePtr = fopen( filename, "r" ) ) == NULL ) {
        fprintf( stderr, " stderr: %s file could not be opened.\n", filename );
        exit( EXIT_FAILURE );
    }
    else {
        /* read genomes from file into leaves */
        for( i = 0; i < rdatasetPtr->numberGenomes; i++ ) {
            /* verify that char symbol '>' exists */
            if ( ( c = fgetc( filePtr ) ) != EOF ) {
                if ( c != '>' ) {
                    fprintf( stderr, " stderr: incorrect format of organism name in %s.\n", filename );
                    exit( EXIT_FAILURE );
                }
            }
            /* read organism name into buffer */
            k = 0;
            while ( ( c = fgetc( filePtr ) ) != EOF ) {
                if ( c == '\n' ) {
                    buffer[ k ] = '\0';
                    break;
                }
                else {
                    buffer[ k ] = c;
                    k++;

                    if ( k + 1 > MAX_STRING_LEN) {
                        fprintf( stderr, " stderr: increment MAX_STRING_LEN value!\n" );
                        exit( EXIT_FAILURE );
                    }
                }
            }

            /* copy buffer into leaf node*/
            rdatasetPtr->rgenomes[ i ]->organism = malloc( ( k + 1 ) * sizeof( char ) );
            strcpy( rdatasetPtr->rgenomes[ i ]->organism, buffer ); 

            /* read genes from next line*/
            k = 0; j = 0; h = 0;
            while ( ( c = fgetc( filePtr ) ) != EOF ) {
                if ( c == ' ' || c == '@' || c == '$' || c == '\n' ){
                    buffer[ k ] = '\0';
                    if ( k > 0 ) {
                        rdatasetPtr->rgenomes[ i ]->genome[ j ] = atoi( buffer );
                        j++;
                        //printf("%d,", atoi(buffer));    
                    }
                    k = 0;

                    if ( c == '@' || c == '$') {
                        rdatasetPtr->rgenomes[ i ]->genome[ j ] = SPLIT;
                        rdatasetPtr->rgenomes[ i ]->chromosomeType[ h ] = 
                                        ( c == '@' ? CIRCULAR_SYM : LINEAR_SYM );
                        j++;
                        h++;
                    }

                    if ( c == '\n' )
                        break;
                }
                else { // c should be a digit
                    buffer[ k ] = c;
                    k++;
                }
            }

        }//end for
    }

    fclose( filePtr );
}

