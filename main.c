/*
 ============================================================================
 Project    : VNS for the Large Phylogeny Problem
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

#include "my_structs.h"
#include "condense.h"
#include "tree.h"
#include "iterate_tree.h"
#include "vns.h"
#include "measure_time.h"
#include "random.h"

#include "dcjdist.h"

void initParameters( ParametersPtr paramsPtr);
void readCommandLine( int argc, char *argv[], ParametersPtr paramsPtr );
int readNumberGenomes( char *filename );
void readNumberGenesAndChromosomes( char *filename, RawDatasetPtr rdatasetPtr );
void readRawGenomes( char *filename, RawDatasetPtr rdatasetPtr );

//int readNumberGenes( ParametersPtr paramsPtr, int numberGenomes );

int main( int argc, char **argv )
{
    struct timeval t_ini, t_fin;
    gettimeofday( &t_ini, NULL );//---------------------------------take start time--
    
    Tree            phyloTree;
    RawDataset      rdataset;
    SetKeys         setkeys;  /* set of keys for condensing genomes */
    int             score;
    Parameters      params;

    initParameters( &params );
    readCommandLine( argc, argv, &params );
    srand( params.seed ); if ( DEBUG ) { printf( "\nseed: %d\n", params.seed ); }

    /* phyloTree.numberLeaves = readNumberGenomes( &params );
    phyloTree.numberGenes = readNumberGenes( &params, phyloTree.numberLeaves );
    allocateMemoryForNodes( &phyloTree, &params );//--from tree.c
    readGenomesFromFile( &phyloTree, &params );//--from tree.c
    
    if ( params.distanceType == INVERSION_DIST ) {
        allocateMemoryForKeys( &phyloTree, &sKeys );//--from condense.c
        findSetKeysForCondensingLeaves( &phyloTree, &sKeys );//--from condense.c
        condenseLeafNodes( &phyloTree, &sKeys );//--from condense.c
	} */

    /* read raw data */
    rdataset.numberGenomes = readNumberGenomes( params.testsetName );
    allocateMemoryForRawData( &rdataset );//--from condense.c
    readNumberGenesAndChromosomes( params.testsetName, &rdataset );
    allocateMemoryForRawGenomes( &rdataset );//--from condense.c
    readRawGenomes( params.testsetName, &rdataset );

    /* condense raw data */
    allocateMemoryForKeys( rdataset.numberGenes, &setkeys );//--from condense.c
    findSetKeysForCondensing( &rdataset, &setkeys );//--from condense.c
    condenseGenomes( &rdataset, &setkeys );//--from condense.c

    /* read genomes from raw data into phylogenetic tree */
    phyloTree.numberLeaves = rdataset.numberGenomes;
    phyloTree.numberGenes = rdataset.numberGenes; 
    allocateMemoryForNodes( &phyloTree, &params );//--from tree.c
    //CONTINUE HERE... (create a function that read)


    return 0;//delete this

    score = createInitialTreeTopology( &phyloTree, &params ); //--from tree.c
    score = VNS( &phyloTree, &params );//--from vns.c

    /* refine the final tree by using an exhaustive mutation procedures */
    score = exhaustiveSubtreeScramble( &phyloTree, &params, score );
    score = exhaustiveLeafSwap( &phyloTree, &params, score );
    
    /*if ( params.distanceType == INVERSION_DIST ) {
        reverseCondensedNodes( &phyloTree, &setkeys, EACH_NODE );//--from condense.c
    }*/

    gettimeofday( &t_fin, NULL );//---------------------------------take final time--
    double timediff = timeval_diff( &t_fin, &t_ini );//--from measure_time.h
    
    /* show results */
    if ( SHOW_JUST_SCORE == TRUE )
        printf( "%d %.2f\n", score, timediff );
    else
        showResults( &phyloTree, params.distanceType, score, timediff );//--from vns.c

    /* free memory */
    if ( params.distanceType == INVERSION_DIST ) {
        freeKeys( rdataset.numberGenes, &setkeys );//--from condense.c
    }
    freeTree( &phyloTree, &params );//--from tree.c

	return 0;
}

void initParameters( ParametersPtr paramsPtr )
{
    paramsPtr->testsetName      = "";
    paramsPtr->distanceType     = INVERSION_DIST;
    paramsPtr->solver           = CAPRARA_INV_MEDIAN;
    paramsPtr->unichromosomes   = TRUE;
    paramsPtr->seed             = time( NULL);
    paramsPtr->useOutgroup      = FALSE;
    paramsPtr->outgroup         = "";
    paramsPtr->circular         = TRUE; //Not used

    paramsPtr->initMethod       = R_LEAF_1BEST_EDGE;
    paramsPtr->opt              = BLANCHETTE; // (*)
    //opt = KOVAC; //super slow, not good results (worst than BLANCHETTE)
    //opt = KOVAC_MODIFIED; // is slow, "almost" the same results as BLANCHETTE
    //opt = MEDIAN_CANDIDATE; // is slow, the same results as BLANCHETTE (candidates dont improve nothing)
}

void readCommandLine( int argc, char *argv[], ParametersPtr paramsPtr )
{   
    int i;
    char option;

    /* show parameter options */
    if ( argc == 1 ) {
        fprintf( stdout, "Parameter Options:\n" );
        fprintf( stdout, "\t-d : distance (rev, dcj) \n" );
        fprintf( stdout, "\t-f : dataset filename\n" );
        fprintf( stdout, "\t-m : multiple-chromosomes [optional] [default: single-chromosomes]\n" );
        fprintf( stdout, "\t-s : seed [optional] [default: seed from system time]\n" );
        fprintf( stdout, "\t-g : outgroup [optional] [default: not use outgroup]\n" );
        
        //fprintf( stdout, " using the default testset: testsets/camp05_cond\n" );
        //fprintf( stdout, " try other >> ./main -t testsets/camp07_cond\n" );
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
                    }
                    else if ( strcmp( argv[ i + 1 ], "dcj" ) == 0 ) {
                        paramsPtr->distanceType = DCJ_DIST;
                        paramsPtr->opt = KOVAC_MODIFIED;
                    }
                    else {
                        fprintf( stderr, " stderr: incorrect distance (-d).\n" );
                        exit( EXIT_FAILURE );
                    }
                    break;
                case 'f':
                    paramsPtr->testsetName = argv[ i + 1 ]; 
                    break;
                case 'm':
                    paramsPtr->unichromosomes = FALSE;
                    break;
                case 's':
                    paramsPtr->seed = atoi( argv[ i + 1 ] );
                    break;
                case 'g': 
                    paramsPtr->outgroup = argv[ i + 1 ]; 
                    paramsPtr->useOutgroup = TRUE;
                    break;
                default:
                    fprintf( stderr, " stderr: incorrect option: %c.\n", option );
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
            paramsPtr->unichromosomes == FALSE ) {
        fprintf( stderr, " stderr: the program does not support using the reversal distance and multiple-chromosome genomes.\n" );
        exit( EXIT_FAILURE );
    }
}

int readNumberGenomes( char *filename )
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
void readNumberGenesAndChromosomes( char *filename, RawDatasetPtr rdatasetPtr )  
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
void readRawGenomes( char *filename, RawDatasetPtr rdatasetPtr )
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
            rdatasetPtr->rgenomePtrArray[ i ]->organism = malloc( ( k + 1 ) * sizeof( char ) );
            strcpy( rdatasetPtr->rgenomePtrArray[ i ]->organism, buffer ); 

            /* read genes from next line*/
            k = 0; j = 0; h = 0;
            while ( ( c = fgetc( filePtr ) ) != EOF ) {
                if ( c == ' ' || c == '@' || c == '$' || c == '\n' ){
                    buffer[ k ] = '\0';
                    if ( k > 0 ) {
                        rdatasetPtr->rgenomePtrArray[ i ]->genome[ j ] = atoi( buffer );
                        j++;
                        //printf("%d,", atoi(buffer));    
                    }
                    k = 0;

                    if ( c == '@' || c == '$') {
                        rdatasetPtr->rgenomePtrArray[ i ]->genome[ j ] = SPLIT;
                        rdatasetPtr->rgenomePtrArray[ i ]->chromosomeType[ h ] = 
                                        ( c == '@' ? CIRCULAR : LINEAR );
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

