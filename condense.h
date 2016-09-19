/*
 ============================================================================
 Authors :
 	Jose Luis Soncco-Alvarez and Mauricio Ayala-Rincon
 	Group of Theory of Computation
 	Universidade de Brasilia (UnB) - Brazil
 ============================================================================
 */

#ifndef CONDENSE_H_
#define CONDENSE_H_

#include "my_structs.h"

#define INCREASING 0
#define DECREASING 1
#define SPLIT 0 	/* split between two chrosomomes represents a $ or @ */
#define LINEAR_SYM 0 	/* $ */
#define CIRCULAR_SYM 1 	/* @ */

struct rawGenome {
	int *genome;
	int numberElements; /* num of elements of genome: number of genes + number of symbols($ and @) */
	int *chromosomeType; /* array of linear($) or circular(@) symbols */
	int numberChromosomes;
	char *organism;
}; 
typedef struct rawGenome RawGenome;
typedef RawGenome *RawGenomePtr;

struct rawDataset {
	int numberGenomes;
	int numberGenes;
	int *numberChromosomesArray;
	RawGenomePtr *rgenomes;	
};
typedef struct rawDataset RawDataset;
typedef RawDataset *RawDatasetPtr;

struct condenseKey {
	int gene;		/* the lower extreme (start or end) */
	int start;		/* start of interval */
	int end;		/* end of interval */
	int diff;		/* difference of extremes of interval */
	int orientation; /* orientation of the sequence: increasing or decreasing */
};
typedef struct condenseKey CondenseKey;
typedef CondenseKey *CondenseKeyPtr;

/* set of keys for condensation */
struct setKeys {
	CondenseKeyPtr *condkeyPtrArray;
	int numberKeys;
};
typedef struct setKeys SetKeys;
typedef SetKeys *SetKeysPtr;

void allocateMemoryForRawData( RawDatasetPtr rdatasetPtr);
void allocateMemoryForRawGenomes( RawDatasetPtr rdatasetPtr );
void allocateMemoryForKeys( int numberGenes, SetKeysPtr setkeysPtr );
void freeKeys( int numberGenes, SetKeysPtr setkeysPtr );
void freeRawDataset( RawDatasetPtr rdatasetPtr );
void findSetKeysForCondensing( RawDatasetPtr rdatasetPtr, SetKeysPtr setkeysPtr );
void condenseGenomes( RawDatasetPtr rdatasetPtr, SetKeysPtr setkeysPtr );
void unpackCondensedGenomes( RawDatasetPtr rdatasetPtr, SetKeysPtr setkeysPtr );
void showKeys( SetKeysPtr setkeysPtr );

#endif