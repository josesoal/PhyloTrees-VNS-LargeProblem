/*
 ============================================================================
 Authors :
 	Jose Luis Soncco-Alvarez and Mauricio Ayala-Rincon
 	Group of Theory of Computation
 	Universidade de Brasilia (UnB) - Brazil
 ============================================================================
 */

#ifndef LABEL_TREE_H_
#define LABEL_TREE_H_

#include "structs.h"
#include "my_structs.h"

int labelOptimizeTree( TreePtr phyloTreePtr, ParametersPtr paramsPtr, 
	MultipleLeafPtr *multiple, TreePtr previousTreePtr, int iteration );
int scoreTree( TreePtr phyloTreePtr, TreeNodePtr nodePtr, ParametersPtr paramsPtr );
int calculateDistance( int *genome1, int *genome2, int numberGenes, ParametersPtr paramsPtr );

void allocateMemForCandidate( TreePtr phyloTreePtr, CandidatePtr *candPtr );
void freeCandidate( TreePtr phyloTreePtr, CandidatePtr *candPtr );

#endif /* LABEL_TREE_H_ */