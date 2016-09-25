/*
 ============================================================================
 Authors :
 	Jose Luis Soncco-Alvarez and Mauricio Ayala-Rincon
 	Group of Theory of Computation
 	Universidade de Brasilia (UnB) - Brazil
 ============================================================================
 */

#ifndef VNS_H_
#define VNS_H_

#include "my_structs.h"

int VNS( TreePtr phyloTreePtr, ParametersPtr paramsPtr, MultipleLeafPtr *multiple );
int exhaustiveLeafSwap( TreePtr phyloTreePtr, ParametersPtr paramsPtr, int pScore, MultipleLeafPtr *multiple );
int exhaustiveSubtreeScramble( TreePtr phyloTreePtr, ParametersPtr paramsPtr, int pScore, MultipleLeafPtr *multiple );
void showResults( TreePtr phyloTreePtr, enum distances dist, int score, double timediff );

#endif /* VNS_H_ */