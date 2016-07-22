/*
 ============================================================================
 Authors :
 	Jose Luis Soncco-Alvarez and Mauricio Ayala-Rincon
 	Group of Theory of Computation
 	Universidade de Brasilia (UnB) - Brazil
 ============================================================================
 */

#ifndef MEDIAN_SOLVERS_H_
#define MEDIAN_SOLVERS_H_

#include "structs.h"

void createTspInstance(int numberGenes, int *genome1, int *genome2, int *genome3, 
    struct adj_struct *adjacencyList);
int callCoalestsp(int numberGenes, int *genome1, int *genome2, int *genome3, 
    struct adj_struct *adjList, int *medianGenome, int *tspCycle, int circular);
int callBbtsp(int numberGenes, int *genome1, int *genome2, int *genome3, 
    struct adj_struct *adjList, int *medianGenome, int *tspCycle, int circular);
void callCapraraInversionMedian(int numberGenes,
    int *genome1, int *genome2, int *genome3, int *medianGenome, int circular);


#endif /* MEDIAN_SOLVERS_H_ */