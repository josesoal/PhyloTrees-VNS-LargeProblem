#ifndef LABEL_TREE_H_
#define LABEL_TREE_H_

#include "structs.h"
#include "my_structs.h"

int labelOptimizeTree( TreePtr phyloTreePtr, ParametersPtr paramsPtr );
int scoreTree( TreePtr phyloTreePtr, TreeNodePtr nodePtr, ParametersPtr paramsPtr );
int calculateDistance( int *genome1, int *genome2, int numberGenes, ParametersPtr paramsPtr );

#endif /* LABEL_TREE_H_ */