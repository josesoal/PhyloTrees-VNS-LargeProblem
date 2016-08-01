/*
 ============================================================================
 Authors :
 	Jose Luis Soncco-Alvarez and Mauricio Ayala-Rincon
 	Group of Theory of Computation
 	Universidade de Brasilia (UnB) - Brazil
 ============================================================================
 */

#ifndef TREE_H_
#define TREE_H_

#include "my_structs.h"
#include "condense.h"

void allocateMemoryForNodes( TreePtr phyloTreePtr, ParametersPtr paramsPtr );
void freeTree( TreePtr phyloTreePtr, ParametersPtr paramsPtr );
void copyTreeInto( TreePtr phyloTree1Ptr, 
	TreePtr phyloTree2Ptr, int copyStructure, ParametersPtr paramsPtr);
void readGenomesFromRawData( 
    TreePtr phyloTreePtr, ParametersPtr paramsPtr, RawDatasetPtr rdatasetPtr );
int createInitialTreeTopology(TreePtr phyloTreePtr, ParametersPtr paramsPtr );

void unlinkNode(TreeNodePtr nodePtr, TreeNodePtr *internalNodePtr);
void relinkNode(TreeNodePtr nodePtr, TreeNodePtr internalNodePtr);

void showTreeNewickFormat(TreeNodePtr currentNodePtr, int option);
void writeNewickFormatToFile( char *filename, TreeNodePtr currentNodePtr, int option );

void showNode( TreePtr phyloTreePtr, TreeNodePtr nodePtr );
void showNodesArray( TreePtr phyloTreePtr );
void showGenomes( TreePtr phyloTreePtr, ParametersPtr paramsPtr );

#endif /* TREE_H_ */