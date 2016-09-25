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

/* this auxiliary structure is used to read genomes from a newick tree format */
struct genomeNode {
	char name[ MAX_STRING_LEN ];
    int visited;
};
typedef struct genomeNode GNode;
typedef GNode *GNodePtr;


void allocateMemoryForNodes( TreePtr phyloTreePtr, ParametersPtr paramsPtr );
void freeTree( TreePtr phyloTreePtr, ParametersPtr paramsPtr );
void copyTreeInto( TreePtr phyloTree1Ptr, 
	TreePtr phyloTree2Ptr, int copyStructure, ParametersPtr paramsPtr);
void readNumberLeavesAndGenes( 
    TreePtr phyloTreePtr, ParametersPtr paramsPtr, RawDatasetPtr rdatasetPtr ); 
void readGenomesFromRawData( TreePtr phyloTreePtr, ParametersPtr paramsPtr, 
            RawDatasetPtr rdatasetPtr, MultipleLeafPtr **multiple );
void allocateMemoryForLeafCandidates( 
                TreePtr phyloTreePtr, MultipleLeafPtr *multiple );
void freeMultipleLeafs( 
	TreePtr phyloTreePtr, MultipleLeafPtr **multiple, ParametersPtr paramsPtr );

void createTopologyFromNewickFormat( TreePtr phyloTreePtr, ParametersPtr paramsPtr );
int createInitialTreeTopology(TreePtr phyloTreePtr, ParametersPtr paramsPtr, MultipleLeafPtr *multiple );

void unlinkNode(TreeNodePtr nodePtr, TreeNodePtr *internalNodePtr);
void relinkNode(TreeNodePtr nodePtr, TreeNodePtr internalNodePtr);

void showTreeNewickFormat(TreeNodePtr currentNodePtr, int option);
void writeNewickFormatToFile( char *filename, TreeNodePtr currentNodePtr, int option );

void showNode( TreePtr phyloTreePtr, TreeNodePtr nodePtr );
void showNodesArray( TreePtr phyloTreePtr );
void showGenomes( TreePtr phyloTreePtr, ParametersPtr paramsPtr );
void showResultsSmallPhylogeny( TreePtr phyloTreePtr, 
                        enum distances dist, int score, double timediff );

#endif /* TREE_H_ */