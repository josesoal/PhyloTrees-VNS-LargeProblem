#include <stdio.h>
#include <stdlib.h>

#include "vns.h"
#include "tree.h"
#include "iterate_tree.h"
#include "random.h"

static void shakingAtNeighborhood(int neighborhood, TreePtr phyloTreePtr);
static void randomNNIforSubtrees(TreePtr phyloTreePtr);
static void randomSingleStep(TreePtr phyloTreePtr);
static void randomSubtreePR(TreePtr phyloTreePtr);
static void randomTreeBR(TreePtr phyloTreePtr);
static void disableSubtree(TreeNodePtr subtreeNodePtr, int *enabled, int *counterPtr);

int VNS(TreePtr phyloTreePtr, ParametersPtr paramsPtr )
{
	Tree temporalTree;
	int kmax = 4; /* number of neighborhood */
	int maxIterations = 30;//100
	int iteration = 1;
	int score, newScore;
	int copy = TRUE;

	/* allocate memory for temporal tree*/
	temporalTree.numberLeaves = phyloTreePtr->numberLeaves;
	temporalTree.numberGenes = phyloTreePtr->numberGenes;
	allocateMemoryForNodes( &temporalTree, paramsPtr );//--method from tree.c

	/*VNS*/
	score = labelOptimizeTree( phyloTreePtr, paramsPtr );//--method from iterate_tree.c
	while ( iteration <= maxIterations ) {
		int k = 1;
		while ( k <= kmax ) {
			/* copy current solution into temporal */
			if ( copy == TRUE ) {
				copyTreeInto( &temporalTree, phyloTreePtr, TRUE, paramsPtr );//--method from tree.c
			}	
			/* shaking (at neighborhood k) */
			shakingAtNeighborhood( k, &temporalTree );

			// dont apply  LS at the moment
			
			/* Move or not */
			newScore = labelOptimizeTree( &temporalTree, paramsPtr );//--method from iterate_tree.c

			if ( newScore < score ){
				score = newScore;
				copyTreeInto( phyloTreePtr, &temporalTree, TRUE, paramsPtr );//--method from tree.c
				copy = FALSE;// phylotree and tempotree has the same copy 
				if (DEBUG){
					printf("[(it:%d, k:%d)improved (score: %d)]", iteration,k,score);
					showTreeNewickFormat(temporalTree.startingNodePtr, SHOW_BY_ID);//--from tree.c
				}
				k = 1;//start again at neighborhood 1 
			}
			else{
				copy = TRUE;//shaking went wrong recover last tree
				k++;// go for the next neighborhood	
			}
		}

		iteration = iteration + 1;
	}

	freeTree( &temporalTree, paramsPtr );//--method from tree.c
			
	return score;
}

static void shakingAtNeighborhood(int neighborhood, TreePtr phyloTreePtr)
{
	if (neighborhood == 1){
		randomNNIforSubtrees(phyloTreePtr);
		//randomSingleStep(phyloTreePtr);
		//randomSubtreePR(phyloTreePtr);
		//randomTreeBR(phyloTreePtr);
	}
	else if (neighborhood == 2){
		randomSingleStep(phyloTreePtr);
	}
	else if (neighborhood == 3){
		randomSubtreePR(phyloTreePtr);
	}
	else if (neighborhood == 4){
		randomTreeBR(phyloTreePtr);
	}
	else{
		fprintf(stderr, " stderr: incorrect neighborhood number\n");
        exit(EXIT_FAILURE);
	} 	
}

/* (NEIGHBORHOOD 1) generates a random "Nearest Neighborhood Interchange" (NNI) for Subtrees
	(select a random "internal" edge and interchange two random subtrees) */
static void randomNNIforSubtrees(TreePtr phyloTreePtr)
{
	/* NOTE:
	* Do not consider edges linked to  a leaf.
	* The input tree must have at least 4 leaves.
	*/

	TreeNodePtr startNodePtr;
	TreeNodePtr endNodePtr;
	int i, low, hi;
	/* select a random edge that begin and end at a internal node*/
	int numberInternalNodes = phyloTreePtr->numberNodes - phyloTreePtr->numberLeaves;
	low = phyloTreePtr->numberLeaves; 
	hi = numberInternalNodes;
	while (TRUE){
		i = low + irand(hi); /* interval: numberLeaves <= i < numberNodes */
		if (phyloTreePtr->nodesPtrArray[i]->ancestorPtr->type == INTERNAL_NODE){
			/* the edge is formed by startNodePtr and endNodePtr */
			startNodePtr = phyloTreePtr->nodesPtrArray[i];
			endNodePtr = startNodePtr->ancestorPtr;
			//endNodePtr = phyloTreePtr->nodesPtrArray[i]->ancestorPtr;
			break;
		}
	}

	/* choose one (first) subtree related with the edge startNode 
	that DOES NOT contain the starting node of the tree */
	TreeNodePtr firstSubtreePtr;
	i = irand(2); /* 0 <= i <= 1 */

	if (startNodePtr->ancestorPtr->id == endNodePtr->id){
		/* if edge startNode ancestor is equal to edge endNode
		then choose the LEFT or RIGHT descendant as  1st subtree */
		if (i == 0){ 
			firstSubtreePtr = startNodePtr->leftDescPtr;
		}
		else{
			firstSubtreePtr = startNodePtr->rightDescPtr;
		}
	}	
	else if (startNodePtr->leftDescPtr->id == endNodePtr->id){
		/* if edge startNode left descendant is equal to edge endNode
		then choose the RIGHT descendant as 1st subtree 
		(not the ancestor because contains the starting node of the tree) */
		firstSubtreePtr = startNodePtr->rightDescPtr;	
	}
	else{ //startNodePtr->rightDescPtr->id == endNodePtr->id
		/* edge startNode right descendant is equal to edge endNode
		then choose the LEFT descendant as 1st subtree 
		(not the ancestor because contains the starting node of the tree) */
		firstSubtreePtr = startNodePtr->leftDescPtr;
	}

	/* choose one (second) subtree related with the edge endNode 
	that DOES NOT the starting node of the tree */
	TreeNodePtr secondSubtreePtr;
	i = irand(2); /* 0 <= i <= 1 */

	if (endNodePtr->ancestorPtr->id == startNodePtr->id){
		if (i == 0){
			secondSubtreePtr = endNodePtr->leftDescPtr;
		}
		else{
			secondSubtreePtr = endNodePtr->rightDescPtr;
		}
	}	
	else if (endNodePtr->leftDescPtr->id == startNodePtr->id){
		secondSubtreePtr = endNodePtr->rightDescPtr;		
	}
	else{ //endNodePtr->rightDescPtr->id == startNodePtr->id
		secondSubtreePtr = endNodePtr->leftDescPtr;
	}

	/* edge startNode now points to secondSubtree */
	if (startNodePtr->leftDescPtr->id == firstSubtreePtr->id){
		startNodePtr->leftDescPtr = secondSubtreePtr;	
	}
	else{//startNodePtr->rightDescPtr->id == firstSubtreePtr->id 
		startNodePtr->rightDescPtr = secondSubtreePtr;
	}

	/* edge endNode now points to firstSubtree */
	if (endNodePtr->leftDescPtr->id == secondSubtreePtr->id){
		endNodePtr->leftDescPtr = firstSubtreePtr;	
	}
	else{//startNodePtr->rightDescPtr->id == secondSubtreePtr->id 
		endNodePtr->rightDescPtr = firstSubtreePtr;
	}

	/* interchange the ancestors of subtree 1 and 2 */
	firstSubtreePtr->ancestorPtr = endNodePtr;
	secondSubtreePtr->ancestorPtr = startNodePtr;
}

/* (NEIGHBORHOOD 2)generates a random  "Single Step" (STEP)
	(select a random leaf and insert into a random edge )*/
static void randomSingleStep(TreePtr phyloTreePtr)
{
	/* NOTE
	The input tree must have at least 4 leaves */

	int i, j;
	TreeNodePtr chosenLeafPtr, edgeNodePtr, isolatedNodePtr; 
	/* choose any leaf except the root node */
	i = irand(phyloTreePtr->numberLeaves);// 0 <= i < numberLeaves
	while (phyloTreePtr->nodesPtrArray[i]->id == phyloTreePtr->startingNodePtr->id){
		i = irand(phyloTreePtr->numberLeaves);// 0 <= i < numberLeaves
	}
	chosenLeafPtr = phyloTreePtr->nodesPtrArray[i];

	/* find other descendant of chosen leaf's parent */
	TreeNodePtr otherDescPtr;
	TreeNodePtr ancestorPtr = phyloTreePtr->nodesPtrArray[i]->ancestorPtr;
	
	if (ancestorPtr->leftDescPtr != NULL && 
			ancestorPtr->leftDescPtr->id == phyloTreePtr->nodesPtrArray[i]->id){
		otherDescPtr = ancestorPtr->rightDescPtr;
	}
	else{//ancestorPtr->rightDescPtr->id == phyloTreePtr->nodesPtrArray[i]->id
		otherDescPtr = ancestorPtr->leftDescPtr;
	}

	/* re-connect other descendant to super-ancestor */
	TreeNodePtr superAncestorPtr = ancestorPtr->ancestorPtr;

	if (superAncestorPtr->leftDescPtr != NULL &&
			superAncestorPtr->leftDescPtr->id == ancestorPtr->id){

		superAncestorPtr->leftDescPtr = otherDescPtr;
		otherDescPtr->ancestorPtr = superAncestorPtr;
	}
	else{//superAncestorPtr->rightDescPtr->id == ancestorPtr->id
		superAncestorPtr->rightDescPtr = otherDescPtr;
		otherDescPtr->ancestorPtr = superAncestorPtr;
	}

	/* make internal node (ancestor of chosen leaf) isolated */
	isolatedNodePtr = ancestorPtr;
	isolatedNodePtr->ancestorPtr = NULL;
	isolatedNodePtr->leftDescPtr = NULL;
	isolatedNodePtr->rightDescPtr = NULL;

	/* choose an edge (a node) to re-link chosen leaf  */
	j = irand(phyloTreePtr->numberNodes); // 0 <= i < numberNodes
	while (TRUE){
		/* if node is not the other descendant, the isolated internal node,
		and the starting node */
		if (phyloTreePtr->nodesPtrArray[j]->id != chosenLeafPtr->id &&
			phyloTreePtr->nodesPtrArray[j]->id != otherDescPtr->id &&
			phyloTreePtr->nodesPtrArray[j]->id != ancestorPtr->id && 
			phyloTreePtr->nodesPtrArray[j]->id != phyloTreePtr->startingNodePtr->id){

			break;
		}
		j = irand(phyloTreePtr->numberNodes); // 0 <= i < numberNodes
	}
	edgeNodePtr = phyloTreePtr->nodesPtrArray[j];
	
	/* re-link the chosen leaf using the isolated internal node */	
	if (edgeNodePtr->ancestorPtr->leftDescPtr != NULL && 
		edgeNodePtr->ancestorPtr->leftDescPtr->id == edgeNodePtr->id){
		/* make isolated internal node, left descendant of ancestor of edge node */
		edgeNodePtr->ancestorPtr->leftDescPtr = isolatedNodePtr;
		isolatedNodePtr->ancestorPtr = edgeNodePtr->ancestorPtr;
		/* make edge node, left descendant of isolated internal node*/
		isolatedNodePtr->leftDescPtr = edgeNodePtr;
		edgeNodePtr->ancestorPtr = isolatedNodePtr;
		/* make chosen leaf, right descendant of isolated internal node*/
		isolatedNodePtr->rightDescPtr = chosenLeafPtr;
		chosenLeafPtr->ancestorPtr = isolatedNodePtr;
	}
	else{//edgeNodePtr->ancestorPtr->rightDescPtr->id == edgeNodePtr->id
		/* make isolated internal node, right descendant of ancestor of edge node */
		edgeNodePtr->ancestorPtr->rightDescPtr = isolatedNodePtr;
		isolatedNodePtr->ancestorPtr = edgeNodePtr->ancestorPtr;
		/* make edge node, right descendant of isolated internal node*/
		isolatedNodePtr->rightDescPtr = edgeNodePtr;
		edgeNodePtr->ancestorPtr = isolatedNodePtr;
		/* make chosen leaf, left descendant of isolated internal node*/
		isolatedNodePtr->leftDescPtr = chosenLeafPtr;
		chosenLeafPtr->ancestorPtr = isolatedNodePtr;
	}
}

/* (NEIGHBORHOOD 3) generates a random "Subtree Pruning and Regrafting" (SPR) 
	(select a random internal node and eliminate its 3 edges,
	then reconnect two random subtress and finally reconnect 
	the last subtree to a random edge) */
static void randomSubtreePR(TreePtr phyloTreePtr)
{
	/* NOTE
	The input tree must have at least 4 leaves */
	
	int i,j,low,hi;
	TreeNodePtr chosenInternalNodePtr, ancestorPtr, leftDescPtr, rightDescPtr;

	/* choose a random internal node  */
	int numberInternalNodes = phyloTreePtr->numberNodes - phyloTreePtr->numberLeaves;
	low = phyloTreePtr->numberLeaves; 
	hi = numberInternalNodes;
	i = low + irand(hi); /* interval: numberLeaves <= i < numberNodes */
	chosenInternalNodePtr = phyloTreePtr->nodesPtrArray[i];

	/* find 3 nodes (subtrees) linked to the chosen internal node*/
	ancestorPtr = chosenInternalNodePtr->ancestorPtr;
	leftDescPtr = chosenInternalNodePtr->leftDescPtr;
	rightDescPtr = chosenInternalNodePtr->rightDescPtr;

	/* find the new edge node to be joined with 
	the ancestor of the chosen internal node */
	i = irand(2); /* interval: 0 <= i < 2*/
	TreeNodePtr newEdgeNodePtr, remainingNodePtr;

	if (i == 0){
		if (leftDescPtr->type == INTERNAL_NODE){
			newEdgeNodePtr = leftDescPtr;
			remainingNodePtr = rightDescPtr;
		}
		else{//rightDescPtr->type == INTERNAL_NODE
			newEdgeNodePtr = rightDescPtr;
			remainingNodePtr = leftDescPtr;
		}
	}
	else{
		if (rightDescPtr->type == INTERNAL_NODE){
			newEdgeNodePtr = rightDescPtr;
			remainingNodePtr = leftDescPtr;
		}
		else{//leftDescPtr->type == INTERNAL_NODE
			
			newEdgeNodePtr = leftDescPtr;
			remainingNodePtr = rightDescPtr;
		}
	}

	/* join the subtree that was ancestor of internal node
	with other subtree whose root is an internal node */
	if (ancestorPtr->leftDescPtr != NULL && 
			ancestorPtr->leftDescPtr->id == chosenInternalNodePtr->id){
		ancestorPtr->leftDescPtr = newEdgeNodePtr;
		newEdgeNodePtr->ancestorPtr = ancestorPtr;
	}
	else{//ancestorPtr->rightDescPtr->id == chosenInternalNodePtr->id
		ancestorPtr->rightDescPtr = newEdgeNodePtr;
		newEdgeNodePtr->ancestorPtr = ancestorPtr;
	}

	/* eliminate 3 edges of chosen internal node */
	chosenInternalNodePtr->ancestorPtr = NULL;
	chosenInternalNodePtr->leftDescPtr = NULL;
	chosenInternalNodePtr->rightDescPtr = NULL;

	/* select the last edge (node) for joining with the remaining subtree */
	int enabled[phyloTreePtr->numberNodes];// malloc?
	int counter = phyloTreePtr->numberNodes;
	for (i=0; i<phyloTreePtr->numberNodes; i++){
		enabled[i] = TRUE;
	}

	enabled[phyloTreePtr->startingNodePtr->index] = FALSE; 
	enabled[newEdgeNodePtr->index] = FALSE;
	enabled[chosenInternalNodePtr->index] = FALSE;
	counter = counter - 3;
	disableSubtree(remainingNodePtr, enabled, &counter);//disable == make false

	i = irand(counter);/* interval: 0 <= i < counter */
	TreeNodePtr lastEdgeNodePtr;
	for (j=0; j<phyloTreePtr->numberNodes; j++){
		if (enabled[j] == TRUE){
			if (i == 0){
				lastEdgeNodePtr = phyloTreePtr->nodesPtrArray[j];
				break;
			}
			i = i - 1;
		}
	}	

	/* join the remaining subtree to the last edge (node) selected,
	except the new edge created*/
	if (lastEdgeNodePtr->ancestorPtr->leftDescPtr != NULL && 
			lastEdgeNodePtr->ancestorPtr->leftDescPtr->id == lastEdgeNodePtr->id){
		lastEdgeNodePtr->ancestorPtr->leftDescPtr = chosenInternalNodePtr;
		chosenInternalNodePtr->ancestorPtr = lastEdgeNodePtr->ancestorPtr;

		chosenInternalNodePtr->leftDescPtr = lastEdgeNodePtr;
		lastEdgeNodePtr->ancestorPtr = chosenInternalNodePtr;

		chosenInternalNodePtr->rightDescPtr = remainingNodePtr;
		remainingNodePtr->ancestorPtr = chosenInternalNodePtr;
	}
	else{//lastEdgeNodePtr->ancestorPtr->rightDescPtr->id == lastEdgeNodePtr->id
		lastEdgeNodePtr->ancestorPtr->rightDescPtr = chosenInternalNodePtr;
		chosenInternalNodePtr->ancestorPtr = lastEdgeNodePtr->ancestorPtr;

		chosenInternalNodePtr->rightDescPtr = lastEdgeNodePtr;
		lastEdgeNodePtr->ancestorPtr = chosenInternalNodePtr;

		chosenInternalNodePtr->leftDescPtr = remainingNodePtr;
		remainingNodePtr->ancestorPtr = chosenInternalNodePtr;
	}

}

/* (NEIGHBORHOOD 4)generated a random "Tree Bisection and Reconnection" (TBR) 
	(eliminate a random "internal" edge and 
	reconnect the two subtrees by a new random edge)
*/
static void randomTreeBR(TreePtr phyloTreePtr)
{
	/* NOTE
	- The input tree must have at least 4 leaves 
	- For a tree with 4 leaves just exists one neighbor */

	if (phyloTreePtr->numberLeaves == 4) return;

	TreeNodePtr startNodePtr,endNodePtr, otherDescendantPtr;
	int i, j, low, hi, numberInternalNodes;
	
	/* select a random edge that begin and end at a internal node */
	numberInternalNodes = phyloTreePtr->numberNodes - phyloTreePtr->numberLeaves;
	low = phyloTreePtr->numberLeaves; 
	hi = numberInternalNodes;
	while (TRUE){
		i = low + irand(hi); /* interval: numberLeaves <= i < numberNodes */
		
		if (phyloTreePtr->nodesPtrArray[i]->ancestorPtr->type == INTERNAL_NODE){
			/* the edge is formed by startNodePtr and endNodePtr */
			startNodePtr = phyloTreePtr->nodesPtrArray[i];
			endNodePtr = startNodePtr->ancestorPtr;
			
			if (endNodePtr->leftDescPtr->id == startNodePtr->id){
				otherDescendantPtr = endNodePtr->rightDescPtr;
			}
			else{//endNodePtr->rightDescPtr->id == startNodePtr->id
				otherDescendantPtr = endNodePtr->leftDescPtr;
			}

			/* avoid the case where ancestor of endNodePtr is the startingNodePtr,
			and the otherDescendant is a leaf node */
			if (endNodePtr->ancestorPtr->id == phyloTreePtr->startingNodePtr->id && 
				otherDescendantPtr->type == LEAF_NODE){
			}
			else{
				break;
			}
		}
	}

	/* join ancestor of end node to the other 
	descendant (otherDescendantPtr) different from start node */
	if (endNodePtr->ancestorPtr->leftDescPtr != NULL && 
			endNodePtr->ancestorPtr->leftDescPtr->id == endNodePtr->id){
		endNodePtr->ancestorPtr->leftDescPtr = otherDescendantPtr;
		
	}
	else{//endNodePtr->ancestorPtr->rightDescPtr->id == endNodePtr->id
		endNodePtr->ancestorPtr->rightDescPtr = otherDescendantPtr;
	}
	otherDescendantPtr->ancestorPtr = endNodePtr->ancestorPtr;

	/* leave the end node (endNodePtr) isolated */
	endNodePtr->ancestorPtr = NULL;
	endNodePtr->leftDescPtr = NULL;
	endNodePtr->rightDescPtr = NULL;

	/* choose a new node (start of new edge node) for inserting 
	the remaining subtree  by the start node (startNodePtr).

	NOTE: the new node must be different from:
	- startingNodePtr
	- otherDescendantPtr 
	- endNodePtr 
	- any node from the remaining subtree (starNodePtr) */

	int enabled[phyloTreePtr->numberNodes];// malloc?
	int counter = phyloTreePtr->numberNodes;

	for (i=0; i<phyloTreePtr->numberNodes; i++){
		enabled[i] = TRUE;
	}

	enabled[phyloTreePtr->startingNodePtr->index] = FALSE; 
	enabled[otherDescendantPtr->index] = FALSE;
	enabled[endNodePtr->index] = FALSE;
	counter = counter - 3;
	disableSubtree(startNodePtr, enabled, &counter);//disable == make false

	i = irand(counter);/* interval: 0 <= i < counter */
	TreeNodePtr newEdgeNodePtr;
	for (j=0; j<phyloTreePtr->numberNodes; j++){
		if (enabled[j] == TRUE){
			if (i == 0){
				newEdgeNodePtr = phyloTreePtr->nodesPtrArray[j];
				break;
			}
			i = i - 1;
		}
	}	

	/* Insert into the new edge node (newEdgeNodePtr),
	the start node (startNodePtr) by using 
	the isolated end node(endNodePtr) as internal node */
	TreeNodePtr internalNodePtr = endNodePtr;

	if (newEdgeNodePtr->ancestorPtr->leftDescPtr != NULL &&
			newEdgeNodePtr->ancestorPtr->leftDescPtr->id == newEdgeNodePtr->id){
		newEdgeNodePtr->ancestorPtr->leftDescPtr = internalNodePtr;
		internalNodePtr->ancestorPtr = newEdgeNodePtr->ancestorPtr;

		internalNodePtr->leftDescPtr = newEdgeNodePtr;
		newEdgeNodePtr->ancestorPtr = internalNodePtr;

		internalNodePtr->rightDescPtr = startNodePtr;
		startNodePtr->ancestorPtr = internalNodePtr;
	}
	else{//newEdgeNodePtr->ancestorPtr->rightDescPtr->id == newEdgeNodePtr->id
		newEdgeNodePtr->ancestorPtr->rightDescPtr = internalNodePtr;
		internalNodePtr->ancestorPtr = newEdgeNodePtr->ancestorPtr;

		internalNodePtr->rightDescPtr = newEdgeNodePtr;
		newEdgeNodePtr->ancestorPtr = internalNodePtr;

		internalNodePtr->leftDescPtr = startNodePtr;
		startNodePtr->ancestorPtr = internalNodePtr;
	}
}

static void disableSubtree(TreeNodePtr subtreeNodePtr, int *enabled, int *counterPtr)
{
	enabled[subtreeNodePtr->index] = FALSE;
	(*counterPtr) = (*counterPtr) -1;
	if (subtreeNodePtr->leftDescPtr != NULL) 
		disableSubtree(subtreeNodePtr->leftDescPtr, enabled, counterPtr);
	if (subtreeNodePtr->rightDescPtr != NULL)
		disableSubtree(subtreeNodePtr->rightDescPtr, enabled, counterPtr);
}

/* (EXHAUSTIVE MUTATION 1) 
return the best score from a exhaustive search by swapping two leaves */
int exhaustiveLeafSwap( TreePtr phyloTreePtr, ParametersPtr paramsPtr, int pScore )
{
	int i,j,score,newScore,copy, restart;
	Tree temporalTree;
	TreeNodePtr node1Ptr, node2Ptr, internalNode1Ptr, internalNode2Ptr;

	/* allocate memory for temporal tree*/
	temporalTree.numberLeaves = phyloTreePtr->numberLeaves;
	temporalTree.numberGenes = phyloTreePtr->numberGenes;
	allocateMemoryForNodes( &temporalTree, paramsPtr );//--from tree.c

	restart = TRUE;
	copy = TRUE;
	//score = labelOptimizeTree( phyloTreePtr, paramsPtr );//from iterate_tree.c
	score = pScore;
	while ( TRUE ) {
		if ( restart == TRUE ) {
			i = 0; 
			j = i + 1;
		}
		else{
			i++; 
			j = i + 1;
			if ( i == phyloTreePtr->numberLeaves-1 ) {
				break;
			}
		}

		/* copy current solution into temporal */
		if ( copy == TRUE ) {
			copyTreeInto( &temporalTree, phyloTreePtr, TRUE, paramsPtr );//--from tree.c
		}
		/* pick two nodes for swapping (avoiding picking the root node and 
		having two leaf nodes with the same ancestor) */
		node1Ptr = temporalTree.nodesPtrArray[ i ];
		node2Ptr = temporalTree.nodesPtrArray[ j ];
		if ( node1Ptr->id == temporalTree.startingNodePtr->id || 
			node2Ptr->id == temporalTree.startingNodePtr->id ||
			node1Ptr->ancestorPtr->id == node2Ptr->ancestorPtr->id ) {

			copy = FALSE;
			restart = FALSE;
			continue;
		}

		/* unlink the node1 and node2 */
		unlinkNode( node1Ptr, &internalNode1Ptr );//--from tree.c
		unlinkNode( node2Ptr, &internalNode2Ptr );//--from tree.c

		/* relink (by swapping) node1 and node2 */
		relinkNode( node1Ptr, internalNode2Ptr );//--from tree.c
		relinkNode( node2Ptr, internalNode1Ptr );//--from tree.c

		/* calculate newScore */ 
		newScore = labelOptimizeTree( &temporalTree, paramsPtr );//from iterate_tree.c

		if ( newScore < score ) {
			score = newScore;
			copyTreeInto( phyloTreePtr, &temporalTree, TRUE, paramsPtr );//--from tree.c
			copy = FALSE;
			restart = TRUE; /* restart the search */
			if ( DEBUG ) {
				printf("[tree refined by SWAP (score: %d)]",score);
				showTreeNewickFormat(temporalTree.startingNodePtr, SHOW_BY_ID);//--from tree.c
			}
		}
		else { 
			copy = TRUE; /* recover previous state */
			restart = FALSE; 
		}
	}
	
	freeTree( &temporalTree, paramsPtr );//--method from tree.c
	return score;
}

/* (EXHAUSTIVE MUTATION 2) 
return the best score from a exhaustive scramble over subtrees */
int exhaustiveSubtreeScramble(TreePtr phyloTreePtr, ParametersPtr paramsPtr, int pScore )
{
	int i, k, score, score1, score2, score3, score4, restart;
	Tree temporalTree1, temporalTree2, temporalTree3, temporalTree4;
	TreeNodePtr subtreeRootPtr, internalNode1Ptr, internalNode2Ptr; 
	TreeNodePtr nodeL_Ptr, nodeRL_Ptr, nodeRR_Ptr, nodeR_Ptr, nodeLL_Ptr, nodeLR_Ptr;

	/* allocate memory for temporal trees */
	temporalTree1.numberLeaves = phyloTreePtr->numberLeaves;
	temporalTree2.numberLeaves = phyloTreePtr->numberLeaves;
	temporalTree3.numberLeaves = phyloTreePtr->numberLeaves;
	temporalTree4.numberLeaves = phyloTreePtr->numberLeaves;
	temporalTree1.numberGenes = phyloTreePtr->numberGenes;
	temporalTree2.numberGenes = phyloTreePtr->numberGenes;
	temporalTree3.numberGenes = phyloTreePtr->numberGenes;
	temporalTree4.numberGenes = phyloTreePtr->numberGenes;
	allocateMemoryForNodes( &temporalTree1, paramsPtr );//--from tree.c
	allocateMemoryForNodes( &temporalTree2, paramsPtr );//--from tree.c
	allocateMemoryForNodes( &temporalTree3, paramsPtr );//--from tree.c
	allocateMemoryForNodes( &temporalTree4, paramsPtr );//--from tree.c

	/* choose an edge that starts at an internal edge */
	//score = labelOptimizeTree( phyloTreePtr, paramsPtr );//from iterate_tree.c
	score = pScore;
	restart = TRUE;
	while (TRUE){
		if (restart == TRUE){
			i = phyloTreePtr->numberLeaves;
		}
		else{
			i++;
			if (i == phyloTreePtr->numberNodes) break;
		}

		/* NOTE: the subtree must have at least one internal node, 
			otherwise the scramble can not be applied */

		/* explore the left descendant of the subtree */
		score1 = score; score2 = score;
		subtreeRootPtr = phyloTreePtr->nodesPtrArray[i];
		if (subtreeRootPtr->leftDescPtr->type == INTERNAL_NODE){
			/* [swap THE RIGHT DESCENDANT of the subtree with THE LEFT DESCENDANT 
			of the left descendant of the subtree] */	
			copyTreeInto(&temporalTree1, phyloTreePtr, TRUE, paramsPtr);
			subtreeRootPtr = temporalTree1.nodesPtrArray[i];
			nodeR_Ptr = subtreeRootPtr->rightDescPtr;
			nodeLL_Ptr = subtreeRootPtr->leftDescPtr->leftDescPtr;
			/* unlink the RIGHT DESCENDANT of subtree root and 
			LEFT DESCENDANT of the left descedant of subtree root */
			unlinkNode(nodeR_Ptr, &internalNode1Ptr);
			unlinkNode(nodeLL_Ptr, &internalNode2Ptr);
			/* relink the nodes by swapping them */
			relinkNode(nodeR_Ptr, internalNode2Ptr);
			relinkNode(nodeLL_Ptr, internalNode1Ptr);
			score1 = labelOptimizeTree( &temporalTree1, paramsPtr );//from iterate_tree.c

			/* [swap THE RIGHT DESCENDANT of the subtree with THE RIGHT DESCENDANT 
			of the left descendant of the subtree] */
			copyTreeInto(&temporalTree2, phyloTreePtr, TRUE, paramsPtr);
			subtreeRootPtr = temporalTree2.nodesPtrArray[i];
			nodeR_Ptr = subtreeRootPtr->rightDescPtr;
			nodeLR_Ptr = subtreeRootPtr->leftDescPtr->rightDescPtr;
			/* unlink the RIGHT DESCENDANT of subtree root and 
			RIGHT DESCENDANT of the left descedant of subtree root */
			unlinkNode(nodeR_Ptr, &internalNode1Ptr);
			unlinkNode(nodeLR_Ptr, &internalNode2Ptr);
			/* relink the nodes by swapping them */
			relinkNode(nodeR_Ptr, internalNode2Ptr);
			relinkNode(nodeLR_Ptr, internalNode1Ptr);
			score2 = labelOptimizeTree( &temporalTree2, paramsPtr );//from iterate_tree.c
		}

		/* explore the right descendant of the subtree */
		score3 = score; score4 = score;
		subtreeRootPtr = phyloTreePtr->nodesPtrArray[i];
		if (subtreeRootPtr->rightDescPtr->type == INTERNAL_NODE){
			/* [swap THE LEFT DESCENDANT of the subtree with THE LEFT DESCENDANT 
			of the right descendant of the subtree] */
			copyTreeInto(&temporalTree3, phyloTreePtr, TRUE, paramsPtr);
			subtreeRootPtr = temporalTree3.nodesPtrArray[i];
			nodeL_Ptr = subtreeRootPtr->leftDescPtr;
			nodeRL_Ptr = subtreeRootPtr->rightDescPtr->leftDescPtr;

			/* unlink the LEFT DESCENDANT of subtree root and 
			LEFT DESCENDANT of the right descedant of subtree root */
			unlinkNode(nodeL_Ptr, &internalNode1Ptr);
			unlinkNode(nodeRL_Ptr, &internalNode2Ptr);

			/* relink the nodes by swapping them */
			relinkNode(nodeL_Ptr, internalNode2Ptr);
			relinkNode(nodeRL_Ptr, internalNode1Ptr);
			score3 = labelOptimizeTree( &temporalTree3, paramsPtr );//from iterate_tree.c

			/* [swap THE LEFT DESCENDANT of the subtree with THE RIGHT DESCENDANT 
			of the right descendant of the subtree] */
			copyTreeInto(&temporalTree4, phyloTreePtr, TRUE, paramsPtr);
			subtreeRootPtr = temporalTree4.nodesPtrArray[i];
			nodeL_Ptr = subtreeRootPtr->leftDescPtr;
			nodeRR_Ptr = subtreeRootPtr->rightDescPtr->rightDescPtr;
			/* unlink the LEFT DESCENDANT of subtree root and 
			RIGHT DESCENDANT of the right descedant of subtree root */
			unlinkNode(nodeL_Ptr, &internalNode1Ptr);
			unlinkNode(nodeRR_Ptr, &internalNode2Ptr);
			/* relink the nodes by swapping them */
			relinkNode(nodeL_Ptr, internalNode2Ptr);
			relinkNode(nodeRR_Ptr, internalNode1Ptr);
			score4 = labelOptimizeTree( &temporalTree4, paramsPtr );//from iterate_tree.c
		}

		/* determine which tree is the best */
		k = 0; restart = FALSE;
		if (score1 < score){
			score = score1; k = 1; restart = TRUE;
		}
		if (score2 < score){
			score = score2; k =2; restart = TRUE;
		}
		if (score3 < score){
			score = score3; k =3; restart = TRUE;
		}
		if (score4 < score){
			score = score4; k =4; restart = TRUE;
		}
		/* update tree */
		if (k == 1) copyTreeInto(phyloTreePtr, &temporalTree1, TRUE, paramsPtr);
		else if (k == 2) copyTreeInto(phyloTreePtr, &temporalTree2, TRUE, paramsPtr);
		else if (k == 3) copyTreeInto(phyloTreePtr, &temporalTree3, TRUE, paramsPtr);
		else if (k == 4) copyTreeInto(phyloTreePtr, &temporalTree4, TRUE, paramsPtr);

		if (DEBUG){
			if (k > 0){
				printf("[tree refined by SCRAMBLE (score: %d)]",score);
				showTreeNewickFormat(phyloTreePtr->startingNodePtr, SHOW_BY_ID);//--from tree.c
			}
		}
	}

	freeTree( &temporalTree1, paramsPtr );//--method from tree.c
	freeTree( &temporalTree2, paramsPtr );//--method from tree.c
	freeTree( &temporalTree3, paramsPtr );//--method from tree.c
	freeTree( &temporalTree4, paramsPtr );//--method from tree.c

	return score;
}

/* (MUTATION 1) 
	perform a SWAP over two random leaves */
/*static*/ void randomLeafSwap(TreePtr phyloTreePtr, 
	enum medianSolvers solver, enum distances distanceType, int circular)
{

}

/* (MUTATION 2)
	perform a SCRAMBLE over a random subtree  */
/*static*/ void randomSubtreeScramble(TreePtr phyloTreePtr, 
	enum medianSolvers solver, enum distances distanceType, int circular)
{

}

void showResults(TreePtr phyloTreePtr, enum distances dist, int score, double timediff)
{

	printf("-----------------------------\n");
    printf("VNS for Large Phylogeny\n");
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
    showTreeNewickFormat(phyloTreePtr->startingNodePtr, SHOW_BY_ID);//--method from tree.c
    showTreeNewickFormat(phyloTreePtr->startingNodePtr, SHOW_BY_NAME);//--method from tree.c
    //showNodesArray(phyloTreePtr);//--method from tree.c
    writeNewickFormatToFile( "newick_tree.txt", phyloTreePtr->startingNodePtr, SHOW_BY_NAME );//--from tree.c 
}

