/*
 ============================================================================
 Authors :
 	Jose Luis Soncco-Alvarez and Mauricio Ayala-Rincon
 	Group of Theory of Computation
 	Universidade de Brasilia (UnB) - Brazil
 ============================================================================
 */

#ifndef MY_STRUCTS_H_
#define MY_STRUCTS_H_

#define MAX_STRING_LEN	20
#define TRUE 			1
#define FALSE 			0
#define LEAF_NODE 		0			/* leaf node */
#define INTERNAL_NODE 	1 			/* internal node */
#define EACH_NODE		2 			/* leaf and internal nodes */
#define ANCESTOR 		0			/* ancestor */
#define LEFT_DESC 		1      		/* left descendant */
#define RIGHT_DESC 		2     		/* right descendant */
#define SHOW_BY_ID 		0
#define SHOW_BY_NAME 	1
#define ADJACENCY		0
#define TELOMERE 		1

#define GO_LEFT			0
#define GO_RIGHT		1
#define DEBUG 			TRUE		/* TRUE for showing the process */
#define SHOW_JUST_SCORE	TRUE 		/* TRUE for showing just a single number (score) */

enum distances {INVERSION_DIST, BREAKPOINT_DIST, DCJ_DIST};
enum medianSolvers {COALESTSP, 			/* Approximate TSP solver */
                     BBTSP, 			/* Exact TSP solver */
                     CAPRARA_INV_MEDIAN /* Branch & Bound Solver for Inv. Median */
                    };
enum initTreeMethod {R_LEAF_R_EDGE, R_LEAF_1BEST_EDGE};
enum optimizer {BLANCHETTE, KOVAC, GREEDY_CANDIDATES};
enum penalty {MULTIPLE_CH, MULT_CIRCULAR_CH, LIN_CH_MULT_CIRC_CH, COMB_LIN_CIRC_CH};

struct parameters {
	unsigned int 			seed;
	char 			       	*testsetName;
	char					*newickFile;
	int 					unichromosomes; /* false if using multichromosomes genomes */
    int                    	circular; 		/* true if the dataset genomes are circular */
    enum distances         	distanceType;
    enum medianSolvers     	solver;
    int                    	useOutgroup;
    char                   	*outgroup;
    enum penalty 			penaltyType;
    enum initTreeMethod    	initMethod;
    enum optimizer         	opt;             /* optimizer for small phylogeny problem */
};
typedef struct parameters 	Parameters;
typedef Parameters 			*ParametersPtr;

/*	Adjacency: x and y are differents
	Telomere: x and y have the same value */
struct pointDCJ {
	int x;
	int y;
	int type; /* ADJACENCY or TELOMERE */ 
};
typedef struct pointDCJ 		PointDCJ;
typedef PointDCJ 				*PointDCJPtr;

struct candidate {
	int 			*genome;
	PointDCJPtr		*genomeDCJ;
	int 			*inverseDCJ;
	int 			numPointsDCJ;
	int 			numLinearCh; /* number of linear chromosomes */
	int 			numCircularCh; /* number of circular chromosomes */	 
};
typedef struct candidate Candidate;
typedef Candidate 		*CandidatePtr;

struct treeNode{
	int 				id;
	int 				index;			/* index at the nodesPtrArray */
	char 				*organism; 		/* Organism name of a LEAF_NODE */
	int 				type; 			/* LEAF_NODE or INTERNAL_NODE */
	int 				*genome;		/* genome formed by integers */
	PointDCJPtr			*genomeDCJ;		/* genome formed by pointers to points (adjaceny ot telomere)*/
	int 				*inverseDCJ;	/* inverse of the genomeDCJ */
	int 				numPointsDCJ; 	/* number of elements of genomeDCJ */
	int 				edgeWeight; 	/* egde between current node and ancestor*/
	struct treeNode 	*ancestorPtr;
	struct treeNode 	*leftDescPtr;
	struct treeNode 	*rightDescPtr;
	int 				avaliable; 		/* TRUE, if node is not being used in the tree */
	int 				extremity;		/* used for the DFS over the tree */
										/* extremity==TRUE denotes a node already labeled */
};
typedef struct treeNode TreeNode;
typedef TreeNode 		*TreeNodePtr;

struct tree{
	TreeNodePtr 		startingNodePtr;
	TreeNodePtr			*nodesPtrArray;	/* array of pointers to the tree nodes */
	int 				numberGenes;	/* number of genes of each genome */
	int 				numberLeaves;	/* equal to number of genomes*/
	int 				numberNodes;	/* equal to (leaves + internal) nodes */

	int 				avaliableLeavesNodes; /* number of avaliable LEAF nodes */
	int 				avaliableInternalNodes; /* number of avaliable internal nodes */
};

typedef struct tree 	Tree;
typedef Tree 			*TreePtr;

/* structs for queue */
struct queueNode {
	struct treeNode 	*dataPtr;  
	struct queueNode 	*nextPtr; 
};
typedef struct queueNode QueueNode;
typedef QueueNode 		*QueueNodePtr;

struct queue {
    QueueNodePtr 		headPtr;
    QueueNodePtr 		tailPtr;
};
typedef struct queue 	Queue;
typedef Queue 			*QueuePtr;

#endif /* MY_STRUCTS_H_ */







