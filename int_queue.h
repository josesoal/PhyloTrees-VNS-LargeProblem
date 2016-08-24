/*
 ============================================================================
 Authors :
 	Jose Luis Soncco-Alvarez and Mauricio Ayala-Rincon
 	Group of Theory of Computation
 	Universidade de Brasilia (UnB) - Brazil
 ============================================================================
 */

 #ifndef INT_QUEUE_H_
 #define INT_QUEUE_H_

struct intQueueNode {
	int 					data;
	struct intQueueNode 	*nextPtr;
};
typedef struct intQueueNode IntQueueNode;
typedef IntQueueNode * IntQueueNodePtr;

struct intQueue {
	IntQueueNodePtr headPtr;
	IntQueueNodePtr tailPtr;
};
typedef struct intQueue IntQueue;
typedef IntQueue *IntQueuePtr; 


void enqueue_i( IntQueuePtr iqueuePtr, int value );
int dequeue_i( IntQueuePtr iqueuePtr );
int isEmpty_i( IntQueuePtr iqueuePtr );


 #endif /* INT_QUEUE_H_ */