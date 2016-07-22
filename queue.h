/*
 ============================================================================
 Authors :
 	Jose Luis Soncco-Alvarez and Mauricio Ayala-Rincon
 	Group of Theory of Computation
 	Universidade de Brasilia (UnB) - Brazil
 ============================================================================
 */
 
#ifndef QUEUE_H_
#define QUEUE_H_ 

#include "my_structs.h"

void enqueue( QueuePtr myQueuePtr, TreeNodePtr valuePtr );
TreeNodePtr dequeue( QueuePtr myQueuePtr );
int isEmpty( QueuePtr myQueuePtr );

#endif /* QUEUE_H_ */