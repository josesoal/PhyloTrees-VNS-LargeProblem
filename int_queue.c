/*
 ============================================================================
 Authors :
 	Jose Luis Soncco-Alvarez and Mauricio Ayala-Rincon
 	Group of Theory of Computation
 	Universidade de Brasilia (UnB) - Brazil
 ============================================================================
 */

 #include <stdio.h>
 #include <stdlib.h>

 #include "int_queue.h"

void enqueue_i( IntQueuePtr iqueuePtr, int value ) 
{
	IntQueueNodePtr newPtr;
	newPtr = malloc( sizeof( IntQueueNode ) );

	if ( newPtr != NULL ) {
		newPtr->data = value;
		newPtr->nextPtr = NULL;
		/* if empty, insert node at head */
		if ( isEmpty_i( iqueuePtr ) ) {
			iqueuePtr->headPtr = newPtr;
		}
		else {
			iqueuePtr->tailPtr->nextPtr = newPtr;
		}
		iqueuePtr->tailPtr = newPtr;
	}
	else {
		fprintf(stderr, 
            "stderr: No memory available for newPtr in enqueue_i().\n");
	}
}

int dequeue_i( IntQueuePtr iqueuePtr )
{
	int value;
	IntQueueNodePtr tempPtr;
	
	value = iqueuePtr->headPtr->data;
	tempPtr = iqueuePtr->headPtr;
	iqueuePtr->headPtr = iqueuePtr->headPtr->nextPtr;

	if ( iqueuePtr->headPtr == NULL ) { /* if queue is empty */
 		iqueuePtr->tailPtr = NULL;
	}

	free( tempPtr );
	return value;
}

/* Return 1 if the list is empty, 0 otherwise */
int isEmpty_i( IntQueuePtr iqueuePtr )
{
	return iqueuePtr->headPtr == NULL;
}