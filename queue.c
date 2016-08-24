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

#include "queue.h"

void enqueue( QueuePtr myQueuePtr, TreeNodePtr valuePtr )
{
    QueueNodePtr newPtr;
    newPtr = malloc( sizeof( QueueNode ) );
    
    if ( newPtr != NULL ) {
        newPtr->dataPtr = valuePtr;
        newPtr->nextPtr = NULL;
        /* if empty, insert node at head */
        if ( isEmpty( myQueuePtr ) ) {
            myQueuePtr->headPtr = newPtr;
        }
        else {
            myQueuePtr->tailPtr->nextPtr = newPtr;
        }
        myQueuePtr->tailPtr = newPtr;
    }
    else {
        fprintf(stderr, 
            "stderr: No memory available for newPtr in enqueue().\n");
    }
}

TreeNodePtr dequeue( QueuePtr myQueuePtr )
{
    TreeNodePtr valuePtr; 
    QueueNodePtr tempPtr;
    
    valuePtr = myQueuePtr->headPtr->dataPtr;
    tempPtr = myQueuePtr->headPtr;
    myQueuePtr->headPtr = myQueuePtr->headPtr->nextPtr;
    
    if ( myQueuePtr->headPtr == NULL ) { /* if queue is empty */
        myQueuePtr->tailPtr = NULL;
    }
    
    free( tempPtr );
    return valuePtr;
}

/* Return 1 if the list is empty, 0 otherwise */
int isEmpty( QueuePtr myQueuePtr )
{
    return myQueuePtr->headPtr == NULL;
}