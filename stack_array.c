/*
 ============================================================================
 Authors :
    Jose Luis Soncco-Alvarez and Mauricio Ayala-Rincon
    Group of Theory of Computation
    Universidade de Brasilia (UnB) - Brazil
 ============================================================================
 */
#include <stdio.h>
#include <string.h>
#include "stack_array.h"

void push( StackPtr stackPtr, char * item ) {
    stackPtr->top++;
    strcpy( stackPtr->s[ stackPtr->top ], item );
}

int isStackEmpty( StackPtr stackPtr ) {
    if ( stackPtr->top <= -1 ) return 1;
    return 0;
}

char* pop( StackPtr stackPtr ) {
    char *item;
    item = stackPtr->s[ stackPtr->top ];
    stackPtr->top--;
    return item;
}

void show( StackPtr stackPtr ) {
    int i;
    if ( isStackEmpty( stackPtr ) )
        printf( "Stack Is Empty!\n" );
    else {
        for ( i = stackPtr->top; i >= 0; i-- )
            printf( "%s,", stackPtr->s[i] );
        printf("\n");
    }
}

/* Before using this data structure init a empty stack:
    Stack st1;
    st1.top = -1;
*/