  /*
 ============================================================================
 Authors :
    Jose Luis Soncco-Alvarez and Mauricio Ayala-Rincon
    Group of Theory of Computation
    Universidade de Brasilia (UnB) - Brazil
 ============================================================================
 */

#ifndef STACK_ARRAY_
#define STACK_ARRAY_

#include "my_structs.h"

/* These constants are used just for the Small-Phylogeny case */
#define STACK_SIZE 100
#define MAX_NEWICK_LEN 10 * MAX_STRING_LEN
#define MAX_NODES 50 /* Max num of nodes (leaves + internal) of a Tree */

struct stack {
    char s[ STACK_SIZE ][ MAX_STRING_LEN ];
    int top;
};
typedef struct stack Stack;
typedef Stack *StackPtr;


void push( StackPtr stackPtr, char * item );
int isStackEmpty( StackPtr stackPtr );
char* pop( StackPtr stackPtr );
void show( StackPtr stackPtr );

#endif