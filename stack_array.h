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

#define STACK_SIZE 100
#define BUFFER_SIZE 50

struct stack {
    char s[STACK_SIZE][BUFFER_SIZE];
    int top;
};
typedef struct stack Stack;
typedef Stack *StackPtr;


void push( StackPtr stackPtr, char * item );
int isEmpty( StackPtr stackPtr );
char* pop( StackPtr stackPtr );
void show( StackPtr stackPtr );

#endif