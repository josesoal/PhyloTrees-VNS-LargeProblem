/*
 ============================================================================
 Authors :
 	Jose Luis Soncco-Alvarez and Mauricio Ayala-Rincon
 	Group of Theory of Computation
 	Universidade de Brasilia (UnB) - Brazil
 ============================================================================
 */

#include <stdio.h>
#include "auxiliary.h"

/* auxiliar functions*/ 
void nomemMessage( char *string ) {
	fprintf( stderr, "stderr: no memory avaliable for %s.\n", string );
}

 int fact( int n ) 
 {
 	int i;
	int f = 1;
    for ( i = 1; i <= n; i++ ) {
          f *= i;
    }
    return f;
}

int min( int a, int b )
{
	int min = a;
	if ( b < a )
		min = b;
	return min;
}

/* NOTE: i must be <= j */
void applyReversal( int *genome, int i, int j ) 
{
	int k, temp;
	if ( i == j ) {
		genome[ i ] = -1 * genome[ i ];
	}
	else {
		for ( k = 0; k < (j - i + 1)/2; k++ ) {
			temp = -1 * genome[ i + k ];
			genome[ i + k ] = -1 * genome[ j - k ];
			genome[ j - k ] = temp;
		}	
	}
}