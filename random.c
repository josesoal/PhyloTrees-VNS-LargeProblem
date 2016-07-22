/* taken from http://www.azillionmonkeys.com/qed/random.html */

#include <stdio.h>
#include <stdlib.h>
#include "random.h"

#define RS_SCALE (1.0 / (1.0 + RAND_MAX))

/* return double number in the range [0,1) */
double drand (void) 
{
    double d;
    do {
       d = (((rand () * RS_SCALE) + rand ()) * RS_SCALE + rand ()) * RS_SCALE;
    } while (d >= 1); /* Round off */
    return d;
}

/* return integer number in the range [0, x) */
/* Ex: irand(2) returns intergers 0 or 1 */
unsigned int irand(unsigned int x)
{
	return (unsigned int) ((x) * drand ());	
}