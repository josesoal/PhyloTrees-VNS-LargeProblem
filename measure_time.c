#include "measure_time.h"

/* Return the time difference in SECONDS */
double timeval_diff(struct timeval *a, struct timeval *b)
{
    return
    (double)(a->tv_sec + (double)a->tv_usec/1000000) -
    (double)(b->tv_sec + (double)b->tv_usec/1000000);
}