#ifndef MEASURE_TIME_H
#define MEASURE_TIME_H

#include <time.h>
#include <sys/time.h>

double timeval_diff(struct timeval *a, struct timeval *b);

#endif