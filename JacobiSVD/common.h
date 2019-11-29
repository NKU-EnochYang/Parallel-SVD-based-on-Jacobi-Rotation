#ifndef JACOBISVD_COMMON_H
#define JACOBISVD_COMMON_H

#include <stdlib.h>
#include <iostream>
#include <sys/time.h>
#include <stdlib.h>
#include <omp.h>
#include <cmath>
#include <algorithm>

using namespace std;

static const int N = 180;
static const double threshold = 1e-5;
static const int NUM_THREADS = 4;

extern double test[N*N];

#endif //JACOBISVD_COMMON_H
