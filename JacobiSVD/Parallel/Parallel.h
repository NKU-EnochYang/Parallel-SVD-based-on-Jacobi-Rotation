#ifndef JACOBISVD_PARALLEL_H
#define JACOBISVD_PARALLEL_H

#include "../common.h"

class Parallel
{
private:
    bool rotate(double *A, double *V, int N, double threshold);

public:
    void SVD(double *A, double *U, double *V, int N, double threshold, double *sigma);
};


#endif //JACOBISVD_PARALLEL_H
