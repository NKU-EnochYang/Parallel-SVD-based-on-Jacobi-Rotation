#ifndef JACOBISVD_BLOCKPARALLEL_H
#define JACOBISVD_BLOCKPARALLEL_H

#include "../common.h"

class BlockParallel
{
private:
    bool rotate(double *A, double *V, int N, double threshold, int M, int n);

public:
    void SVD(double *A, double *U, double *V, int N, double threshold, double *sigma);
};


#endif //JACOBISVD_BLOCKPARALLEL_H
