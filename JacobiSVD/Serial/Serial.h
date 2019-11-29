#ifndef JACOBISVD_SERIAL_H
#define JACOBISVD_SERIAL_H

#include "../common.h"

class Serial
{
private:
    bool rotate(double *A, double *V, int N, double threshold);

public:
    void SVD(double *A, double *U, double *V, int N, double threshold, double *sigma);
};


#endif //JACOBISVD_SERIAL_H
