#ifndef JACOBISVD_MATRIX_H
#define JACOBISVD_MATRIX_H

#include "../common.h"

class Matrix
{
public:
    static void print_vec(double vec[]);
    static void random_vec(double vec[]);
    static void reset_vec(double src[], double dst[]);
    static void identity_vec(double vec[]);
};


#endif //JACOBISVD_MATRIX_H
