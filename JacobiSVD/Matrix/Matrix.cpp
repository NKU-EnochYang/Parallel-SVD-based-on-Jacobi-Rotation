#include "Matrix.h"

double test[N * N] = {0};

void Matrix::random_vec(double *vec)
{
    srand((unsigned) time(NULL));
    for (int i = 0; i < N * N; i++)
        vec[i] = rand() / 100;
}

void Matrix::reset_vec(double *src, double *dst)
{
    for (int i = 0; i < N * N; i++)
        dst[i] = src[i];
}

void Matrix::print_vec(double vec[])
{
    if (N > 16)
        return;
    cout << endl;
    for (int i = 0; i < N; i++)
        cout << vec[i] << " ";
    cout << endl;
}

void Matrix::identity_vec(double *vec)
{
    for (int i = 0; i < N; i++)
    {
        for (int j = 0; j < N; j++)
        {
            if (i == j)
                vec[i * N + j] = 1;
            else
                vec[i * N + j] = 0;
        }
    }
}