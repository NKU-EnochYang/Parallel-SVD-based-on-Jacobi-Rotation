#include "common.h"
#include "BlockParallel/BlockParallel.h"
#include "Serial/Serial.h"
#include "Parallel/Parallel.h"
#include "Matrix/Matrix.h"

int main()
{
    timeval t_start, t_end;
    Matrix::random_vec(test);

    double A[N*N];
    double V[N*N];
    double U[N*N];
    double sigma[N];

    Matrix::reset_vec(test, A);
    Matrix::reset_vec(A, U);
    Matrix::identity_vec(V);

    Serial serial;
    gettimeofday(&t_start, NULL);
    serial.SVD(A, U, V, N, threshold, sigma);
    gettimeofday(&t_end, NULL);

    cout << "Serial SVD time cost: "
         << 1000 * (t_end.tv_sec - t_start.tv_sec) +
            0.001 * (t_end.tv_usec - t_start.tv_usec) << "ms" << endl;

    Matrix::reset_vec(test, A);
    Matrix::reset_vec(A, U);
    Matrix::identity_vec(V);

    Parallel parallel;
    gettimeofday(&t_start, NULL);
    parallel.SVD(A, U, V, N, threshold, sigma);
    gettimeofday(&t_end, NULL);

    cout << "Parallel SVD time cost: "
         << 1000 * (t_end.tv_sec - t_start.tv_sec) +
            0.001 * (t_end.tv_usec - t_start.tv_usec) << "ms" << endl;

    Matrix::reset_vec(test, A);
    Matrix::reset_vec(A, U);
    Matrix::identity_vec(V);

    BlockParallel block_parallel;
    gettimeofday(&t_start, NULL);
    block_parallel.SVD(A, U, V, N, threshold, sigma);
    gettimeofday(&t_end, NULL);

    cout << "Block Parallel SVD time cost: "
         << 1000 * (t_end.tv_sec - t_start.tv_sec) +
            0.001 * (t_end.tv_usec - t_start.tv_usec) << "ms" << endl;
    return 0;
}