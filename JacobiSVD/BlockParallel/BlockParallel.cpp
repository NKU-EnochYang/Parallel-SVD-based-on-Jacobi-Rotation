#include "BlockParallel.h"

bool BlockParallel::rotate(double *A, double *V, int N, double threshold, int M, int n)
{
    int *I = new int[(int) ceil(0.5 * (M - 1))];
    int *J = new int[(int) ceil(0.5 * (M - 1))];
    int solved = 0;

#pragma omp parallel for num_threads(NUM_THREADS) reduction(+:solved) schedule(static) shared(A, V, n, threshold, N, M)
    for (int a = 0; a < M; a++)
    {
        for (int i = a * n; i < a * n + min(n, N - a * n) - 1; i++)
        {
            for (int j = i + 1; j < a * n + min(n, N - a * n); j++)
            {
                double alpha = 0;
                double beta = 0;
                double gamma = 0;
                double cosi, sine, q, c;
                for (int h = 0; h < N; h++)
                {
                    alpha += A[h + N * i] * A[h + N * j];
                    beta += A[h + N * i] * A[h + N * i];
                    gamma += A[h + N * j] * A[h + N * j];
                }
                q = beta - gamma;
                if ((alpha * alpha / (beta * gamma) < threshold) || (beta < threshold) || (gamma < threshold))
                    solved++;
                else
                {
                    c = sqrt(4 * alpha * alpha + q * q);
                    if (q >= 0)
                    {
                        cosi = sqrt((c + q) / (2 * c));
                        sine = alpha / (c * cosi);
                    }
                    else
                    {
                        sine = (alpha >= 0) ? sqrt((c - q) / 2 / c) : -sqrt((c - q) / 2 / c);
                        cosi = alpha / (c * sine);
                    }
                    double tA, tV;
                    for (int h = 0; h < N; h++)
                    {
                        tA = A[h + N * i];
                        A[h + N * i] = cosi * A[h + N * i] + sine * A[h + N * j];
                        A[h + N * j] = -sine * tA + cosi * A[h + N * j];
                    }
                    for (int h = 0; h < N; h++)
                    {
                        tV = V[h + N * i];
                        V[h + N * i] = cosi * V[h + N * i] + sine * V[h + N * j];
                        V[h + N * j] = -sine * tV + cosi * V[h + N * j];
                    }
                }

            }
        }
    }

    for (int h = 0; h < M; h++)
    {
        int k = h + 1;
        if (k == M)
        {
            for (int i = 2; i <= (int) ceil(0.5 * M); i++)
            {
                I[i - 2] = i - 1;
                J[i - 2] = M + 1 - i;
            }
        }
        else
        {
            for (int i = 1; i <= (int) ceil(0.5 * (M - k)); i++)
                I[i - 1] = i - 1;
            for (int i = 1; i <= (int) ceil(0.5 * (M - k)); i++)
                J[i - 1] = M - k + 1 - i;
            if (k > 2)
            {
                int j = (int) ceil(0.5 * (M - k));
                for (int i = M - k + 2; i <= M - (int) floor(0.5 * k); i++)
                {
                    I[j] = i - 1;
                    J[j] = 2 * M - k + 1 - i;
                    j++;
                }
            }

        }
        int ng = (k % 2 == 0) ? (int) floor(0.5 * (M - 1)) : (int) floor(0.5 * M);

#pragma omp parallel for num_threads(NUM_THREADS) reduction(+:solved) schedule(static) shared(A, V, I, J, n, threshold, N, ng)
        for (int g = 0; g < ng; g++)
        {
            int block_i = I[g];
            int block_j = J[g];
            for (int i = block_i * n; i < block_i * n + min(n, N - block_i * n); i++)
            {
                for (int j = block_j * n; j < block_j * n + min(n, N - block_j * n); j++)
                {
                    double alpha = 0;
                    double beta = 0;
                    double gamma = 0;
                    double cosi, sine, q, c;
                    int h;
                    for (h = 0; h < N; h++)
                    {
                        alpha += A[h + N * i] * A[h + N * j];
                        beta += A[h + N * i] * A[h + N * i];
                        gamma += A[h + N * j] * A[h + N * j];
                    }
                    q = beta - gamma;
                    if ((alpha * alpha / (beta * gamma) < threshold) || (beta < threshold) || (gamma < threshold))
                        solved++;
                    else
                    {
                        c = sqrt(4 * alpha * alpha + q * q);
                        if (q >= 0)
                        {
                            cosi = sqrt((c + q) / (2 * c));
                            sine = alpha / (c * cosi);
                        }
                        else
                        {
                            sine = (alpha >= 0) ? sqrt((c - q) / 2 / c) : -sqrt((c - q) / 2 / c);
                            cosi = alpha / (c * sine);
                        }
                        double tA, tV;
                        for (h = 0; h < N; h++)
                        {
                            tA = A[h + N * i];
                            A[h + N * i] = cosi * A[h + N * i] + sine * A[h + N * j];
                            A[h + N * j] = -sine * tA + cosi * A[h + N * j];
                        }
                        for (h = 0; h < N; h++)
                        {
                            tV = V[h + N * i];
                            V[h + N * i] = cosi * V[h + N * i] + sine * V[h + N * j];
                            V[h + N * j] = -sine * tV + cosi * V[h + N * j];
                        }
                    }

                }
            }
        }
    }
    return (2 * solved) == (N * (N - 1));
}

void BlockParallel::SVD(double *A, double *U, double *V, int N, double threshold, double *sigma)
{
    int n = 20;
    int M = N / n + int(((N % n) != 0) ? 1 : 0);
    bool converged = false;
    while (!converged)
    {
        converged = rotate(A, V, N, threshold, M, n);
    }
#pragma omp parallel for schedule(static) shared(sigma, A, U, N)
    for (int i = 0; i < N; i++)
    {
        double si = 0;
        for (int j = 0; j < N; j++)
            si += A[j + N * i] * A[j + N * i];
        si = sqrt(si);
        for (int k = 0; k < N; k++)
            U[k + N * i] = A[k + N * i] / si;
        sigma[i] = si;
    }
}
