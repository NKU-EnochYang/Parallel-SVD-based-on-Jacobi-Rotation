#include "Parallel.h"

bool Parallel::rotate(double *A, double *V, int N, double threshold)
{
    int *I = new int[(int) ceil(0.5 * (N - 1))];
    int *J = new int[(int) ceil(0.5 * (N - 1))];
    int solved = 0;
    for (int r = 0; r < N; r++)
    {
        int k = r + 1;
        if (k == N)
        {
            for (int i = 2; i <= (int) ceil(0.5 * N); i++)
            {
                I[i - 2] = i - 1;
                J[i - 2] = N + 1 - i;
            }
        }
        else
        {
            for (int i = 1; i <= (int) ceil(0.5 * (N - k)); i++)
                I[i - 1] = i - 1;
            for (int i = 1; i <= (int) ceil(0.5 * (N - k)); i++)
                J[i - 1] = N - k + 1 - i;
            if (k > 2)
            {
                int j = (int) ceil(0.5 * (N - k));
                for (int i = N - k + 2; i <= N - (int) floor(0.5 * k); i++)
                {
                    I[j] = i - 1;
                    J[j] = 2 * N - k + 1 - i;
                    j++;
                }
            }
        }
        int n = (k % 2 == 0) ? (int) floor(0.5 * (N - 1)) : (int) floor(0.5 * N);

#pragma omp parallel for schedule(static) reduction(+:solved) default(none) shared(I, J, A, V, N, n, threshold)
        for (int g = 0; g < n; g++)
        {
            int i = I[g];
            int j = J[g];

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
            if (alpha * alpha / (beta * gamma) < threshold)
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
    return (2 * solved) == (N * (N - 1));
}

void Parallel::SVD(double *A, double *U, double *V, int N, double threshold, double *sigma)
{
    bool converged = false;
    while (!converged)
    {
        converged = rotate(A, V, N, threshold);
    }
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
