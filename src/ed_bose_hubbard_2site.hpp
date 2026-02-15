#pragma once
#include <vector>
#include <cmath>
#include <algorithm>

namespace pimc
{

    // Simple struct to hold ED results
    struct EDResult2Site
    {
        double E0;                 // ground state energy
        std::vector<double> n_avg; // <n_i> on each site
    };

    // Onsite energy: U/2 n(n-1) - mu n
    inline double onsite_energy(int n, double U, double mu)
    {
        return 0.5 * U * n * (n - 1) - mu * n;
    }

    // ED for 2-site Bose-Hubbard with fixed total N = 1 or 2
    inline EDResult2Site exactDiagonalization2Site(int N, double t, double U, double mu)
    {
        EDResult2Site res;
        res.n_avg.assign(2, 0.0);

        if (N == 1)
        {
            // Basis: |1,0>, |0,1>
            // H = [ E10   -t
            //       -t    E01 ]
            double E10 = onsite_energy(1, U, mu) + onsite_energy(0, U, mu);
            double E01 = onsite_energy(0, U, mu) + onsite_energy(1, U, mu);

            double a = E10;
            double d = E01;
            double b = -t;

            double tr = a + d;
            double det = a * d - b * b;
            double disc = std::sqrt(std::max(0.0, tr * tr / 4.0 - det));

            double E_minus = tr / 2.0 - disc;
            double E_plus = tr / 2.0 + disc;

            res.E0 = E_minus;

            // Ground state vector components (normalized)
            // Solve (H - E0 I) v = 0
            double v0 = b;
            double v1 = E_minus - a; // correct ground-state eigenvalue

            double norm = std::sqrt(v0 * v0 + v1 * v1);
            if (norm < 1e-14)
            {
                // symmetric case E10 = E01
                v0 = 1.0 / std::sqrt(2.0);
                v1 = 1.0 / std::sqrt(2.0);
            }
            else
            {
                v0 /= norm;
                v1 /= norm;
            }

            // <n0> = |v0|^2 * 1 + |v1|^2 * 0
            // <n1> = |v0|^2 * 0 + |v1|^2 * 1
            res.n_avg[0] = v0 * v0;
            res.n_avg[1] = v1 * v1;

            return res;
        }

        if (N == 2)
        {
            // Basis: |2,0>, |1,1>, |0,2>
            // Diagonal:
            double E20 = onsite_energy(2, U, mu) + onsite_energy(0, U, mu);
            double E11 = onsite_energy(1, U, mu) + onsite_energy(1, U, mu);
            double E02 = onsite_energy(0, U, mu) + onsite_energy(2, U, mu);

            // Off-diagonal hopping:
            // <2,0|H|1,1> = -t * sqrt(2)
            // <1,1|H|0,2> = -t * sqrt(2)
            double v = -t * std::sqrt(2.0);

            // 3x3 symmetric matrix:
            // [ E20   v    0 ]
            // [  v   E11   v ]
            // [  0    v   E02 ]
            double a = E20, d = E11, f = E02;
            double b = v, c = 0.0, e = v;

            // We can just do a simple numeric diagonalization for 3x3.
            // For simplicity, use a few Jacobi sweeps (good enough here).

            double H[3][3] = {
                {a, b, c},
                {b, d, e},
                {c, e, f}};

            double V[3][3] = {
                {1, 0, 0},
                {0, 1, 0},
                {0, 0, 1}};

            auto jacobi_rotate = [](double A[3][3], double V[3][3], int p, int q)
            {
                if (A[p][q] == 0.0)
                    return;
                double app = A[p][p];
                double aqq = A[q][q];
                double apq = A[p][q];

                double phi = 0.5 * std::atan2(2.0 * apq, aqq - app);
                double c = std::cos(phi);
                double s = std::sin(phi);

                for (int k = 0; k < 3; ++k)
                {
                    double Apk = A[p][k];
                    double Aqk = A[q][k];
                    A[p][k] = c * Apk - s * Aqk;
                    A[q][k] = s * Apk + c * Aqk;
                }
                for (int k = 0; k < 3; ++k)
                {
                    double Akp = A[k][p];
                    double Akq = A[k][q];
                    A[k][p] = c * Akp - s * Akq;
                    A[k][q] = s * Akp + c * Akq;
                }
                A[p][p] = c * c * app - 2.0 * s * c * apq + s * s * aqq;
                A[q][q] = s * s * app + 2.0 * s * c * apq + c * c * aqq;
                A[p][q] = A[q][p] = 0.0;

                for (int k = 0; k < 3; ++k)
                {
                    double Vkp = V[k][p];
                    double Vkq = V[k][q];
                    V[k][p] = c * Vkp - s * Vkq;
                    V[k][q] = s * Vkp + c * Vkq;
                }
            };

            // A few sweeps
            for (int sweep = 0; sweep < 10; ++sweep)
            {
                jacobi_rotate(H, V, 0, 1);
                jacobi_rotate(H, V, 0, 2);
                jacobi_rotate(H, V, 1, 2);
            }

            // Eigenvalues are on the diagonal of H
            double evals[3] = {H[0][0], H[1][1], H[2][2]};
            int idx_min = 0;
            if (evals[1] < evals[idx_min])
                idx_min = 1;
            if (evals[2] < evals[idx_min])
                idx_min = 2;

            res.E0 = evals[idx_min];

            // Ground state vector is column idx_min of V
            double v0 = V[0][idx_min];
            double v1 = V[1][idx_min];
            double v2 = V[2][idx_min];

            double norm = std::sqrt(v0 * v0 + v1 * v1 + v2 * v2);
            if (norm < 1e-14)
            {
                v0 = 1.0;
                v1 = 0.0;
                v2 = 0.0;
                norm = 1.0;
            }
            v0 /= norm;
            v1 /= norm;
            v2 /= norm;

            // <n0> = |v0|^2 * 2 + |v1|^2 * 1 + |v2|^2 * 0
            // <n1> = |v0|^2 * 0 + |v1|^2 * 1 + |v2|^2 * 2
            res.n_avg[0] = v0 * v0 * 2.0 + v1 * v1 * 1.0 + v2 * v2 * 0.0;
            res.n_avg[1] = v0 * v0 * 0.0 + v1 * v1 * 1.0 + v2 * v2 * 2.0;

            return res;
        }

        // Fallback: unsupported N
        res.E0 = 0.0;
        res.n_avg[0] = res.n_avg[1] = 0.0;
        return res;
    }

} // namespace pimc
