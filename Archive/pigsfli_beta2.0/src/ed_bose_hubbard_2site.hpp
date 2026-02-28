#pragma once
#include <vector>
#include <cmath>
#include <algorithm>

namespace pimc
{

    struct EDResult2Site
    {
        double E0;
        std::vector<double> n_avg;
    };

    inline double onsite_energy(int n, double U, double mu)
    {
        return 0.5 * U * n * (n - 1) - mu * n;
    }

    inline EDResult2Site exactDiagonalization2Site(int N, double t, double U, double mu)
    {
        EDResult2Site res;
        res.n_avg.assign(2, 0.0);

        // ============================================================
        // N = 1
        // ============================================================
        if (N == 1)
        {
            double E10 = onsite_energy(1, U, mu) + onsite_energy(0, U, mu);
            double E01 = onsite_energy(0, U, mu) + onsite_energy(1, U, mu);

            double a = E10;
            double d = E01;
            double b = -t;

            double tr = a + d;
            double det = a * d - b * b;
            double disc = std::sqrt(std::max(0.0, tr * tr / 4.0 - det));

            double E_minus = tr / 2.0 - disc;
            res.E0 = E_minus;

            double v0 = b;
            double v1 = E_minus - a;

            double norm = std::sqrt(v0 * v0 + v1 * v1);
            if (norm < 1e-14)
            {
                v0 = v1 = 1.0 / std::sqrt(2.0);
            }
            else
            {
                v0 /= norm;
                v1 /= norm;
            }

            res.n_avg[0] = v0 * v0;
            res.n_avg[1] = v1 * v1;
            return res;
        }

        // ============================================================
        // N = 2
        // ============================================================
        if (N == 2)
        {
            double E20 = onsite_energy(2, U, mu) + onsite_energy(0, U, mu);
            double E11 = onsite_energy(1, U, mu) + onsite_energy(1, U, mu);
            double E02 = onsite_energy(0, U, mu) + onsite_energy(2, U, mu);

            double v = -t * std::sqrt(2.0);

            double H[3][3] = {
                {E20, v, 0.0},
                {v, E11, v},
                {0.0, v, E02}};

            double V[3][3] = {
                {1, 0, 0},
                {0, 1, 0},
                {0, 0, 1}};

            auto jacobi = [&](int p, int q)
            {
                if (H[p][q] == 0.0)
                    return;
                double app = H[p][p];
                double aqq = H[q][q];
                double apq = H[p][q];

                double phi = 0.5 * std::atan2(2.0 * apq, aqq - app);
                double c = std::cos(phi);
                double s = std::sin(phi);

                for (int k = 0; k < 3; ++k)
                {
                    double Apk = H[p][k];
                    double Aqk = H[q][k];
                    H[p][k] = c * Apk - s * Aqk;
                    H[q][k] = s * Apk + c * Aqk;
                }
                for (int k = 0; k < 3; ++k)
                {
                    double Akp = H[k][p];
                    double Akq = H[k][q];
                    H[k][p] = c * Akp - s * Akq;
                    H[k][q] = s * Akp + c * Akq;
                }

                H[p][q] = H[q][p] = 0.0;

                for (int k = 0; k < 3; ++k)
                {
                    double Vkp = V[k][p];
                    double Vkq = V[k][q];
                    V[k][p] = c * Vkp - s * Vkq;
                    V[k][q] = s * Vkp + c * Vkq;
                }
            };

            for (int sweep = 0; sweep < 10; ++sweep)
            {
                jacobi(0, 1);
                jacobi(0, 2);
                jacobi(1, 2);
            }

            double evals[3] = {H[0][0], H[1][1], H[2][2]};
            int idx = 0;
            if (evals[1] < evals[idx])
                idx = 1;
            if (evals[2] < evals[idx])
                idx = 2;

            res.E0 = evals[idx];

            double v0 = V[0][idx];
            double v1 = V[1][idx];
            double v2 = V[2][idx];

            double norm = std::sqrt(v0 * v0 + v1 * v1 + v2 * v2);
            v0 /= norm;
            v1 /= norm;
            v2 /= norm;

            res.n_avg[0] = 2 * v0 * v0 + 1 * v1 * v1 + 0 * v2 * v2;
            res.n_avg[1] = 0 * v0 * v0 + 1 * v1 * v1 + 2 * v2 * v2;
            return res;
        }

        // ============================================================
        // N = 3
        // ============================================================
        if (N == 3)
        {
            double E30 = onsite_energy(3, U, mu) + onsite_energy(0, U, mu);
            double E21 = onsite_energy(2, U, mu) + onsite_energy(1, U, mu);
            double E12 = onsite_energy(1, U, mu) + onsite_energy(2, U, mu);
            double E03 = onsite_energy(0, U, mu) + onsite_energy(3, U, mu);

            double v1 = -t * std::sqrt(3.0);
            double v2 = -t * 2.0;
            double v3 = -t * std::sqrt(3.0);

            double H[4][4] = {
                {E30, v1, 0.0, 0.0},
                {v1, E21, v2, 0.0},
                {0.0, v2, E12, v3},
                {0.0, 0.0, v3, E03}};

            double V[4][4] = {
                {1, 0, 0, 0},
                {0, 1, 0, 0},
                {0, 0, 1, 0},
                {0, 0, 0, 1}};

            auto jacobi4 = [&](int p, int q)
            {
                if (H[p][q] == 0.0)
                    return;
                double app = H[p][p];
                double aqq = H[q][q];
                double apq = H[p][q];

                double phi = 0.5 * std::atan2(2.0 * apq, aqq - app);
                double c = std::cos(phi);
                double s = std::sin(phi);

                for (int k = 0; k < 4; ++k)
                {
                    double Apk = H[p][k];
                    double Aqk = H[q][k];
                    H[p][k] = c * Apk - s * Aqk;
                    H[q][k] = s * Apk + c * Aqk;
                }
                for (int k = 0; k < 4; ++k)
                {
                    double Akp = H[k][p];
                    double Akq = H[k][q];
                    H[k][p] = c * Akp - s * Akq;
                    H[k][q] = s * Akp + c * Akq;
                }

                H[p][q] = H[q][p] = 0.0;

                for (int k = 0; k < 4; ++k)
                {
                    double Vkp = V[k][p];
                    double Vkq = V[k][q];
                    V[k][p] = c * Vkp - s * Vkq;
                    V[k][q] = s * Vkp + c * Vkq;
                }
            };

            for (int sweep = 0; sweep < 15; ++sweep)
            {
                jacobi4(0, 1);
                jacobi4(0, 2);
                jacobi4(0, 3);
                jacobi4(1, 2);
                jacobi4(1, 3);
                jacobi4(2, 3);
            }

            double evals[4] = {H[0][0], H[1][1], H[2][2], H[3][3]};
            int idx = 0;
            for (int k = 1; k < 4; ++k)
                if (evals[k] < evals[idx])
                    idx = k;

            res.E0 = evals[idx];

            double v0 = V[0][idx];
            double v1c = V[1][idx];
            double v2c = V[2][idx];
            double v3c = V[3][idx];

            double norm = std::sqrt(v0 * v0 + v1c * v1c + v2c * v2c + v3c * v3c);
            v0 /= norm;
            v1c /= norm;
            v2c /= norm;
            v3c /= norm;

            res.n_avg[0] = 3 * v0 * v0 + 2 * v1c * v1c + 1 * v2c * v2c + 0 * v3c * v3c;
            res.n_avg[1] = 0 * v0 * v0 + 1 * v1c * v1c + 2 * v2c * v2c + 3 * v3c * v3c;
            return res;
        }

        // fallback
        res.E0 = 0.0;
        res.n_avg[0] = res.n_avg[1] = 0.0;
        return res;
    }

} // namespace pimc
