if (N == 3)
{
    // Basis: |3,0>, |2,1>, |1,2>, |0,3>
    // Diagonal:
    double E30 = onsite_energy(3, U, mu) + onsite_energy(0, U, mu);
    double E21 = onsite_energy(2, U, mu) + onsite_energy(1, U, mu);
    double E12 = onsite_energy(1, U, mu) + onsite_energy(2, U, mu);
    double E03 = onsite_energy(0, U, mu) + onsite_energy(3, U, mu);

    // Off-diagonal hopping:
    // <3,0|H|2,1> = -t * sqrt(3)
    // <2,1|H|1,2> = -t * 2
    // <1,2|H|0,3> = -t * sqrt(3)
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

    auto jacobi_rotate4 = [](double A[4][4], double V[4][4], int p, int q)
    {
        if (A[p][q] == 0.0)
            return;
        double app = A[p][p];
        double aqq = A[q][q];
        double apq = A[p][q];

        double phi = 0.5 * std::atan2(2.0 * apq, aqq - app);
        double c = std::cos(phi);
        double s = std::sin(phi);

        for (int k = 0; k < 4; ++k)
        {
            double Apk = A[p][k];
            double Aqk = A[q][k];
            A[p][k] = c * Apk - s * Aqk;
            A[q][k] = s * Apk + c * Aqk;
        }
        for (int k = 0; k < 4; ++k)
        {
            double Akp = A[k][p];
            double Akq = A[k][q];
            A[k][p] = c * Akp - s * Akq;
            A[k][q] = s * Akp + c * Akq;
        }
        A[p][p] = c * c * app - 2.0 * s * c * apq + s * s * aqq;
        A[q][q] = s * s * app + 2.0 * s * c * apq + c * c * aqq;
        A[p][q] = A[q][p] = 0.0;

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
        jacobi_rotate4(H, V, 0, 1);
        jacobi_rotate4(H, V, 0, 2);
        jacobi_rotate4(H, V, 0, 3);
        jacobi_rotate4(H, V, 1, 2);
        jacobi_rotate4(H, V, 1, 3);
        jacobi_rotate4(H, V, 2, 3);
    }

    double evals[4] = {H[0][0], H[1][1], H[2][2], H[3][3]};
    int idx_min = 0;
    for (int k = 1; k < 4; ++k)
        if (evals[k] < evals[idx_min])
            idx_min = k;

    res.E0 = evals[idx_min];

    double v0 = V[0][idx_min];
    double v1c = V[1][idx_min];
    double v2c = V[2][idx_min];
    double v3c = V[3][idx_min];

    double norm = std::sqrt(v0 * v0 + v1c * v1c + v2c * v2c + v3c * v3c);
    if (norm < 1e-14)
    {
        v0 = 1.0;
        v1c = v2c = v3c = 0.0;
        norm = 1.0;
    }
    v0 /= norm;
    v1c /= norm;
    v2c /= norm;
    v3c /= norm;

    // <n0> = 3|v0|^2 + 2|v1|^2 + 1|v2|^2 + 0|v3|^2
    // <n1> = 0|v0|^2 + 1|v1|^2 + 2|v2|^2 + 3|v3|^2
    res.n_avg[0] = 3.0 * v0 * v0 + 2.0 * v1c * v1c + 1.0 * v2c * v2c + 0.0 * v3c * v3c;
    res.n_avg[1] = 0.0 * v0 * v0 + 1.0 * v1c * v1c + 2.0 * v2c * v2c + 3.0 * v3c * v3c;

    return res;
}
