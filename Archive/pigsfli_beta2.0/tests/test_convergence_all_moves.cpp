#include "../src/pimc.hpp"
#include "../src/moves_mc.hpp"
#include "../src/estimators.hpp"
#include "../src/ed_bose_hubbard_2site.hpp"

#include <iostream>
#include <vector>
#include <cmath>

using namespace pimc;

// ============================================================
// Convergence test using both estimators
// ============================================================

int main()
{
    std::cout << "Running test_convergence_midtime...\n";

    int N = 2;
    int L = 2;
    double t = 0.5;
    double U = 1.0;
    double mu = 0.0;

    EDResult2Site ed = exactDiagonalization2Site(N, t, U, mu);
    double E_exact = ed.E0;

    std::cout << "Exact ground state energy = " << E_exact << "\n\n";
    std::cout << "beta\tE_total\t\tE_mid\n";

    std::vector<double> betas = {0.5, 1.0, 2.0, 4.0, 8.0, 16.0};

    for (double beta : betas)
    {

        SimulationParameters params(beta);

        System sys(L, N, t, U, mu);
        Lattice lat(L, 1);
        Configuration config(sys, lat, 1);

        BoseHubbardHamiltonian H(sys, lat);
        config.setHamiltonian(&H);

        std::vector<int> fock = {N, 0};
        config.initialize(fock);

        RNG rng(12345 + int(beta * 100));

        // Moves
        KinkAntikinkInsertionMC insertMove(params);
        KinkAntikinkRemovalMC removeMove(params);
        BoundaryKinkAntikinkInsertion0MC bInsert0(params);
        BoundaryKinkAntikinkRemoval0MC bRemove0(params);
        BoundaryKinkAntikinkInsertionBetaMC bInsertB(params);
        BoundaryKinkAntikinkRemovalBetaMC bRemoveB(params);

        // Estimators
        TotalEnergyEstimator Etot(beta);
        MidTimeEnergyEstimator Emid(beta);

        int nSteps = 2000000;
        int therm = 200000;

        double accumE_tot = 0.0;
        double accumE_mid = 0.0;
        int count = 0;

        for (int step = 0; step < nSteps; ++step)
        {
            double r = rng.uniform();

            if (r < 0.3)
                insertMove.attempt(config, rng);
            else if (r < 0.6)
                removeMove.attempt(config, rng);
            else if (r < 0.75)
                bInsert0.attempt(config, rng);
            else if (r < 0.9)
                bRemove0.attempt(config, rng);
            else if (r < 0.95)
                bInsertB.attempt(config, rng);
            else
                bRemoveB.attempt(config, rng);

            if (step >= therm)
            {
                accumE_tot += Etot.measure(config);
                accumE_mid += Emid.measure(config);
                count++;
            }
        }

        double E_mc_tot = accumE_tot / std::max(1, count);
        double E_mc_mid = accumE_mid / std::max(1, count);

        std::cout << beta << "\t" << E_mc_tot << "\t" << E_mc_mid << "\n";
    }

    std::cout << "\nDone.\n";
    return 0;
}
