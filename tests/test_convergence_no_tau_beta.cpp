#include "../src/pimc.hpp"
#include "../src/moves_mc.hpp"
#include "../src/estimators.hpp"
#include "../src/ed_bose_hubbard_2site.hpp"

#include <iostream>
#include <vector>
#include <cmath>

using namespace pimc;

int main()
{
    std::cout << "Running test_convergence_no_tau_beta...\n";

    int N = 2;
    int L = 2;
    double t = 0.5;
    double U = 1.0;
    double mu = 0.0;

    EDResult2Site ed = exactDiagonalization2Site(N, t, U, mu);
    double E_exact = ed.E0;

    std::cout << "Exact ground state energy = " << E_exact << "\n\n";
    std::cout << "beta\tE_MC\n";

    std::vector<double> betas = {0.5, 1.0, 2.0, 4.0, 8.0};

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

        RNG rng(54321 + int(beta * 100));

        KinkAntikinkInsertionMC insertMove(params);
        KinkAntikinkRemovalMC removeMove(params);
        BoundaryKinkAntikinkInsertion0MC bInsert0(params);
        BoundaryKinkAntikinkRemoval0MC bRemove0(params);

        TotalEnergyEstimator Etot(beta);

        int nSteps = 2000000;
        int therm = 200000;

        double accumE = 0.0;
        int count = 0;

        for (int step = 0; step < nSteps; ++step)
        {
            double r = rng.uniform();

            if (r < 0.4)
                insertMove.attempt(config, rng);
            else if (r < 0.8)
                removeMove.attempt(config, rng);
            else if (r < 0.9)
                bInsert0.attempt(config, rng);
            else
                bRemove0.attempt(config, rng);

            if (step >= therm)
            {
                accumE += Etot.measure(config);
                count++;
            }
        }

        double E_mc = accumE / std::max(1, count);
        std::cout << beta << "\t" << E_mc << "\n";
    }

    std::cout << "\nDone.\n";
    return 0;
}
