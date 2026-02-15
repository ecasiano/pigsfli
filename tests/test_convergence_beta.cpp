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
    std::cout << "Running test_convergence_beta...\n";

    double t = 0.5;
    double U = 1.0;
    double mu = 0.0;
    int N = 2; // try N=1,2,3 as you like

    EDResult2Site ed = exactDiagonalization2Site(N, t, U, mu);
    double E_exact = ed.E0;

    std::cout << "Exact ground state energy (N=" << N << ") = " << E_exact << "\n\n";
    std::cout << "beta\tE_MC\n";

    std::vector<double> betas = {0.5, 1.0, 2.0, 4.0, 8.0};

    for (double beta : betas)
    {
        System sys(2, N, t, U, mu);
        Lattice lat(2, 1);
        Configuration config(sys, lat, 1);

        BoseHubbardHamiltonian H(sys, lat);
        config.setHamiltonian(&H);

        // simple initial Fock state: all bosons on site 0
        std::vector<int> fock = {N, 0};
        config.initialize(fock);

        RNG rng(12345 + int(beta * 100));

        KinkAntikinkInsertionMC insertMove(beta);
        KinkAntikinkRemovalMC removeMove(beta);

        TotalEnergyEstimator Etot(beta);

        int nSteps = 40000;
        int therm = 5000;

        double accumE = 0.0;
        int count = 0;

        for (int step = 0; step < nSteps; ++step)
        {
            if (rng.uniform() < 0.5)
                insertMove.attempt(config, rng);
            else
                removeMove.attempt(config, rng);

            if (step >= therm)
            {
                double e = Etot.measure(config);
                accumE += e;
                count++;
            }
        }

        double E_mc = accumE / std::max(1, count);
        std::cout << beta << "\t" << E_mc << "\n";
    }

    std::cout << "\nDone.\n";
    return 0;
}
