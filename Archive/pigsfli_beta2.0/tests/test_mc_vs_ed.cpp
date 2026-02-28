#include "../src/pimc.hpp"
#include "../src/moves_mc.hpp"
#include "../src/estimators.hpp"
#include "../src/ed_bose_hubbard_2site.hpp"

#include <iostream>
#include <cassert>
#include <cmath>

using namespace pimc;

int main()
{
    std::cout << "Running test_mc_vs_ed...\n";

    // Physical parameters
    double beta = 4.0; // large enough to project onto ground state
    double t = 0.5;
    double U = 1.0;
    double mu = 0.0;
    int N = 1; // total boson number

    // Exact diagonalization
    EDResult2Site ed = exactDiagonalization2Site(N, t, U, mu);
    double E_exact = ed.E0;

    std::cout << "Exact ground state energy = " << E_exact << "\n";

    // Build MC system
    System sys(2, N, t, U, mu);
    Lattice lat(2, 1);
    Configuration config(sys, lat, 1);

    BoseHubbardHamiltonian H(sys, lat);
    config.setHamiltonian(&H);

    // Initial Fock state: 1 boson on site 0
    std::vector<int> fock = {1, 0};
    config.initialize(fock);

    RNG rng(12345);

    // Moves
    KinkAntikinkInsertionMC insertMove(beta);
    KinkAntikinkRemovalMC removeMove(beta);

    // Estimator
    DiagonalEnergyEstimator Ediag(beta);
    KineticEnergyEstimator Ekin(beta);
    TotalEnergyEstimator Etot(beta);

    int nSteps = 30000;
    int therm = 5000;

    double accumE = 0.0;
    int count = 0;

    for (int step = 0; step < nSteps; ++step)
    {
        // Simple move selection
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

    std::cout << "MC diagonal energy = " << E_mc << "\n";

    // Exact diagonal energy for N=1, U=0, mu=0 is 0
    double E_exact_diag = 0.0;

    double diff = std::abs(E_mc - E_exact);
    std::cout << "Difference = " << diff << "\n";
    assert(diff < 0.05);

    std::cout << "test_mc_vs_ed passed.\n";
    return 0;
}
