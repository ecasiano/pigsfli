#include "../src/pimc.hpp"
#include "../src/moves_mc.hpp"
#include "../src/estimators.hpp"

#include <cassert>
#include <iostream>
#include <cmath>

using namespace pimc;

int main()
{
    std::cout << "Running test_mc_diagonal_energy...\n";

    double beta = 1.0;

    // 2-site system, U = 0, mu = 0 => E_onsite(n) = 0 for all n
    // So exact diagonal energy = 0, regardless of configuration.
    System sys(2, 1, 0.5, 0.0, 0.0); // t = 0.5, U = 0, mu = 0
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

    int nSteps = 20000;
    int therm = 2000;

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
            double e = Ediag.measure(config);
            accumE += e;
            count++;
        }
    }

    double avgE = accumE / std::max(1, count);

    std::cout << "Average diagonal energy (U=0, mu=0) = " << avgE << "\n";

    // Should be close to 0
    assert(std::abs(avgE) < 1e-2);

    std::cout << "test_mc_diagonal_energy passed.\n";
    return 0;
}
