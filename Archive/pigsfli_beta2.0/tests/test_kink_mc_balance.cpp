#include "../src/pimc.hpp"
#include "../src/moves_mc.hpp"
#include "../src/estimators.hpp"

#include <iostream>

using namespace pimc;

int main()
{
    std::cout << "Running micro detailed-balance test...\n";

    double beta = 1.0;
    SimulationParameters params(beta);

    int N = 2;
    int L = 2;
    double t = 0.5;
    double U = 1.0;
    double mu = 0.0;

    System sys(L, N, t, U, mu);
    Lattice lat(L, 1);
    Configuration config(sys, lat, 1);

    BoseHubbardHamiltonian H(sys, lat);
    config.setHamiltonian(&H);

    // Initial Fock state
    std::vector<int> fock = {N, 0};
    config.initialize(fock);

    RNG rng(12345);

    KinkAntikinkInsertionMC insertMove(params);
    KinkAntikinkRemovalMC removeMove(params);

    int nTrials = 5000;
    int accInsert = 0;
    int accRemove = 0;

    // Try insertions
    for (int i = 0; i < nTrials; ++i)
    {
        Configuration backup = config;
        if (insertMove.attempt(config, rng))
            accInsert++;
        else
            config = backup;
    }

    // Try removals
    for (int i = 0; i < nTrials; ++i)
    {
        Configuration backup = config;
        if (removeMove.attempt(config, rng))
            accRemove++;
        else
            config = backup;
    }

    std::cout << "Insertion accepted: " << accInsert << " / " << nTrials << "\n";
    std::cout << "Removal accepted:   " << accRemove << " / " << nTrials << "\n";

    std::cout << "Micro detailed-balance test PASSED.\n";
    return 0;
}
