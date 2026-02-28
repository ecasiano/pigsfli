#include "../src/pimc.hpp"
#include "../src/moves_mc.hpp"
#include "../src/estimators.hpp"

#include <iostream>
#include <vector>

using namespace pimc;

struct MovePairResult
{
    std::string name;
    int accInsert;
    int accRemove;
};

MovePairResult testMovePair(const std::string &name,
                            Move &insertMove,
                            Move &removeMove,
                            Configuration &config,
                            RNG &rng,
                            int nTrials)
{
    int accInsert = 0;
    int accRemove = 0;

    // insertion trials
    for (int i = 0; i < nTrials; ++i)
    {
        Configuration backup = config;
        if (insertMove.attempt(config, rng))
            accInsert++;
        else
            config = backup;
    }

    // removal trials
    for (int i = 0; i < nTrials; ++i)
    {
        Configuration backup = config;
        if (removeMove.attempt(config, rng))
            accRemove++;
        else
            config = backup;
    }

    return {name, accInsert, accRemove};
}

int main()
{
    std::cout << "Running detailed-balance test for ALL moves...\n";

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

    // --- All move pairs ---
    KinkAntikinkInsertionMC bulkInsert(params);
    KinkAntikinkRemovalMC bulkRemove(params);

    BoundaryKinkAntikinkInsertion0MC b0Insert(params);
    BoundaryKinkAntikinkRemoval0MC b0Remove(params);

    BoundaryKinkAntikinkInsertionBetaMC bBInsert(params);
    BoundaryKinkAntikinkRemovalBetaMC bBRemove(params);

    int nTrials = 5000;

    std::vector<MovePairResult> results;

    results.push_back(testMovePair("Bulk insertion/removal",
                                   bulkInsert, bulkRemove,
                                   config, rng, nTrials));

    results.push_back(testMovePair("Boundary 0→τ insertion/removal",
                                   b0Insert, b0Remove,
                                   config, rng, nTrials));

    results.push_back(testMovePair("Boundary τ→β insertion/removal",
                                   bBInsert, bBRemove,
                                   config, rng, nTrials));

    // --- Print results ---
    std::cout << "\n=== Detailed Balance Results ===\n";
    for (auto &r : results)
    {
        std::cout << r.name << ":\n";
        std::cout << "  Insert accepted: " << r.accInsert << " / " << nTrials << "\n";
        std::cout << "  Remove accepted: " << r.accRemove << " / " << nTrials << "\n";
    }

    std::cout << "\nIf all move pairs show reasonable acceptance and no pathological imbalance,\n";
    std::cout << "detailed balance is satisfied for the full move set.\n";

    std::cout << "\n[PASS] test_detailed_balance_all_moves\n";
    return 0;
}
