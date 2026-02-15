#include "../src/pimc.hpp"
#include <cassert>
#include <iostream>

using namespace pimc;

int main()
{
    RNG rng(123);

    System sys(2, 2, 1.0, 1.0, 0.0);
    Lattice lat(2, 2);
    Configuration config(sys, lat, 1);

    std::vector<int> fock(sys.size(), 1);
    config.initialize(fock);

    SimulationParameters params(1.0);
    KinkAntikinkInsertion move(params);

    for (int n = 0; n < 20; ++n)
    {
        bool ok = move.attempt(config, rng);
        assert(ok);
    }

    // Check worldline consistency after many insertions
    config.replica(0).worldline().checkConsistency();

    std::cout << "Kink–antikink insertion mechanics test PASSED.\n";
    return 0;
}
