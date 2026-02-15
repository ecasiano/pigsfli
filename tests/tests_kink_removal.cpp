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
    KinkAntikinkInsertion ins(params);
    KinkAntikinkRemoval rem(params);

    // Insert a bunch of hops
    for (int n = 0; n < 50; ++n)
    {
        bool ok = ins.attempt(config, rng);
        assert(ok);
    }

    Worldline &wl = config.replica(0).worldline();
    wl.checkConsistency();

    // Now try to remove until no pairs remain
    for (int n = 0; n < 100; ++n)
    {
        rem.attempt(config, rng);
        wl.checkConsistency();
    }

    // Verify no paired kinks remain
    int M = config.latticeSize();
    for (int site = 0; site < M; ++site)
    {
        int idx = wl.firstKink(site);
        while (idx != -1)
        {
            const Kink &k = wl[idx];
            assert(k.partner == -1); // only flat or unpaired kinks
            idx = k.next;
        }
    }

    std::cout << "Kink–antikink removal mechanics test PASSED.\n";
    return 0;
}
