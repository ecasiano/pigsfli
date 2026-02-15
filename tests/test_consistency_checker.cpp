#include "../src/pimc.hpp"
#include <iostream>
#include <cassert>

using namespace pimc;

int main()
{
    System sys(2, 1, 1.0, 1.0, 0.0);
    Lattice lat(2, 1);
    Configuration config(sys, lat, 1);

    std::vector<int> fock = {1, 0};
    config.initialize(fock);

    Worldline &wl = config.replica(0).worldline();

    wl.insertHop(
        Kink{0.3, 0, 0, 1, -1, -1, 0, 0, -1, 0},
        Kink{0.3, 0, 0, 1, -1, -1, 0, 0, -1, 1});

    wl.checkConsistency(); // should pass

    // break a pointer
    wl[wl.firstKink(0)].next = 9999;

    bool threw = false;
    try
    {
        wl.checkConsistency();
    }
    catch (...)
    {
        threw = true;
    }

    assert(threw);
    std::cout << "test_consistency_checker passed.\n";
}
