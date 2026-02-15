#include "../src/pimc.hpp"
#include <cassert>
#include <iostream>

using namespace pimc;

int main()
{
    System sys(2, 1, 1.0, 1.0, 0.0);
    Lattice lat(2, 1);
    Configuration config(sys, lat, 1);

    std::vector<int> fock = {1, 0};
    config.initialize(fock);

    Worldline &wl = config.replica(0).worldline();

    // hop 0 -> 1 at tau=0.4
    Kink k1{0.4, 0, 0, 1, -1, -1, 0, 0, -1, 0};
    Kink k2{0.4, 0, 0, 1, -1, -1, 0, 0, -1, 1};
    wl.insertHop(k1, k2);

    wl.checkConsistency();

    assert(wl.occupationAt(0, 0.3) == 1);
    assert(wl.occupationAt(1, 0.3) == 0);

    assert(wl.occupationAt(0, 0.5) == 0);
    assert(wl.occupationAt(1, 0.5) == 1);

    std::cout << "test_single_hop passed.\n";
}
