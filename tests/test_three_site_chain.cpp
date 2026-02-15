#include "../src/pimc.hpp"
#include <cassert>
#include <iostream>

using namespace pimc;

int main()
{
    System sys(3, 1, 1.0, 1.0, 0.0);
    Lattice lat(3, 1);
    Configuration config(sys, lat, 1);

    std::vector<int> fock = {0, 1, 0};
    config.initialize(fock);

    Worldline &wl = config.replica(0).worldline();

    // 1 -> 2 at 0.2
    wl.insertHop(
        Kink{0.2, 0, 1, 2, -1, -1, 0, 0, -1, 1},
        Kink{0.2, 0, 1, 2, -1, -1, 0, 0, -1, 2});

    // 2 -> 1 at 0.4
    wl.insertHop(
        Kink{0.4, 0, 2, 1, -1, -1, 0, 0, -1, 2},
        Kink{0.4, 0, 2, 1, -1, -1, 0, 0, -1, 1});

    // 1 -> 0 at 0.6
    wl.insertHop(
        Kink{0.6, 0, 1, 0, -1, -1, 0, 0, -1, 1},
        Kink{0.6, 0, 1, 0, -1, -1, 0, 0, -1, 0});

    wl.checkConsistency();

    assert(wl.occupationAt(0, 0.1) == 0);
    assert(wl.occupationAt(1, 0.1) == 1);
    assert(wl.occupationAt(2, 0.1) == 0);

    assert(wl.occupationAt(2, 0.3) == 1);

    assert(wl.occupationAt(1, 0.5) == 1);

    assert(wl.occupationAt(0, 0.7) == 1);

    std::cout << "test_three_site_chain passed.\n";
}
