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

    // 0 -> 1 at 0.3
    wl.insertHop(
        Kink{0.3, 0, 0, 1, -1, -1, 0, 0, -1, 0},
        Kink{0.3, 0, 0, 1, -1, -1, 0, 0, -1, 1});

    // 1 -> 0 at 0.7
    wl.insertHop(
        Kink{0.7, 0, 1, 0, -1, -1, 0, 0, -1, 1},
        Kink{0.7, 0, 1, 0, -1, -1, 0, 0, -1, 0});

    wl.checkConsistency();

    assert(wl.occupationAt(0, 0.2) == 1);
    assert(wl.occupationAt(1, 0.2) == 0);

    assert(wl.occupationAt(0, 0.5) == 0);
    assert(wl.occupationAt(1, 0.5) == 1);

    assert(wl.occupationAt(0, 0.9) == 1);
    assert(wl.occupationAt(1, 0.9) == 0);

    std::cout << "test_round_trip passed.\n";
}
