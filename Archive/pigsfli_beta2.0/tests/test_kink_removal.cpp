#include "../src/pimc.hpp"
#include <cassert>
#include <iostream>

using namespace pimc;

int main()
{
    std::cout << "Running test_kink_removal...\n";

    // Simple 2-site system
    System sys(2, 1, 1.0, 1.0, 0.0);
    Lattice lat(2, 1);
    Configuration config(sys, lat, 1);

    // Initial Fock state: 1 boson on site 0
    std::vector<int> fock = {1, 0};
    config.initialize(fock);

    Worldline &wl = config.replica(0).worldline();

    // Insert hop 0 -> 1 at tau = 0.4
    Kink k1{0.4, 0, 0, 1, -1, -1, 0, 0, -1, 0}; // lives on site 0 list
    Kink k2{0.4, 0, 0, 1, -1, -1, 0, 0, -1, 1}; // lives on site 1 list

    auto [idxA, idxB] = wl.insertHop(k1, k2);

    wl.checkConsistency();

    // After insertion:
    // tau < 0.4 → (1,0)
    // tau > 0.4 → (0,1)
    assert(wl.occupationAt(0, 0.2) == 1);
    assert(wl.occupationAt(1, 0.2) == 0);

    assert(wl.occupationAt(0, 0.5) == 0);
    assert(wl.occupationAt(1, 0.5) == 1);

    // Now remove the hop
    wl.deleteHop(idxA);

    wl.checkConsistency();

    // After removal, worldline must return to original state
    assert(wl.occupationAt(0, 0.2) == 1);
    assert(wl.occupationAt(1, 0.2) == 0);

    assert(wl.occupationAt(0, 0.8) == 1);
    assert(wl.occupationAt(1, 0.8) == 0);

    std::cout << "test_kink_removal passed.\n";
    return 0;
}
