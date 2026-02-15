#include "../src/pimc.hpp"
#include <cassert>
#include <iostream>

using namespace pimc;

int main()
{
    System sys(4, 1, 1.0, 1.0, 0.0);
    Lattice lat(4, 1);
    Configuration config(sys, lat, 1);

    std::vector<int> fock = {2, 0, 1, 3};
    config.initialize(fock);

    Worldline &wl = config.replica(0).worldline();

    for (int site = 0; site < sys.size(); site++)
    {
        for (double tau : {0.0, 0.25, 0.5, 0.75, 0.99})
        {
            assert(wl.occupationAt(site, tau) == fock[site]);
        }
    }

    wl.checkConsistency();
    std::cout << "test_fock_init passed.\n";
}
