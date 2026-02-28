#include "../src/pimc.hpp"
#include <iostream>
#include <cassert>

using namespace pimc;

void dump_worldline(const Worldline &wl)
{
    int M = wl.latticeSize();
    std::cout << "Worldline dump:\n";
    for (int site = 0; site < M; ++site)
    {
        std::cout << "  site " << site << ":\n";
        int idx = wl.firstKink(site);
        while (idx != -1)
        {
            const Kink &k = wl[idx];
            std::cout << "    idx=" << idx
                      << " tau=" << k.tau
                      << " n=" << k.n
                      << " src=" << k.src
                      << " dest=" << k.dest
                      << " site=" << k.site
                      << " partner=" << k.partner
                      << " prev=" << k.prev
                      << " next=" << k.next
                      << "\n";
            idx = k.next;
        }
    }
}

int main()
{
    std::cout << "Running occupation reconstruction test...\n";

    RNG rng(123);

    // Tiny system: 2 sites, 1 boson
    System sys(2, 1, 1.0, 1.0, 0.0);
    Lattice lat(2, 1);
    Configuration config(sys, lat, 1);

    // Attach Hamiltonian
    BoseHubbardHamiltonian H(sys, lat);
    config.setHamiltonian(&H);

    // Initial Fock state: one boson on site 0
    std::vector<int> fock(sys.size(), 0);
    fock[0] = 1;
    config.initialize(fock);

    Worldline &wl = config.replica(0).worldline();

    // ------------------------------------------------------------
    // Insert mechanical hop 0 -> 1 at tau = 0.3
    // ------------------------------------------------------------
    {
        // kink on site 0 list
        Kink k1{0.3, 0, 0, 1, -1, -1, 0, 0, -1, 0};
        // kink on site 1 list
        Kink k2{0.3, 0, 0, 1, -1, -1, 0, 0, -1, 1};
        wl.insertHop(k1, k2);
        wl.checkConsistency();
    }

    // ------------------------------------------------------------
    // Insert mechanical hop 1 -> 0 at tau = 0.7
    // ------------------------------------------------------------
    {
        // kink on site 1 list
        Kink k1{0.7, 0, 1, 0, -1, -1, 0, 0, -1, 1};
        // kink on site 0 list
        Kink k2{0.7, 0, 1, 0, -1, -1, 0, 0, -1, 0};
        wl.insertHop(k1, k2);
        wl.checkConsistency();
    }

    dump_worldline(wl);

    // ------------------------------------------------------------
    // Query occupations at several tau values
    // ------------------------------------------------------------
    double taus[] = {0.0, 0.2, 0.5, 0.9};

    for (double tau : taus)
    {
        int n0 = wl.occupationAt(0, tau);
        int n1 = wl.occupationAt(1, tau);

        std::cout << "tau = " << tau
                  << "  n0 = " << n0
                  << "  n1 = " << n1 << "\n";
    }

    std::cout << "Occupation reconstruction test completed.\n";
    return 0;
}
