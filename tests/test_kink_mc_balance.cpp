#include "../src/pimc.hpp"
#include "../src/moves_mc.hpp"
#include <iostream>
#include <cassert>

using namespace pimc;

int main()
{
    std::cout << "Running micro detailed-balance test...\n";

    RNG rng(123);

    // Tiny system: 2 sites, 1D, 1 boson
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

    // MC moves
    SimulationParameters params(1.0); // β = 1
    KinkAntikinkInsertionMC ins(params);
    KinkAntikinkRemovalMC rem(params);

    Worldline &wl = config.replica(0).worldline();

    int accepted_ins = 0;
    int accepted_rem = 0;
    int attempts_ins = 0;
    int attempts_rem = 0;

    // Alternate insertion and removal
    for (int n = 0; n < 5000; ++n)
    {
        // Try insertion
        attempts_ins++;
        bool acc_ins = ins.attempt(config, rng);
        if (acc_ins)
            accepted_ins++;
        wl.checkConsistency();

        // Try removal
        attempts_rem++;
        bool acc_rem = rem.attempt(config, rng);
        if (acc_rem)
            accepted_rem++;
        wl.checkConsistency();
    }

    std::cout << "Insertion accepted: " << accepted_ins
              << " / " << attempts_ins << "\n";
    std::cout << "Removal accepted:   " << accepted_rem
              << " / " << attempts_rem << "\n";

    // Sanity checks
    assert(accepted_ins > 0);
    assert(accepted_rem > 0);

    std::cout << "Micro detailed-balance test PASSED.\n";
    return 0;
}
