#include "../src/pimc.hpp"
#include <cassert>
#include <iostream>

using namespace pimc;

int main()
{
    std::cout << "Running test_diagonal_action...\n";

    double beta = 1.0;
    System sys(2, 1, 0.0, 1.0, 0.0); // t=0 so only diagonal matters
    Lattice lat(2, 1);
    Configuration config(sys, lat, 1);

    BoseHubbardHamiltonian H(sys, lat);
    config.setHamiltonian(&H);

    // initial: one boson on site 0
    std::vector<int> fock = {1, 0};
    config.initialize(fock);

    Worldline &wl = config.replica(0).worldline();

    // 0 -> 1 at tau=0.3
    wl.insertHop(
        Kink{0.3, 0, 0, 1, -1, -1, 0, 0, -1, 0},
        Kink{0.3, 0, 0, 1, -1, -1, 0, 0, -1, 1});

    // 1 -> 0 at tau=0.7
    wl.insertHop(
        Kink{0.7, 0, 1, 0, -1, -1, 0, 0, -1, 1},
        Kink{0.7, 0, 1, 0, -1, -1, 0, 0, -1, 0});

    wl.checkConsistency();

    // Analytic S_diag:
    // U=1, mu=0
    // n=1 => E = 0.5*1*0 - 0 = 0
    // so E_onsite(1) = 0 everywhere, S_diag = 0
    double S = diagonalAction(config, 0, beta);
    assert(std::abs(S - 0.0) < 1e-12);

    std::cout << "test_diagonal_action passed.\n";
    return 0;
}
