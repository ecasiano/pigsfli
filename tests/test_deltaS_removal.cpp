#include "../src/pimc.hpp"
#include <cassert>
#include <iostream>
#include <cmath>

using namespace pimc;

int main()
{
    std::cout << "Running test_deltaS_removal...\n";

    double beta = 1.0;

    // Simple 2-site Bose-Hubbard system
    System sys(2, 1, 0.5, 1.0, 0.0); // t=0.5, U=1, mu=0
    Lattice lat(2, 1);
    Configuration config(sys, lat, 1);

    BoseHubbardHamiltonian H(sys, lat);
    config.setHamiltonian(&H);

    // Initial Fock state: 1 boson on site 0
    std::vector<int> fock = {1, 0};
    config.initialize(fock);

    Worldline &wl = config.replica(0).worldline();

    // Build a known worldline: 0->1 at 0.3, 1->0 at 0.7
    int i = 0, j = 1;

    auto [idx1, idx2] = wl.insertHop(
        Kink{0.3, 0, i, j, -1, -1, 0, 0, -1, i},
        Kink{0.3, 0, i, j, -1, -1, 0, 0, -1, j});

    auto [idx3, idx4] = wl.insertHop(
        Kink{0.7, 0, i, j, -1, -1, 0, 0, -1, i},
        Kink{0.7, 0, i, j, -1, -1, 0, 0, -1, j});

    wl.checkConsistency();

    // Compute brute-force diagonal action BEFORE
    double S_before = diagonalAction(config, 0, beta);

    // Compute site-local BEFORE (for MC move)
    double S_i_before = diagonalActionSite(wl, H, i, beta);
    double S_j_before = diagonalActionSite(wl, H, j, beta);

    // Remove one hop pair mechanically
    // Choose the first pair (idx1, idx2)
    wl.deleteHop(idx1);

    wl.checkConsistency();

    // Compute brute-force diagonal action AFTER
    double S_after = diagonalAction(config, 0, beta);

    // Compute site-local AFTER (for MC move)
    double S_i_after = diagonalActionSite(wl, H, i, beta);
    double S_j_after = diagonalActionSite(wl, H, j, beta);

    // Brute-force ΔS
    double dS_bruteforce = S_after - S_before;

    // MC-move ΔS (diagonal part only)
    double dS_move = (S_i_after + S_j_after) - (S_i_before + S_j_before);

    // They must match
    double diff = std::abs(dS_bruteforce - dS_move);
    std::cout << "ΔS_bruteforce = " << dS_bruteforce << "\n";
    std::cout << "ΔS_move       = " << dS_move << "\n";
    std::cout << "difference    = " << diff << "\n";

    assert(diff < 1e-12);

    std::cout << "test_deltaS_removal passed.\n";
    return 0;
}
