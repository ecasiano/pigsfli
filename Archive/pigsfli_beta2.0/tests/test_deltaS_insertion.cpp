#include "../src/pimc.hpp"
#include <cassert>
#include <iostream>
#include <cmath>

using namespace pimc;

int main()
{
    std::cout << "Running test_deltaS_insertion...\n";

    double beta = 1.0;

    // Simple 2-site Bose-Hubbard system
    // U = 1, mu = 0, t = 0.5 (arbitrary)
    System sys(2, 1, 0.5, 1.0, 0.0);
    Lattice lat(2, 1);
    Configuration config(sys, lat, 1);

    BoseHubbardHamiltonian H(sys, lat);
    config.setHamiltonian(&H);

    // Initial Fock state: 1 boson on site 0
    std::vector<int> fock = {1, 0};
    config.initialize(fock);

    Worldline &wl = config.replica(0).worldline();

    // Choose a bond 0 <-> 1
    int i = 0;
    int j = 1;

    // Choose insertion times
    double tau1 = 0.3;
    double tau2 = 0.7;

    // Compute brute-force diagonal action BEFORE
    double S_before = diagonalAction(config, 0, beta);

    // Compute site-local BEFORE (for MC move)
    double S_i_before = diagonalActionSite(wl, H, i, beta);
    double S_j_before = diagonalActionSite(wl, H, j, beta);

    // Perform mechanical insertion (same as MC move)
    Kink k1{tau1, 0, i, j, -1, -1, 0, 0, -1, i};
    Kink k2{tau1, 0, i, j, -1, -1, 0, 0, -1, j};
    auto [idx1, idx2] = wl.insertHop(k1, k2);

    Kink k3{tau2, 0, i, j, -1, -1, 0, 0, -1, i};
    Kink k4{tau2, 0, i, j, -1, -1, 0, 0, -1, j};
    auto [idx3, idx4] = wl.insertHop(k3, k4);

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

    std::cout << "test_deltaS_insertion passed.\n";
    return 0;
}
