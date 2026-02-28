#include "../src/pimc.hpp"
#include <cassert>
#include <iostream>

using namespace pimc;

// ------------------------------------------------------------
// Helper: build a simple configuration with a flat worldline
// ------------------------------------------------------------
Configuration build_config(int L, int D, int N, RNG &rng)
{
    System sys(L, D, 1.0, 1.0, 0.0);
    Lattice lat(L, D);
    Configuration config(sys, lat, 1);

    // Attach Hamiltonian
    auto *H = new BoseHubbardHamiltonian(sys, lat);
    config.setHamiltonian(H);

    // Random Fock state
    int M = config.latticeSize();
    auto fock = random_fock_state(M, N, rng);
    config.initialize(fock);

    return config;
}

// ------------------------------------------------------------
// Test 1: Hamiltonian attachment
// ------------------------------------------------------------
void test_hamiltonian_attachment()
{
    RNG rng(123);
    System sys(2, 2, 1.0, 1.0, 0.0);
    Lattice lat(2, 2);
    Configuration config(sys, lat, 1);

    BoseHubbardHamiltonian H(sys, lat);
    config.setHamiltonian(&H);

    assert(&config.hamiltonian() == &H);
    std::cout << "[OK] Hamiltonian attachment\n";
}

// ------------------------------------------------------------
// Test 2: Onsite energy formula
// ------------------------------------------------------------
void test_onsite_energy()
{
    System sys(2, 2, 1.0, 1.0, 0.0);
    Lattice lat(2, 2);
    BoseHubbardHamiltonian H(sys, lat);

    assert(H.onsiteEnergy(0) == 0.0);
    assert(H.onsiteEnergy(1) == 0.0);
    assert(H.onsiteEnergy(2) == 1.0);
    assert(H.onsiteEnergy(3) == 3.0);

    std::cout << "[OK] Onsite energy\n";
}

// ------------------------------------------------------------
// Test 3: Hopping amplitude
// ------------------------------------------------------------
void test_hopping_amplitude()
{
    System sys(2, 2, 1.0, 1.0, 0.0);
    Lattice lat(2, 2);
    BoseHubbardHamiltonian H(sys, lat);

    assert(H.hoppingAmplitude(0, 1) == -1.0);
    assert(H.hoppingAmplitude(0, 2) == -1.0);
    assert(H.hoppingAmplitude(0, 3) == 0.0);

    std::cout << "[OK] Hopping amplitude\n";
}

// ------------------------------------------------------------
// Test 4: Diagonal energy
// ------------------------------------------------------------
void test_diagonal_energy()
{
    RNG rng(7);
    Configuration config = build_config(2, 2, 3, rng);
    auto &H = config.hamiltonian();

    const Worldline &wl = config.replica(0).worldline();
    int M = config.latticeSize();

    double E_manual = 0.0;
    for (int site = 0; site < M; ++site)
    {
        int idx = wl.firstKink(site);
        E_manual += H.onsiteEnergy(wl[idx].n);
    }

    assert(std::abs(E_manual - H.diagonalEnergy(config)) < 1e-12);
    std::cout << "[OK] Diagonal energy\n";
}

// ------------------------------------------------------------
// Test 5: Local diagonal energy
// ------------------------------------------------------------
void test_local_diagonal_energy()
{
    RNG rng(11);
    Configuration config = build_config(2, 2, 4, rng);
    auto &H = config.hamiltonian();

    const Worldline &wl = config.replica(0).worldline();
    int M = config.latticeSize();

    for (int site = 0; site < M; ++site)
    {
        int idx = wl.firstKink(site);
        double E_local = H.localDiagonalEnergy(config, 0, site, 0.0);
        double E_expected = H.onsiteEnergy(wl[idx].n);
        assert(std::abs(E_local - E_expected) < 1e-12);
    }

    std::cout << "[OK] Local diagonal energy\n";
}

// ------------------------------------------------------------
// Main test runner
// ------------------------------------------------------------
int main()
{
    test_hamiltonian_attachment();
    test_onsite_energy();
    test_hopping_amplitude();
    test_diagonal_energy();
    test_local_diagonal_energy();

    std::cout << "All Hamiltonian tests PASSED.\n";
    return 0;
}
