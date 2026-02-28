#include "../src/pimc.hpp"
#include <iostream>
#include <cassert>

using namespace pimc;

// Helper: check that every site has exactly one kink at tau=0
bool check_flat_worldline(const Worldline &wl, int M, const std::vector<int> &fock)
{
    for (int site = 0; site < M; site++)
    {
        int idx = wl.firstKink(site);
        if (idx < 0)
            return false;

        const Kink &k = wl[idx];

        if (k.tau != 0.0)
            return false;
        if (k.n != fock[site])
            return false;
        if (k.prev != -1 || k.next != -1)
            return false;
        if (k.partner != -1)
            return false;

        // Ensure no extra kinks
        if (k.next != -1)
            return false;
    }
    return true;
}

void test_initialization()
{
    std::cout << "Running initialization test...\n";

    int L = 2, D = 2;
    int M = std::pow(L, D);
    int N = 5;

    RNG rng(1234);
    System sys(L, D, 1.0, 1.0, 0.0);
    Lattice lat(L, D);
    Configuration config(sys, lat, 1);

    auto fock = random_fock_state(M, N, rng);

    // Check Fock state sum
    int sum = 0;
    for (int x : fock)
        sum += x;
    assert(sum == N);

    config.initialize(fock);

    const Worldline &wl = config.replica(0).worldline();

    // Check structural correctness
    assert(check_flat_worldline(wl, M, fock));

    // Check consistency
    wl.checkConsistency();

    std::cout << "Initialization test PASSED.\n\n";
}

void test_extreme_fock_states()
{
    std::cout << "Running extreme Fock state tests...\n";

    int L = 2, D = 2;
    int M = std::pow(L, D);

    RNG rng(999);
    System sys(L, D, 1.0, 1.0, 0.0);
    Lattice lat(L, D);
    Configuration config(sys, lat, 1);

    // Case 1: All bosons on one site
    {
        std::vector<int> fock(M, 0);
        fock[0] = 20;

        config.initialize(fock);
        const Worldline &wl = config.replica(0).worldline();

        assert(check_flat_worldline(wl, M, fock));
        wl.checkConsistency();
    }

    // Case 2: One boson per site
    {
        std::vector<int> fock(M, 1);

        config.initialize(fock);
        const Worldline &wl = config.replica(0).worldline();

        assert(check_flat_worldline(wl, M, fock));
        wl.checkConsistency();
    }

    // Case 3: Zero bosons
    {
        std::vector<int> fock(M, 0);

        config.initialize(fock);
        const Worldline &wl = config.replica(0).worldline();

        assert(check_flat_worldline(wl, M, fock));
        wl.checkConsistency();
    }

    std::cout << "Extreme Fock state tests PASSED.\n\n";
}

void test_stress()
{
    std::cout << "Running stress test (10,000 initializations)...\n";

    int L = 2, D = 2;
    int M = std::pow(L, D);
    int N = 10;

    RNG rng(42);
    System sys(L, D, 1.0, 1.0, 0.0);
    Lattice lat(L, D);
    Configuration config(sys, lat, 1);

    for (int i = 0; i < 10000; i++)
    {
        auto fock = random_fock_state(M, N, rng);
        config.initialize(fock);

        const Worldline &wl = config.replica(0).worldline();
        wl.checkConsistency();
    }

    std::cout << "Stress test PASSED.\n\n";
}

int main()
{
    test_initialization();
    test_extreme_fock_states();
    test_stress();

    std::cout << "All tests PASSED.\n";
    return 0;
}
