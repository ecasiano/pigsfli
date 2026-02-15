#include "../src/pimc.hpp"
#include <cassert>
#include <iostream>

using namespace pimc;

void test_flat_worldline()
{
    RNG rng(123);
    int M = 4;
    std::vector<int> fock = {1, 0, 2, 1};

    Worldline wl(M);
    pimc::initialize_worldline_from_fock(wl, fock);

    for (int site = 0; site < M; ++site)
    {
        int idx = wl.firstKink(site);
        assert(idx >= 0);
        assert(wl[idx].n == fock[site]);
        assert(wl[idx].prev == -1);
        assert(wl[idx].next == -1);
    }

    std::cout << "[OK] Flat worldline initialization\n";
}

void test_insert_and_delete_hop()
{
    RNG rng(42);
    Worldline wl(2);

    pimc::initialize_worldline_from_fock(wl, {1, 1});

    Kink k1{0.3, 0, 0, 1, -1, -1, 0, 0, -1};
    Kink k2{0.3, 0, 1, 0, -1, -1, 0, 0, -1};

    auto [i1, i2] = wl.insertHop(k1, k2);

    assert(wl[i1].partner == i2);
    assert(wl[i2].partner == i1);

    wl.deleteHop(i1);

    std::cout << "[OK] Hop insertion/deletion\n";
}

void test_time_ordering()
{
    RNG rng(1);
    Worldline wl(1);

    pimc::initialize_worldline_from_fock(wl, {1});

    // Insert hops at different times
    Kink k1{0.7, 0, 0, 0, -1, -1, 0, 0, -1};
    Kink k2{0.7, 0, 0, 0, -1, -1, 0, 0, -1};

    Kink k3{0.3, 0, 0, 0, -1, -1, 0, 0, -1};
    Kink k4{0.3, 0, 0, 0, -1, -1, 0, 0, -1};

    wl.insertHop(k1, k2);
    wl.insertHop(k3, k4);

    int cur = wl.firstKink(0);
    double last_tau = -1.0;
    int steps = 0;

    while (cur != -1)
    {
        assert(wl[cur].tau >= last_tau);
        last_tau = wl[cur].tau;
        cur = wl[cur].next;
        steps++;
        assert(steps < 20);
    }

    std::cout << "[OK] Time-ordered insertion\n";
}

void test_random_insertions_consistency()
{
    RNG rng(99);
    Worldline wl(3);

    pimc::initialize_worldline_from_fock(wl, {1, 1, 1});

    for (int n = 0; n < 50; ++n)
    {
        int site = rng.randint(0, 2);
        double tau = rng.uniform();

        Kink k1{tau, 0, site, site, -1, -1, 0, 0, -1};
        Kink k2{tau, 0, site, site, -1, -1, 0, 0, -1};

        wl.insertHop(k1, k2);
    }

    wl.checkConsistency();

    std::cout << "[OK] Random insertion consistency\n";
}

int main()
{
    test_flat_worldline();
    test_insert_and_delete_hop();
    test_time_ordering();
    test_random_insertions_consistency();

    std::cout << "All worldline tests PASSED.\n";
    return 0;
}
