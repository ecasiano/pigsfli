#include "../src/pimc.hpp"
#include <cassert>
#include <iostream>
#include <map>

using namespace pimc;

int main()
{
    RNG rng(999);

    System sys(3, 1, 1.0, 1.0, 0.0);
    Lattice lat(3, 1);
    Configuration config(sys, lat, 1);

    std::vector<int> occ = {1, 0, 0}; // brute force
    config.initialize(occ);

    Worldline &wl = config.replica(0).worldline();

    struct Hop
    {
        double tau;
        int i, j;
    };
    std::vector<Hop> hops;

    // generate 20 random hops
    for (int k = 0; k < 20; k++)
    {
        int i = rng.randint(0, 2);
        const auto &neigh = lat.neighbors(i);
        int j = neigh[rng.randint(0, neigh.size() - 1)];
        double tau = rng.uniform();

        hops.push_back({tau, i, j});

        // brute force update
        for (int s = 0; s < 3; s++)
        {
            if (s == i)
                occ[s]--;
            if (s == j)
                occ[s]++;
        }

        // insert into worldline
        wl.insertHop(
            Kink{tau, 0, i, j, -1, -1, 0, 0, -1, i},
            Kink{tau, 0, i, j, -1, -1, 0, 0, -1, j});
    }

    wl.checkConsistency();

    // sample random taus
    for (double tau : {0.1, 0.3, 0.5, 0.7, 0.9})
    {
        for (int s = 0; s < 3; s++)
        {
            int n = wl.occupationAt(s, tau);

            // brute force recompute at tau
            int brute[3] = {1, 0, 0};
            for (auto &h : hops)
            {
                if (h.tau <= tau)
                {
                    brute[h.i]--;
                    brute[h.j]++;
                }
            }

            assert(n == brute[s]);
        }
    }

    std::cout << "test_random_sequence passed.\n";
}
