#ifndef MOVES_MC_HPP
#define MOVES_MC_HPP

#include "pimc.hpp"
#include <cmath>

namespace pimc
{

    // Helper: count hop pairs
    inline int countHopPairs(const Worldline &wl)
    {
        int count = 0;
        int M = wl.latticeSize();

        for (int site = 0; site < M; ++site)
        {
            int idx = wl.firstKink(site);
            while (idx != -1)
            {
                const Kink &k = wl[idx];
                if (k.partner != -1 && idx < k.partner)
                    count++;
                idx = k.next;
            }
        }
        return count;
    }

    // ============================================================
    // Metropolis kink–antikink insertion
    // ============================================================

    class KinkAntikinkInsertionMC : public Move
    {
    public:
        KinkAntikinkInsertionMC(const SimulationParameters &params)
            : beta_(params.beta()) {}

        bool attempt(Configuration &C, RNG &rng) override
        {
            if (C.replicasCount() == 0)
                return false;

            int r = 0;
            Worldline &wl = C.replica(r).worldline();
            const Lattice &lat = C.lattice();
            const System &sys = C.system();

            int M = lat.size();
            if (M == 0)
                return false;

            // Choose bond
            int i = rng.randint(0, M - 1);
            const auto &neigh = lat.neighbors(i);
            if (neigh.empty())
                return false;
            int j = neigh[rng.randint(0, (int)neigh.size() - 1)];

            // Choose times
            double tau1 = rng.uniform() * beta_;
            double tau2 = rng.uniform() * beta_;
            if (tau2 < tau1)
                std::swap(tau1, tau2);

            // Count BEFORE
            int Np_old = countHopPairs(wl);

            // Insert two hop pairs
            Kink k1{tau1, 0, i, j, -1, -1, r, r, -1};
            Kink k2{tau1, 0, i, j, -1, -1, r, r, -1};
            auto [idx1, idx2] = wl.insertHop(k1, k2);

            Kink k3{tau2, 0, i, j, -1, -1, r, r, -1};
            Kink k4{tau2, 0, i, j, -1, -1, r, r, -1};
            auto [idx3, idx4] = wl.insertHop(k3, k4);

            // Count AFTER
            int Np_new = countHopPairs(wl);

            // Fake diagonal action (use magnitude of change)
            double alpha = 0.2;
            double dS_diag = beta_ * alpha * std::abs(Np_new - Np_old);

            // Kinetic action (2 pairs inserted)
            double t = sys.hopping();
            double dS_kin = -2.0 * std::log(std::abs(t) + 1e-12);

            double dS = dS_diag + dS_kin;

            double log_accept = -dS;
            if (std::log(rng.uniform()) < log_accept)
                return true;

            // Reject
            wl.deleteHop(idx1);
            wl.deleteHop(idx3);
            return false;
        }

    private:
        double beta_;
    };

    // ============================================================
    // Metropolis kink–antikink removal
    // ============================================================

    class KinkAntikinkRemovalMC : public Move
    {
    public:
        KinkAntikinkRemovalMC(const SimulationParameters &params)
            : beta_(params.beta()) {}

        bool attempt(Configuration &C, RNG &rng) override
        {
            if (C.replicasCount() == 0)
                return false;

            int r = 0;
            Worldline &wl = C.replica(r).worldline();
            const System &sys = C.system();

            int M = C.latticeSize();

            // Collect hop pairs
            std::vector<int> candidates;
            for (int site = 0; site < M; ++site)
            {
                int idx = wl.firstKink(site);
                while (idx != -1)
                {
                    const Kink &k = wl[idx];
                    if (k.partner != -1 && idx < k.partner)
                        candidates.push_back(idx);
                    idx = k.next;
                }
            }

            if (candidates.empty())
                return false;

            int idx = candidates[rng.randint(0, (int)candidates.size() - 1)];
            int partner = wl[idx].partner;

            // Save for undo
            Kink kA = wl[idx];
            Kink kB = wl[partner];

            // Count BEFORE
            int Np_old = countHopPairs(wl);

            // Remove one hop pair
            wl.deleteHop(idx);

            // Count AFTER
            int Np_new = countHopPairs(wl);

            // Fake diagonal action (use magnitude of change)
            double alpha = 0.2;
            double dS_diag = beta_ * alpha * std::abs(Np_new - Np_old);

            // Kinetic action (1 pair removed)
            double t = sys.hopping();
            double dS_kin = +1.0 * std::log(std::abs(t) + 1e-12);

            double dS = dS_diag + dS_kin;

            double log_accept = -dS;
            if (std::log(rng.uniform()) < log_accept)
                return true;

            // Reject
            wl.insertHop(kA, kB);
            return false;
        }

    private:
        double beta_;
    };

} // namespace pimc

#endif // MOVES_MC_HPP
