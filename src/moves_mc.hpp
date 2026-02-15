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
            const Hamiltonian &H = C.hamiltonian();

            int M = lat.size();
            if (M == 0)
                return false;

            int i = rng.randint(0, M - 1);
            const auto &neigh = lat.neighbors(i);
            if (neigh.empty())
                return false;
            int j = neigh[rng.randint(0, (int)neigh.size() - 1)];

            double tau1 = rng.uniform() * beta_;
            double tau2 = rng.uniform() * beta_;
            if (tau2 < tau1)
                std::swap(tau1, tau2);

            // Diagonal action BEFORE on affected sites
            double S_i_before = diagonalActionSite(wl, H, i, beta_);
            double S_j_before = diagonalActionSite(wl, H, j, beta_);

            // Insert two hop pairs mechanically
            Kink k1{tau1, 0, i, j, -1, -1, r, r, -1, i};
            Kink k2{tau1, 0, i, j, -1, -1, r, r, -1, j};
            auto [idx1, idx2] = wl.insertHop(k1, k2);

            Kink k3{tau2, 0, i, j, -1, -1, r, r, -1, i};
            Kink k4{tau2, 0, i, j, -1, -1, r, r, -1, j};
            auto [idx3, idx4] = wl.insertHop(k3, k4);

            // Diagonal action AFTER
            double S_i_after = diagonalActionSite(wl, H, i, beta_);
            double S_j_after = diagonalActionSite(wl, H, j, beta_);

            double dS_diag = (S_i_after + S_j_after) - (S_i_before + S_j_before);

            // Kinetic part: 2 pairs = 4 hops
            double t = sys.hopping();
            double dS_kin = -4.0 * std::log(std::abs(t) + 1e-12);

            double dS = dS_diag + dS_kin;
            double log_accept = -dS;

            if (std::log(rng.uniform()) < log_accept)
                return true;

            // Reject: undo
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
            const Hamiltonian &H = C.hamiltonian();

            int M = C.latticeSize();

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

            int site_i = wl[idx].site;
            int site_j = wl[partner].site;

            // Save for undo
            Kink kA = wl[idx];
            Kink kB = wl[partner];

            // Diagonal action BEFORE
            double S_i_before = diagonalActionSite(wl, H, site_i, beta_);
            double S_j_before = diagonalActionSite(wl, H, site_j, beta_);

            // Remove one hop pair
            wl.deleteHop(idx);

            // Diagonal action AFTER
            double S_i_after = diagonalActionSite(wl, H, site_i, beta_);
            double S_j_after = diagonalActionSite(wl, H, site_j, beta_);

            double dS_diag = (S_i_after + S_j_after) - (S_i_before + S_j_before);

            // Kinetic: removing 2 hops
            double t = sys.hopping();
            double dS_kin = +2.0 * std::log(std::abs(t) + 1e-12);

            double dS = dS_diag + dS_kin;
            double log_accept = -dS;

            if (std::log(rng.uniform()) < log_accept)
                return true;

            // Reject: undo
            wl.insertHop(kA, kB);
            return false;
        }

    private:
        double beta_;
    };

} // namespace pimc

#endif // MOVES_MC_HPP
