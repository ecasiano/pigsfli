#pragma once
#include "pimc.hpp"
#include <cmath>

namespace pimc
{

    class BoundaryInsertionMC : public Move
    {
    public:
        explicit BoundaryInsertionMC(const SimulationParameters &params)
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

            // small epsilon to avoid exact 0 or beta
            double eps = 1e-9 * beta_;

            // 1. Pick a bond (i,j)
            std::vector<std::pair<int, int>> bonds;
            for (int i = 0; i < M; ++i)
                for (int j : lat.neighbors(i))
                    if (i < j)
                        bonds.emplace_back(i, j);

            if (bonds.empty())
                return false;

            auto [i, j] = bonds[rng.randint(0, (int)bonds.size() - 1)];

            // 2. Choose boundary: near 0 or near beta
            bool at_beta = (rng.uniform() < 0.5);
            double tau = at_beta ? (beta_ - eps) : eps;

            // 3. Occupations just before tau
            int n_i = wl.occupationAt(i, tau);
            int n_j = wl.occupationAt(j, tau);

            if (n_i <= 0)
                return false;

            double t = sys.t();
            double abs_t = std::abs(t);

            double me = abs_t * std::sqrt((double)n_i * (double)(n_j + 1));
            if (me == 0.0)
                return false;

            // 4. Diagonal BEFORE
            double S_i_before = diagonalActionSite(wl, H, i, beta_);
            double S_j_before = diagonalActionSite(wl, H, j, beta_);

            // 5. Insert hop pair at boundary tau
            Kink k1{tau, 0, i, j, -1, -1, r, r, -1, i};
            Kink k2{tau, 0, i, j, -1, -1, r, r, -1, j};
            auto [idx1, idx2] = wl.insertHop(k1, k2);

            wl.checkConsistency();

            // 6. Diagonal AFTER
            double S_i_after = diagonalActionSite(wl, H, i, beta_);
            double S_j_after = diagonalActionSite(wl, H, j, beta_);

            double dS = (S_i_after + S_j_after) - (S_i_before + S_j_before);

            double log_ratio = std::log(me) - dS;

            if (std::log(rng.uniform()) < log_ratio)
                return true;

            wl.deleteHop(idx1);
            wl.checkConsistency();
            return false;
        }

    private:
        double beta_;
    };

    class BoundaryRemovalMC : public Move
    {
    public:
        explicit BoundaryRemovalMC(const SimulationParameters &params)
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

            double eps = 1e-9 * beta_;

            // 1. Find hop pairs whose tau is near 0 or near beta
            std::vector<int> boundary_pairs;
            for (int site = 0; site < M; ++site)
            {
                int idx = wl.firstKink(site);
                while (idx != -1)
                {
                    const Kink &k = wl[idx];
                    if (k.partner != -1 && idx < k.partner)
                    {
                        double tau = k.tau;
                        if (tau < eps || tau > beta_ - eps)
                            boundary_pairs.push_back(idx);
                    }
                    idx = k.next;
                }
            }

            if (boundary_pairs.empty())
                return false;

            int idx = boundary_pairs[rng.randint(0, (int)boundary_pairs.size() - 1)];
            Kink k = wl[idx];
            int partner_idx = k.partner;
            Kink k_partner = wl[partner_idx];

            int i = k.src;
            int j = k.dest;
            double tau = k.tau;

            // 3. Diagonal BEFORE
            double S_i_before = diagonalActionSite(wl, H, i, beta_);
            double S_j_before = diagonalActionSite(wl, H, j, beta_);

            int n_i = wl.occupationAt(i, tau);
            int n_j = wl.occupationAt(j, tau);

            double t = sys.t();
            double abs_t = std::abs(t);

            double me = abs_t * std::sqrt((double)n_i * (double)(n_j + 1));
            if (me == 0.0)
                return false;

            wl.deleteHop(idx);
            wl.checkConsistency();

            double S_i_after = diagonalActionSite(wl, H, i, beta_);
            double S_j_after = diagonalActionSite(wl, H, j, beta_);

            double dS = (S_i_after + S_j_after) - (S_i_before + S_j_before);

            double log_ratio = -std::log(me) - dS;

            if (std::log(rng.uniform()) < log_ratio)
                return true;

            wl.insertHop(k, k_partner);
            wl.checkConsistency();
            return false;
        }

    private:
        double beta_;
    };

} // namespace pimc
