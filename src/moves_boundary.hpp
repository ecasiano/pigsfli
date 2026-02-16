#pragma once
#include "pimc.hpp"
#include "moves_mc.hpp"
#include <cmath>

namespace pimc
{

    struct BoundaryPair0
    {
        int idx_fwd;
        int idx_bwd;
        int site_i;
        int site_j;
    };

    // enumerate hop pairs where the earlier kink is at tau=0 and the later in (0,beta)
    inline std::vector<BoundaryPair0> enumerateBoundaryPairs0(const Worldline &wl, double beta)
    {
        auto pairs = enumerateHopPairs(wl);
        std::vector<BoundaryPair0> out;

        double eps = 1e-12;
        for (const auto &hp : pairs)
        {
            double tau_f = wl[hp.idx_fwd].tau;
            double tau_b = wl[hp.idx_bwd].tau;

            double tau_min = std::min(tau_f, tau_b);
            double tau_max = std::max(tau_f, tau_b);

            if (std::fabs(tau_min) < eps && tau_max > eps && tau_max < beta - eps)
            {
                out.push_back({hp.idx_fwd, hp.idx_bwd, hp.site_i, hp.site_j});
            }
        }
        return out;
    }

    // ------------------------------------------------------------
    // Boundary kink pair insertion at tau = 0 (partner in bulk)
    // ------------------------------------------------------------
    class BoundaryKinkPairInsertion0MC : public Move
    {
    public:
        explicit BoundaryKinkPairInsertion0MC(const SimulationParameters &params)
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

            // 1. bonds
            std::vector<std::pair<int, int>> bonds;
            bonds.reserve(M * 2);
            for (int i = 0; i < M; ++i)
                for (int j : lat.neighbors(i))
                    if (i < j)
                        bonds.emplace_back(i, j);

            if (bonds.empty())
                return false;

            int Nb = (int)bonds.size();
            auto [i, j] = bonds[rng.randint(0, Nb - 1)];

            // 2. times: tau1 = 0, tau2 uniform in (0,beta)
            double tau1 = 0.0;
            double tau2 = rng.uniform() * beta_;
            if (tau2 <= 0.0 || tau2 >= beta_)
                return false;

            // 3. occupations at tau=0
            int n_i0 = wl.occupationAt(i, 0.0);
            int n_j0 = wl.occupationAt(j, 0.0);
            if (n_i0 <= 0)
                return false;

            double t = sys.t();
            double abs_t = std::abs(t);

            double me1 = abs_t * std::sqrt((double)n_i0 * (double)(n_j0 + 1));
            double me2 = abs_t * std::sqrt((double)(n_j0 + 1) * (double)n_i0);
            if (me1 == 0.0 || me2 == 0.0)
                return false;

            // 4. diagonal BEFORE
            double S_i_before = diagonalActionSite(wl, H, i, beta_);
            double S_j_before = diagonalActionSite(wl, H, j, beta_);

            // 5. insert kinks: i->j at 0, j->i at tau2
            Kink k1{tau1, 0, i, j, -1, -1, r, r, -1, i};
            Kink k2{tau1, 0, i, j, -1, -1, r, r, -1, j};
            auto [idx1, idx2] = wl.insertHop(k1, k2);

            Kink k3{tau2, 0, j, i, -1, -1, r, r, -1, j};
            Kink k4{tau2, 0, j, i, -1, -1, r, r, -1, i};
            auto [idx3, idx4] = wl.insertHop(k3, k4);

            wl.checkConsistency();

            // 6. diagonal AFTER
            double S_i_after = diagonalActionSite(wl, H, i, beta_);
            double S_j_after = diagonalActionSite(wl, H, j, beta_);

            double dS_diag = (S_i_after + S_j_after) - (S_i_before + S_j_before);
            double logW_ratio = std::log(me1) + std::log(me2) - dS_diag;

            // 7. proposal ratio
            auto boundary_pairs_new = enumerateBoundaryPairs0(wl, beta_);
            int N_pairs_new = (int)boundary_pairs_new.size();
            if (N_pairs_new == 0)
            {
                wl.deleteHop(idx1);
                wl.deleteHop(idx3);
                wl.checkConsistency();
                return false;
            }

            double logP_ins = -std::log((double)Nb) - std::log(beta_); // bond * tau2
            double logP_rem = -std::log((double)N_pairs_new);

            double log_ratio = logW_ratio + (logP_rem - logP_ins);

            if (std::log(rng.uniform()) < log_ratio)
                return true;

            // reject
            wl.deleteHop(idx1);
            wl.deleteHop(idx3);
            wl.checkConsistency();
            return false;
        }

    private:
        double beta_;
    };

    // ------------------------------------------------------------
    // Boundary kink pair removal at tau = 0 (partner in bulk)
    // ------------------------------------------------------------
    class BoundaryKinkPairRemoval0MC : public Move
    {
    public:
        explicit BoundaryKinkPairRemoval0MC(const SimulationParameters &params)
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

            auto boundary_pairs = enumerateBoundaryPairs0(wl, beta_);
            int N_pairs_old = (int)boundary_pairs.size();
            if (N_pairs_old == 0)
                return false;

            int k = rng.randint(0, N_pairs_old - 1);
            BoundaryPair0 bp = boundary_pairs[k];

            int i = bp.site_i;
            int j = bp.site_j;

            // diagonal BEFORE
            double S_i_before = diagonalActionSite(wl, H, i, beta_);
            double S_j_before = diagonalActionSite(wl, H, j, beta_);

            // save kinks for undo
            Kink k_fwd = wl[bp.idx_fwd];
            Kink k_fwd_partner = wl[k_fwd.partner];
            Kink k_bwd = wl[bp.idx_bwd];
            Kink k_bwd_partner = wl[k_bwd.partner];

            // occupations at tau=0 (old config)
            int n_i0 = wl.occupationAt(i, 0.0);
            int n_j0 = wl.occupationAt(j, 0.0);

            double t = sys.t();
            double abs_t = std::abs(t);

            double me1 = abs_t * std::sqrt((double)n_i0 * (double)(n_j0 + 1));
            double me2 = abs_t * std::sqrt((double)(n_j0 + 1) * (double)n_i0);
            if (me1 == 0.0 || me2 == 0.0)
                return false;

            // remove hops
            wl.deleteHop(bp.idx_fwd);
            wl.deleteHop(bp.idx_bwd);
            wl.checkConsistency();

            // diagonal AFTER
            double S_i_after = diagonalActionSite(wl, H, i, beta_);
            double S_j_after = diagonalActionSite(wl, H, j, beta_);

            double dS_diag = (S_i_after + S_j_after) - (S_i_before + S_j_before);
            double logW_ratio = -std::log(me1) - std::log(me2) - dS_diag;

            // proposal ratio
            std::vector<std::pair<int, int>> bonds;
            bonds.reserve(M * 2);
            for (int s = 0; s < M; ++s)
                for (int nb : lat.neighbors(s))
                    if (s < nb)
                        bonds.emplace_back(s, nb);

            int Nb = (int)bonds.size();
            if (Nb == 0)
            {
                // undo
                wl.insertHop(k_fwd, k_fwd_partner);
                wl.insertHop(k_bwd, k_bwd_partner);
                wl.checkConsistency();
                return false;
            }

            double logP_rem = -std::log((double)N_pairs_old);
            double logP_ins = -std::log((double)Nb) - std::log(beta_);

            double log_ratio = logW_ratio + (logP_ins - logP_rem);

            if (std::log(rng.uniform()) < log_ratio)
                return true;

            // reject: undo
            wl.insertHop(k_fwd, k_fwd_partner);
            wl.insertHop(k_bwd, k_bwd_partner);
            wl.checkConsistency();
            return false;
        }

    private:
        double beta_;
    };

} // namespace pimc
