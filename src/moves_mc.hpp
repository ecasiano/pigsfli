#ifndef MOVES_MC_HPP
#define MOVES_MC_HPP

#include "pimc.hpp"
#include <cmath>

namespace pimc
{

    // ============================================================
    // Helper: count hop pairs (forward+backward on a bond)
    // ============================================================

    struct HopPair
    {
        int idx_fwd; // index of forward hop (i -> j)
        int idx_bwd; // index of backward hop (j -> i)
        int site_i;  // one site of the bond
        int site_j;  // the other site
    };

    inline std::vector<HopPair> enumerateHopPairs(const Worldline &wl)
    {
        std::vector<HopPair> pairs;
        int M = wl.latticeSize();

        struct HopInfo
        {
            int idx;
            int src;
            int dest;
            double tau;
        };
        std::vector<HopInfo> hops;

        for (int site = 0; site < M; ++site)
        {
            int idx = wl.firstKink(site);
            while (idx != -1)
            {
                const Kink &k = wl[idx];
                if (k.partner != -1 && idx < k.partner)
                {
                    hops.push_back({idx, k.src, k.dest, k.tau});
                }
                idx = k.next;
            }
        }

        int H = (int)hops.size();
        for (int a = 0; a < H; ++a)
        {
            for (int b = a + 1; b < H; ++b)
            {
                const auto &ha = hops[a];
                const auto &hb = hops[b];

                bool sameBond =
                    (ha.src == hb.dest && ha.dest == hb.src);

                if (!sameBond)
                    continue;

                if (ha.tau < hb.tau)
                {
                    pairs.push_back({ha.idx, hb.idx, ha.src, ha.dest});
                }
                else
                {
                    pairs.push_back({hb.idx, ha.idx, hb.src, hb.dest});
                }
            }
        }

        return pairs;
    }

    // ============================================================
    // Helper: local occupations just before tau
    // ============================================================

    inline int occupation_before_tau(const Worldline &wl, int site, double tau)
    {
        // occupationAt uses kinks up to <= tau; we want just before tau
        // so we query at tau - epsilon; numerically, tau - 1e-12 is enough
        double eps = 1e-12;
        double t_query = tau - eps;
        if (t_query < 0.0)
            t_query += 0.0; // effectively 0; initial kink at tau=0 carries n
        return wl.occupationAt(site, t_query);
    }

    // ============================================================
    // Metropolis kink–antikink insertion (bulk, with bosonic factors)
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

            // --- 1. Build bond list and pick a bond uniformly
            std::vector<std::pair<int, int>> bonds;
            bonds.reserve(M * 2);
            for (int i = 0; i < M; ++i)
            {
                for (int j : lat.neighbors(i))
                {
                    if (i < j)
                        bonds.emplace_back(i, j);
                }
            }
            if (bonds.empty())
                return false;

            int Nb = (int)bonds.size();
            auto [i, j] = bonds[rng.randint(0, Nb - 1)];

            // --- 2. Propose two times uniformly in [0, beta), order them
            double tau1 = rng.uniform() * beta_;
            double tau2 = rng.uniform() * beta_;
            if (tau2 < tau1)
                std::swap(tau1, tau2);

            // --- 3. Occupations just before tau1 and tau2
            int n_i_tau1 = occupation_before_tau(wl, i, tau1);
            int n_j_tau1 = occupation_before_tau(wl, j, tau1);

            if (n_i_tau1 <= 0)
                return false;

            double t = sys.t();
            double abs_t = std::abs(t);

            // Use occupations before tau1 for both hops (|2,0> -> |1,1> -> |2,0>)
            double me1 = abs_t * std::sqrt((double)n_i_tau1 * (double)(n_j_tau1 + 1));
            double me2 = abs_t * std::sqrt((double)(n_j_tau1 + 1) * (double)n_i_tau1);

            if (me1 == 0.0 || me2 == 0.0)
                return false;

            // --- 4. Diagonal action BEFORE (only sites i and j)
            double S_i_before = diagonalActionSite(wl, H, i, beta_);
            double S_j_before = diagonalActionSite(wl, H, j, beta_);

            // --- 5. Insert one hop pair: i->j at tau1, j->i at tau2
            Kink k1{tau1, 0, i, j, -1, -1, r, r, -1, i};
            Kink k2{tau1, 0, i, j, -1, -1, r, r, -1, j};
            auto [idx1, idx2] = wl.insertHop(k1, k2);

            Kink k3{tau2, 0, j, i, -1, -1, r, r, -1, j};
            Kink k4{tau2, 0, j, i, -1, -1, r, r, -1, i};
            auto [idx3, idx4] = wl.insertHop(k3, k4);

            wl.checkConsistency();

            // --- 6. Diagonal action AFTER
            double S_i_after = diagonalActionSite(wl, H, i, beta_);
            double S_j_after = diagonalActionSite(wl, H, j, beta_);

            double dS_diag = (S_i_after + S_j_after) - (S_i_before + S_j_before);

            // --- 7. Weight ratio: |me1 * me2| * exp(-dS_diag)
            double logW_ratio = std::log(me1) + std::log(me2) - dS_diag;

            // --- 8. Proposal ratio: P_rem / P_ins
            // P_ins = (1/Nb) * (2 / beta^2)
            // P_rem = 1 / N_pairs_new
            auto pairs_new = enumerateHopPairs(wl);
            int N_pairs_new = (int)pairs_new.size();
            if (N_pairs_new == 0)
            {
                wl.deleteHop(idx1);
                wl.deleteHop(idx3);
                wl.checkConsistency();
                return false;
            }

            double logP_ins = -std::log((double)Nb) + std::log(2.0) - 2.0 * std::log(beta_);
            double logP_rem = -std::log((double)N_pairs_new);

            double log_ratio = logW_ratio + (logP_rem - logP_ins);

            if (std::log(rng.uniform()) < log_ratio)
                return true;

            // reject: undo
            wl.deleteHop(idx1);
            wl.deleteHop(idx3);
            wl.checkConsistency();
            return false;
        }

    private:
        double beta_;
    };

    // ============================================================
    // Metropolis kink–antikink removal (bulk, with bosonic factors)
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
            const Lattice &lat = C.lattice();
            const System &sys = C.system();
            const Hamiltonian &H = C.hamiltonian();

            int M = C.latticeSize();
            if (M == 0)
                return false;

            // --- 1. Enumerate hop pairs
            auto pairs = enumerateHopPairs(wl);
            int N_pairs_old = (int)pairs.size();
            if (N_pairs_old == 0)
                return false;

            // --- 2. Choose one pair uniformly
            int k = rng.randint(0, N_pairs_old - 1);
            HopPair hp = pairs[k];

            int i = hp.site_i;
            int j = hp.site_j;

            // --- 3. Diagonal BEFORE
            double S_i_before = diagonalActionSite(wl, H, i, beta_);
            double S_j_before = diagonalActionSite(wl, H, j, beta_);

            // Save kinks for undo
            Kink k_fwd = wl[hp.idx_fwd];
            Kink k_fwd_partner = wl[k_fwd.partner];
            Kink k_bwd = wl[hp.idx_bwd];
            Kink k_bwd_partner = wl[k_bwd.partner];

            double tau1 = k_fwd.tau;
            double tau2 = k_bwd.tau;

            // --- 4. Occupations just before tau1 and tau2 (in the old config)
            int n_i_tau1 = occupation_before_tau(wl, i, tau1);
            int n_j_tau1 = occupation_before_tau(wl, j, tau1);

            int n_j_tau2 = occupation_before_tau(wl, j, tau2);
            int n_i_tau2 = occupation_before_tau(wl, i, tau2);

            double t = sys.t();
            double abs_t = std::abs(t);

            double me1 = abs_t * std::sqrt((double)n_i_tau1 * (double)(n_j_tau1 + 1));
            double me2 = abs_t * std::sqrt((double)n_j_tau2 * (double)(n_i_tau2 + 1));

            if (me1 == 0.0 || me2 == 0.0)
                return false;

            // --- 5. Remove the two hops
            wl.deleteHop(hp.idx_fwd);
            wl.deleteHop(hp.idx_bwd);
            wl.checkConsistency();

            // --- 6. Diagonal AFTER
            double S_i_after = diagonalActionSite(wl, H, i, beta_);
            double S_j_after = diagonalActionSite(wl, H, j, beta_);

            double dS_diag = (S_i_after + S_j_after) - (S_i_before + S_j_before);

            // --- 7. Weight ratio: W_new / W_old = 1 / (|me1 * me2|) * exp(-dS_diag)
            double logW_ratio = -std::log(me1) - std::log(me2) - dS_diag;

            // --- 8. Proposal ratio: P_ins / P_rem
            // P_rem = 1 / N_pairs_old
            // P_ins = (1/Nb) * (2 / beta^2)
            std::vector<std::pair<int, int>> bonds;
            bonds.reserve(M * 2);
            for (int s = 0; s < M; ++s)
            {
                for (int nb : lat.neighbors(s))
                {
                    if (s < nb)
                        bonds.emplace_back(s, nb);
                }
            }
            int Nb = (int)bonds.size();
            if (Nb == 0)
            {
                // undo and return false
                wl.insertHop(k_fwd, k_fwd_partner);
                wl.insertHop(k_bwd, k_bwd_partner);
                wl.checkConsistency();
                return false;
            }

            double logP_rem = -std::log((double)N_pairs_old);
            double logP_ins = -std::log((double)Nb) + std::log(2.0) - 2.0 * std::log(beta_);

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

    // ============================================================
    // NOTE: boundary moves can be reintroduced later, but for now
    // we keep only bulk moves for clean ED benchmarking.
    // ============================================================

} // namespace pimc

#endif // MOVES_MC_HPP
