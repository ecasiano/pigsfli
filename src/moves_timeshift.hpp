#pragma once
#include "pimc.hpp"
#include "moves_mc.hpp" // for HopPair, enumerateHopPairs
#include <cmath>

namespace pimc
{

    class HopPairTimeShiftMC : public Move
    {
    public:
        explicit HopPairTimeShiftMC(const SimulationParameters &params)
            : beta_(params.beta()) {}

        bool attempt(Configuration &C, RNG &rng) override
        {
            if (C.replicasCount() == 0)
                return false;

            int r = 0;
            Worldline &wl = C.replica(r).worldline();
            const Hamiltonian &H = C.hamiltonian();

            // 1. Enumerate hop pairs
            auto pairs = enumerateHopPairs(wl);
            int N_pairs = (int)pairs.size();
            if (N_pairs == 0)
                return false;

            // 2. Choose one pair uniformly
            int k = rng.randint(0, N_pairs - 1);
            HopPair hp = pairs[k];

            int i = hp.site_i;
            int j = hp.site_j;

            // Original times
            double tau1_old = wl[hp.idx_fwd].tau;
            double tau2_old = wl[hp.idx_bwd].tau;

            // 3. Propose new times uniformly in [0, beta), ordered
            double tau1_new = rng.uniform() * beta_;
            double tau2_new = rng.uniform() * beta_;
            if (tau2_new < tau1_new)
                std::swap(tau1_new, tau2_new);

            // 4. Diagonal action BEFORE (sites i and j)
            double S_i_before = diagonalActionSite(wl, H, i, beta_);
            double S_j_before = diagonalActionSite(wl, H, j, beta_);

            // 5. Temporarily move the kinks
            Kink k_fwd_old = wl[hp.idx_fwd];
            Kink k_bwd_old = wl[hp.idx_bwd];

            wl[hp.idx_fwd].tau = tau1_new;
            wl[hp.idx_bwd].tau = tau2_new;
            wl.checkConsistency();

            // 6. Diagonal action AFTER
            double S_i_after = diagonalActionSite(wl, H, i, beta_);
            double S_j_after = diagonalActionSite(wl, H, j, beta_);

            double dS_diag = (S_i_after + S_j_after) - (S_i_before + S_j_before);

            // 7. Proposal is symmetric (uniform in [0,beta)^2 with ordering),
            // so acceptance is just exp(-dS_diag)
            double log_ratio = -dS_diag;

            if (std::log(rng.uniform()) < log_ratio)
                return true;

            // 8. Reject: restore old times
            wl[hp.idx_fwd].tau = tau1_old;
            wl[hp.idx_bwd].tau = tau2_old;
            wl.checkConsistency();
            return false;
        }

    private:
        double beta_;
    };

} // namespace pimc
