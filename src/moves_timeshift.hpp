#pragma once
#include "pimc.hpp"
#include "moves_mc.hpp"
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

            auto pairs = enumerateHopPairs(wl);
            int N_pairs = (int)pairs.size();
            if (N_pairs == 0)
                return false;

            int k = rng.randint(0, N_pairs - 1);
            HopPair hp = pairs[k];

            int i = hp.site_i;
            int j = hp.site_j;

            // Save old kinks for undo
            Kink k_fwd_old = wl[hp.idx_fwd];
            Kink k_fwd_partner_old = wl[k_fwd_old.partner];
            Kink k_bwd_old = wl[hp.idx_bwd];
            Kink k_bwd_partner_old = wl[k_bwd_old.partner];

            double tau1_old = k_fwd_old.tau;
            double tau2_old = k_bwd_old.tau;

            // Propose new times
            double tau1_new = rng.uniform() * beta_;
            double tau2_new = rng.uniform() * beta_;
            if (tau2_new < tau1_new)
                std::swap(tau1_new, tau2_new);

            double S_i_before = diagonalActionSite(wl, H, i, beta_);
            double S_j_before = diagonalActionSite(wl, H, j, beta_);

            // Remove old hops
            wl.deleteHop(hp.idx_fwd);
            wl.deleteHop(hp.idx_bwd);

            // Build new kinks with updated times
            Kink k_fwd_new = k_fwd_old;
            Kink k_fwd_partner_new = k_fwd_partner_old;
            Kink k_bwd_new = k_bwd_old;
            Kink k_bwd_partner_new = k_bwd_partner_old;

            k_fwd_new.tau = tau1_new;
            k_fwd_partner_new.tau = tau1_new;
            k_bwd_new.tau = tau2_new;
            k_bwd_partner_new.tau = tau2_new;

            auto [idx1_new, idx2_new] = wl.insertHop(k_fwd_new, k_fwd_partner_new);
            auto [idx3_new, idx4_new] = wl.insertHop(k_bwd_new, k_bwd_partner_new);

            wl.checkConsistency();

            double S_i_after = diagonalActionSite(wl, H, i, beta_);
            double S_j_after = diagonalActionSite(wl, H, j, beta_);

            double dS_diag = (S_i_after + S_j_after) - (S_i_before + S_j_before);
            double log_ratio = -dS_diag; // symmetric proposal

            if (std::log(rng.uniform()) < log_ratio)
                return true;

            // Reject: undo
            wl.deleteHop(idx1_new);
            wl.deleteHop(idx3_new);

            auto [idx1_old, idx2_old] = wl.insertHop(k_fwd_old, k_fwd_partner_old);
            auto [idx3_old, idx4_old] = wl.insertHop(k_bwd_old, k_bwd_partner_old);

            wl.checkConsistency();
            return false;
        }

    private:
        double beta_;
    };

} // namespace pimc
