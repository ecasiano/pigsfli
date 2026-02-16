#ifndef MOVES_MC_HPP
#define MOVES_MC_HPP

#include "pimc.hpp"
#include <cmath>

namespace pimc
{

    // ============================================================
    // Helper: count hop pairs
    // ============================================================

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
    // Metropolis kink–antikink insertion (bulk)
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

            // pick random bond i-j
            int i = rng.randint(0, M - 1);
            const auto &neigh = lat.neighbors(i);
            if (neigh.empty())
                return false;
            int j = neigh[rng.randint(0, (int)neigh.size() - 1)];

            // pick two times
            double tau1 = rng.uniform() * beta_;
            double tau2 = rng.uniform() * beta_;
            if (tau2 < tau1)
                std::swap(tau1, tau2);

            // diagonal action BEFORE
            double S_i_before = diagonalActionSite(wl, H, i, beta_);
            double S_j_before = diagonalActionSite(wl, H, j, beta_);

            // mechanical insertion: two hops at tau1, two at tau2
            Kink k1{tau1, 0, i, j, -1, -1, r, r, -1, i};
            Kink k2{tau1, 0, i, j, -1, -1, r, r, -1, j};
            auto [idx1, idx2] = wl.insertHop(k1, k2);

            Kink k3{tau2, 0, i, j, -1, -1, r, r, -1, i};
            Kink k4{tau2, 0, i, j, -1, -1, r, r, -1, j};
            auto [idx3, idx4] = wl.insertHop(k3, k4);

            wl.checkConsistency();

            // diagonal action AFTER
            double S_i_after = diagonalActionSite(wl, H, i, beta_);
            double S_j_after = diagonalActionSite(wl, H, j, beta_);

            double dS_diag = (S_i_after + S_j_after) - (S_i_before + S_j_before);

            // kinetic: 4 hops added
            double t = sys.hopping();
            double dS_kin = -4.0 * std::log(std::abs(t) + 1e-12);

            double dS = dS_diag + dS_kin;
            double log_accept = -dS;

            if (std::log(rng.uniform()) < log_accept)
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
    // Metropolis kink–antikink removal (bulk)
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

            // collect candidate hop pairs
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

            // diagonal BEFORE
            double S_i_before = diagonalActionSite(wl, H, site_i, beta_);
            double S_j_before = diagonalActionSite(wl, H, site_j, beta_);

            // save for undo
            Kink kA = wl[idx];
            Kink kB = wl[partner];

            // remove pair
            wl.deleteHop(idx);
            wl.checkConsistency();

            // diagonal AFTER
            double S_i_after = diagonalActionSite(wl, H, site_i, beta_);
            double S_j_after = diagonalActionSite(wl, H, site_j, beta_);

            double dS_diag = (S_i_after + S_j_after) - (S_i_before + S_j_before);

            // kinetic: 2 hops removed
            double t = sys.hopping();
            double dS_kin = +2.0 * std::log(std::abs(t) + 1e-12);

            double dS = dS_diag + dS_kin;
            double log_accept = -dS;

            if (std::log(rng.uniform()) < log_accept)
                return true;

            // reject: undo
            wl.insertHop(kA, kB);
            wl.checkConsistency();
            return false;
        }

    private:
        double beta_;
    };

    // ============================================================
    // Boundary kink–antikink insertion at tau = 0
    // ============================================================

    class BoundaryKinkAntikinkInsertion0MC : public Move
    {
    public:
        BoundaryKinkAntikinkInsertion0MC(const SimulationParameters &params)
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

            // pick random bond i-j
            int i = rng.randint(0, M - 1);
            const auto &neigh = lat.neighbors(i);
            if (neigh.empty())
                return false;
            int j = neigh[rng.randint(0, (int)neigh.size() - 1)];

            // times: tau1 = 0, tau2 in (0, beta)
            double tau1 = 0.0;
            double tau2 = rng.uniform() * beta_;
            if (tau2 <= 0.0 || tau2 >= beta_)
                return false;

            // diagonal BEFORE
            double S_i_before = diagonalActionSite(wl, H, i, beta_);
            double S_j_before = diagonalActionSite(wl, H, j, beta_);

            // mechanical insertion
            Kink k1{tau1, 0, i, j, -1, -1, r, r, -1, i};
            Kink k2{tau1, 0, i, j, -1, -1, r, r, -1, j};
            auto [idx1, idx2] = wl.insertHop(k1, k2);

            Kink k3{tau2, 0, i, j, -1, -1, r, r, -1, i};
            Kink k4{tau2, 0, i, j, -1, -1, r, r, -1, j};
            auto [idx3, idx4] = wl.insertHop(k3, k4);

            wl.checkConsistency();

            // diagonal AFTER
            double S_i_after = diagonalActionSite(wl, H, i, beta_);
            double S_j_after = diagonalActionSite(wl, H, j, beta_);

            double dS_diag = (S_i_after + S_j_after) - (S_i_before + S_j_before);

            // kinetic: 4 hops added
            double t = sys.hopping();
            double dS_kin = -4.0 * std::log(std::abs(t) + 1e-12);

            double dS = dS_diag + dS_kin;
            double log_accept = -dS;

            if (std::log(rng.uniform()) < log_accept)
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

    // ============================================================
    // Boundary kink–antikink removal at tau = 0
    // ============================================================

    class BoundaryKinkAntikinkRemoval0MC : public Move
    {
    public:
        BoundaryKinkAntikinkRemoval0MC(const SimulationParameters &params)
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

            // collect hop pairs at tau = 0
            std::vector<int> candidates;
            for (int site = 0; site < M; ++site)
            {
                int idx = wl.firstKink(site);
                while (idx != -1)
                {
                    const Kink &k = wl[idx];
                    if (k.partner != -1 &&
                        idx < k.partner &&
                        std::abs(k.tau - 0.0) < 1e-12)
                    {
                        candidates.push_back(idx);
                    }
                    idx = k.next;
                }
            }

            if (candidates.empty())
                return false;

            int idx = candidates[rng.randint(0, (int)candidates.size() - 1)];
            int partner = wl[idx].partner;

            int site_i = wl[idx].site;
            int site_j = wl[partner].site;

            // diagonal BEFORE
            double S_i_before = diagonalActionSite(wl, H, site_i, beta_);
            double S_j_before = diagonalActionSite(wl, H, site_j, beta_);

            // save for undo
            Kink kA = wl[idx];
            Kink kB = wl[partner];

            // remove
            wl.deleteHop(idx);
            wl.checkConsistency();

            // diagonal AFTER
            double S_i_after = diagonalActionSite(wl, H, site_i, beta_);
            double S_j_after = diagonalActionSite(wl, H, site_j, beta_);

            double dS_diag = (S_i_after + S_j_after) - (S_i_before + S_j_before);

            // kinetic: 2 hops removed
            double t = sys.hopping();
            double dS_kin = +2.0 * std::log(std::abs(t) + 1e-12);

            double dS = dS_diag + dS_kin;
            double log_accept = -dS;

            if (std::log(rng.uniform()) < log_accept)
                return true;

            // reject
            wl.insertHop(kA, kB);
            wl.checkConsistency();
            return false;
        }

    private:
        double beta_;
    };

    // ============================================================
    // Boundary kink–antikink insertion at tau = beta
    // ============================================================

    class BoundaryKinkAntikinkInsertionBetaMC : public Move
    {
    public:
        BoundaryKinkAntikinkInsertionBetaMC(const SimulationParameters &params)
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

            // pick random bond i-j
            int i = rng.randint(0, M - 1);
            const auto &neigh = lat.neighbors(i);
            if (neigh.empty())
                return false;
            int j = neigh[rng.randint(0, (int)neigh.size() - 1)];

            // times: tau1 in (0, beta), tau2 = beta
            double tau1 = rng.uniform() * beta_;
            double tau2 = beta_;

            if (tau1 <= 0.0 || tau1 >= beta_)
                return false;

            // diagonal BEFORE
            double S_i_before = diagonalActionSite(wl, H, i, beta_);
            double S_j_before = diagonalActionSite(wl, H, j, beta_);

            // mechanical insertion
            Kink k1{tau1, 0, i, j, -1, -1, r, r, -1, i};
            Kink k2{tau1, 0, i, j, -1, -1, r, r, -1, j};
            auto [idx1, idx2] = wl.insertHop(k1, k2);

            Kink k3{tau2, 0, i, j, -1, -1, r, r, -1, i};
            Kink k4{tau2, 0, i, j, -1, -1, r, r, -1, j};
            auto [idx3, idx4] = wl.insertHop(k3, k4);

            wl.checkConsistency();

            // diagonal AFTER
            double S_i_after = diagonalActionSite(wl, H, i, beta_);
            double S_j_after = diagonalActionSite(wl, H, j, beta_);

            double dS_diag = (S_i_after + S_j_after) - (S_i_before + S_j_before);

            // kinetic: 4 hops added
            double t = sys.hopping();
            double dS_kin = -4.0 * std::log(std::abs(t) + 1e-12);

            double dS = dS_diag + dS_kin;
            double log_accept = -dS;

            if (std::log(rng.uniform()) < log_accept)
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

    // ============================================================
    // Boundary kink–antikink removal at tau = beta
    // ============================================================

    class BoundaryKinkAntikinkRemovalBetaMC : public Move
    {
    public:
        BoundaryKinkAntikinkRemovalBetaMC(const SimulationParameters &params)
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

            // collect hop pairs at tau = beta
            std::vector<int> candidates;
            for (int site = 0; site < M; ++site)
            {
                int idx = wl.firstKink(site);
                while (idx != -1)
                {
                    const Kink &k = wl[idx];
                    if (k.partner != -1 &&
                        idx < k.partner &&
                        std::abs(k.tau - beta_) < 1e-12)
                    {
                        candidates.push_back(idx);
                    }
                    idx = k.next;
                }
            }

            if (candidates.empty())
                return false;

            int idx = candidates[rng.randint(0, (int)candidates.size() - 1)];
            int partner = wl[idx].partner;

            int site_i = wl[idx].site;
            int site_j = wl[partner].site;

            // diagonal BEFORE
            double S_i_before = diagonalActionSite(wl, H, site_i, beta_);
            double S_j_before = diagonalActionSite(wl, H, site_j, beta_);

            // save for undo
            Kink kA = wl[idx];
            Kink kB = wl[partner];

            // remove
            wl.deleteHop(idx);
            wl.checkConsistency();

            // diagonal AFTER
            double S_i_after = diagonalActionSite(wl, H, site_i, beta_);
            double S_j_after = diagonalActionSite(wl, H, site_j, beta_);

            double dS_diag = (S_i_after + S_j_after) - (S_i_before + S_j_before);

            // kinetic: 2 hops removed
            double t = sys.hopping();
            double dS_kin = +2.0 * std::log(std::abs(t) + 1e-12);

            double dS = dS_diag + dS_kin;
            double log_accept = -dS;

            if (std::log(rng.uniform()) < log_accept)
                return true;

            // reject
            wl.insertHop(kA, kB);
            wl.checkConsistency();
            return false;
        }

    private:
        double beta_;
    };

} // namespace pimc

#endif // MOVES_MC_HPP
