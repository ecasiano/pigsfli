#ifndef PIMC_HPP
#define PIMC_HPP

#include <iostream>
#include <vector>
#include <cmath>
#include <chrono>
#include <iomanip>
#include <fstream>
#include <cstdlib>
#include <sstream>
#include <cstring>
#include <random>
#include <stdexcept>

namespace pimc
{

    // ============================================================
    // SECTION 1 — Utilities (RNG, helpers)
    // ============================================================

    class RNG
    {
    public:
        RNG() : engine(std::random_device{}()) {}
        RNG(uint64_t seed) : engine(seed) {}

        double uniform()
        {
            return std::uniform_real_distribution<double>(0.0, 1.0)(engine);
        }

        int randint(int a, int b)
        {
            return std::uniform_int_distribution<int>(a, b)(engine);
        }

    private:
        std::mt19937_64 engine;
    };

    inline std::vector<int> random_fock_state(int M, int N, RNG &rng)
    {
        std::vector<int> fock(M, 0);
        for (int p = 0; p < N; p++)
        {
            int site = rng.randint(0, M - 1);
            fock[site] += 1;
        }
        return fock;
    }

    // ============================================================
    // SECTION 2 — Static Model (System, Lattice)
    // ============================================================

    class System
    {
    public:
        System(int L, int D, double t, double U, double mu)
            : L_(L), D_(D), M_(std::pow(L, D)), t_(t), U_(U), mu_(mu) {}

        int size() const { return M_; }
        int dimension() const { return D_; }
        int linearSize() const { return L_; }

        double hopping() const { return t_; }
        double interaction() const { return U_; }
        double chemicalPotential() const { return mu_; }

        double t() const { return t_; }
        double U() const { return U_; }
        double mu() const { return mu_; }

    private:
        int L_, D_, M_;
        double t_, U_, mu_;
    };

    class Lattice
    {
    public:
        Lattice(int L, int D)
            : L(L), D(D), M(std::pow(L, D))
        {
            build_adjacency();
        }

        int size() const { return M; }
        int dimension() const { return D; }
        int linearSize() const { return L; }

        const std::vector<int> &neighbors(int site) const
        {
            return adjacency[site];
        }

    private:
        int L, D, M;
        std::vector<std::vector<int>> adjacency;

        void build_adjacency()
        {
            adjacency.assign(M, {});
            for (int site = 0; site < M; site++)
            {
                int x = site % L;
                int y = (site / L) % L;
                int z = site / (L * L);

                if (D >= 1)
                {
                    adjacency[site].push_back((x + 1) % L + L * y + L * L * z);
                    adjacency[site].push_back((x - 1 + L) % L + L * y + L * L * z);
                }
                if (D >= 2)
                {
                    adjacency[site].push_back(x + L * ((y + 1) % L) + L * L * z);
                    adjacency[site].push_back(x + L * ((y - 1 + L) % L) + L * L * z);
                }
                if (D == 3)
                {
                    adjacency[site].push_back(x + L * y + L * L * ((z + 1) % L));
                    adjacency[site].push_back(x + L * y + L * L * ((z - 1 + L) % L));
                }
            }
        }
    };

    // ============================================================
    // SECTION 3 — Worldline Data Structures (Kink, Worldline)
    // ============================================================

    struct Kink
    {
        double tau;
        int n;
        int src, dest; // physical hop direction
        int prev, next;
        int src_replica;
        int dest_replica;
        int partner;
        int site; // which site list this kink belongs to
    };

    class Worldline
    {
    public:
        Worldline(int M)
            : M(M), head(M, -1) {}

        int latticeSize() const { return M; }

        int kinkCount() const
        {
            return kinks.size();
        }

        const Kink &operator[](int idx) const { return kinks.at(idx); }
        Kink &operator[](int idx) { return kinks.at(idx); }

        int firstKink(int site) const { return head.at(site); }

        void setHead(int site, int idx) { head.at(site) = idx; }

        int addKink(const Kink &k)
        {
            kinks.push_back(k);
            return static_cast<int>(kinks.size()) - 1;
        }

        // ------------------------------------------------------------
        // Query occupation n_i(tau)
        // ------------------------------------------------------------
        int occupationAt(int site, double tau) const
        {
            int idx = head.at(site);
            if (idx == -1)
                return 0;

            int n = kinks[idx].n;

            idx = kinks[idx].next;
            while (idx != -1 && kinks[idx].tau <= tau)
            {
                const Kink &k = kinks[idx];

                if (k.site == site && k.src != k.dest)
                {
                    if (site == k.src)
                        n -= 1;
                    else if (site == k.dest)
                        n += 1;
                }

                idx = k.next;
            }

            return n;
        }

    private:
        int insertKinkOrdered(int site, const Kink &k)
        {
            int idx = addKink(k);
            Kink &nk = kinks[idx];
            nk.site = site; // bookkeeping

            int h = head[site];

            if (h == -1)
            {
                nk.prev = -1;
                nk.next = -1;
                head[site] = idx;
                return idx;
            }

            if (nk.tau < kinks[h].tau)
            {
                nk.prev = -1;
                nk.next = h;
                kinks[h].prev = idx;
                head[site] = idx;
                return idx;
            }

            int cur = h;
            while (kinks[cur].next != -1 &&
                   kinks[kinks[cur].next].tau < nk.tau)
            {
                cur = kinks[cur].next;
            }

            int nxt = kinks[cur].next;
            nk.prev = cur;
            nk.next = nxt;
            kinks[cur].next = idx;

            if (nxt != -1)
                kinks[nxt].prev = idx;

            return idx;
        }

    public:
        std::pair<int, int> insertHop(const Kink &k1, const Kink &k2)
        {
            int idx1 = insertKinkOrdered(k1.site, k1);
            int idx2 = insertKinkOrdered(k2.site, k2);

            kinks[idx1].partner = idx2;
            kinks[idx2].partner = idx1;

            return {idx1, idx2};
        }

        void removeKink(int idx)
        {
            Kink &k = kinks.at(idx);
            int site = k.site;

            if (k.prev != -1)
                kinks[k.prev].next = k.next;
            else
                head[site] = k.next;

            if (k.next != -1)
                kinks[k.next].prev = k.prev;

            k.prev = k.next = k.partner = -1;
        }

        void deleteHop(int idx)
        {
            int p = kinks.at(idx).partner;
            removeKink(p);
            removeKink(idx);
        }

        void checkConsistency() const
        {
            for (int site = 0; site < M; ++site)
            {
                int current = head[site];
                int steps = 0;

                while (current != -1)
                {
                    if (current < 0 || current >= (int)kinks.size())
                        throw std::runtime_error("Invalid kink index");

                    const Kink &k = kinks[current];

                    if (k.partner != -1)
                    {
                        if (k.partner < 0 || k.partner >= (int)kinks.size())
                            throw std::runtime_error("Invalid partner index");
                        if (kinks[k.partner].partner != current)
                            throw std::runtime_error("Partner mismatch");
                    }

                    if (k.next == current)
                        throw std::runtime_error("Self-loop detected");

                    current = k.next;

                    if (++steps > (int)kinks.size())
                        throw std::runtime_error("Cycle detected");
                }
            }
        }

    private:
        int M;
        std::vector<int> head;
        std::vector<Kink> kinks;
    };

    // ============================================================
    // Forward declaration
    // ============================================================

    inline void initialize_worldline_from_fock(
        Worldline &wl,
        const std::vector<int> &fock_state);

    // ============================================================
    // SECTION 4 — Replica and Configuration
    // ============================================================

    class Replica
    {
    public:
        explicit Replica(int M)
            : wl(M) {}

        Worldline &worldline() { return wl; }
        const Worldline &worldline() const { return wl; }

    private:
        Worldline wl;
    };

    class Hamiltonian;

    class Configuration
    {
    public:
        Configuration(const System &sys, const Lattice &lat, int num_replicas)
            : sys(sys), lat(lat), H(nullptr)
        {
            int M = sys.size();
            replicas.reserve(num_replicas);
            for (int r = 0; r < num_replicas; r++)
                replicas.emplace_back(M);
        }

        void setHamiltonian(Hamiltonian *h) { H = h; }

        Hamiltonian &hamiltonian() { return *H; }
        const Hamiltonian &hamiltonian() const { return *H; }

        int replicasCount() const { return replicas.size(); }
        int latticeSize() const { return sys.size(); }

        Replica &replica(int r) { return replicas.at(r); }
        const Replica &replica(int r) const { return replicas.at(r); }

        const System &system() const { return sys; }
        const Lattice &lattice() const { return lat; }

        void initialize(const std::vector<int> &fock_state)
        {
            for (auto &rep : replicas)
                initialize_worldline_from_fock(rep.worldline(), fock_state);
        }

    private:
        System sys;
        Lattice lat;
        Hamiltonian *H;
        std::vector<Replica> replicas;
    };

    // ============================================================
    // SECTION 5 — Hamiltonian Abstraction Layer
    // ============================================================

    class Hamiltonian
    {
    public:
        virtual ~Hamiltonian() {}

        virtual double onsiteEnergy(int n) const = 0;
        virtual double hoppingAmplitude(int src, int dest) const = 0;

        virtual double localDiagonalEnergy(const Configuration &C,
                                           int replica,
                                           int site,
                                           double tau) const = 0;

        virtual double diagonalEnergy(const Configuration &C,
                                      int replica = 0) const = 0;
    };

    class BoseHubbardHamiltonian : public Hamiltonian
    {
    public:
        BoseHubbardHamiltonian(const System &sys, const Lattice &lat)
            : sys_(sys), lat_(lat),
              t_(sys.hopping()), U_(sys.interaction()), mu_(sys.chemicalPotential()) {}

        double onsiteEnergy(int n) const override
        {
            return 0.5 * U_ * n * (n - 1) - mu_ * n;
        }

        double hoppingAmplitude(int src, int dest) const override
        {
            const auto &neigh = lat_.neighbors(src);
            for (int s : neigh)
                if (s == dest)
                    return -t_;
            return 0.0;
        }

        double localDiagonalEnergy(const Configuration &C,
                                   int replica,
                                   int site,
                                   double /*tau*/) const override
        {
            const Worldline &wl = C.replica(replica).worldline();
            int idx = wl.firstKink(site);
            if (idx < 0)
                return 0.0;
            return onsiteEnergy(wl[idx].n);
        }

        double diagonalEnergy(const Configuration &C,
                              int replica = 0) const override
        {
            const Worldline &wl = C.replica(replica).worldline();
            int M = C.latticeSize();
            double E = 0.0;

            for (int site = 0; site < M; ++site)
            {
                int idx = wl.firstKink(site);
                if (idx < 0)
                    continue;
                E += onsiteEnergy(wl[idx].n);
            }
            return E;
        }

    private:
        const System &sys_;
        const Lattice &lat_;
        double t_, U_, mu_;
    };

    // ============================================================
    // Diagonal action for a single site (projector, open in tau)
    // ============================================================

    inline double diagonalActionSite(
        const Worldline &wl,
        const Hamiltonian &H,
        int site,
        double beta)
    {
        int idx = wl.firstKink(site);
        if (idx == -1)
            return 0.0;

        double S = 0.0;

        const Kink *cur = &wl[idx];
        int n = cur->n;
        double prevTau = 0.0;

        idx = cur->next;
        while (idx != -1)
        {
            const Kink &k = wl[idx];

            double dt = k.tau - prevTau;
            if (dt < 0.0)
                throw std::runtime_error("Non-monotonic tau in worldline");

            S += dt * H.onsiteEnergy(n);

            if (k.site == site && k.src != k.dest)
            {
                if (site == k.src)
                    n -= 1;
                else if (site == k.dest)
                    n += 1;
            }

            prevTau = k.tau;
            idx = k.next;
        }

        double dt = beta - prevTau;
        if (dt < 0.0)
            throw std::runtime_error("beta smaller than last kink tau");

        S += dt * H.onsiteEnergy(n);

        return S;
    }

    inline double diagonalAction(
        const Configuration &C,
        int replica,
        double beta)
    {
        const Hamiltonian &H = C.hamiltonian();
        const Worldline &wl = C.replica(replica).worldline();
        int M = C.latticeSize();

        double S = 0.0;
        for (int site = 0; site < M; ++site)
            S += diagonalActionSite(wl, H, site, beta);

        return S;
    }

    // ============================================================
    // SECTION 6 — Simulation parameters
    // ============================================================

    class SimulationParameters
    {
    public:
        SimulationParameters(double beta)
            : beta_(beta), eta_(0.0) {}

        double beta() const { return beta_; }
        double eta() const { return eta_; }

        void setEta(double eta) { eta_ = eta; }

    private:
        double beta_;
        double eta_;
    };

    // ============================================================
    // SECTION 7 — Monte Carlo Framework (MECHANICAL ONLY)
    // ============================================================

    class Move
    {
    public:
        virtual ~Move() {}
        virtual bool attempt(Configuration &C, RNG &rng) = 0;
    };

    // ------------------------------
    // Pure mechanical insertion
    // ------------------------------
    class KinkAntikinkInsertion : public Move
    {
    public:
        KinkAntikinkInsertion(const SimulationParameters &params)
            : beta_(params.beta()) {}

        bool attempt(Configuration &C, RNG &rng) override
        {
            if (C.replicasCount() == 0)
                return false;

            int r = 0;
            Worldline &wl = C.replica(r).worldline();
            const Lattice &lat = C.lattice();
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

            // hop i -> j at tau1
            Kink k1{tau1, 0, i, j, -1, -1, r, r, -1, i}; // stored on site i
            Kink k2{tau1, 0, i, j, -1, -1, r, r, -1, j}; // stored on site j
            wl.insertHop(k1, k2);

            // hop i -> j at tau2
            Kink k3{tau2, 0, i, j, -1, -1, r, r, -1, i};
            Kink k4{tau2, 0, i, j, -1, -1, r, r, -1, j};
            wl.insertHop(k3, k4);

            return true;
        }

    private:
        double beta_;
    };

    // ------------------------------
    // Pure mechanical removal
    // ------------------------------
    class KinkAntikinkRemoval : public Move
    {
    public:
        KinkAntikinkRemoval(const SimulationParameters &params)
            : beta_(params.beta()) {}

        bool attempt(Configuration &C, RNG &rng) override
        {
            if (C.replicasCount() == 0)
                return false;

            int r = 0;
            Worldline &wl = C.replica(r).worldline();
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

            wl.deleteHop(idx);

            return true;
        }

    private:
        double beta_;
    };

    // ============================================================
    // SECTION 8 — Estimators and Simulation
    // ============================================================

    class Estimator
    {
    public:
        virtual ~Estimator() {}
        virtual double measure(const Configuration &C) = 0;
    };

    class Simulation
    {
    public:
        Simulation(System sys, Configuration conf, RNG rng, SimulationParameters params)
            : system(sys), config(conf), rng(rng), params(params) {}

        void addMove(Move *move) { moves.push_back(move); }
        void addEstimator(Estimator *est) { estimators.push_back(est); }

        void run(int steps)
        {
            for (int s = 0; s < steps; s++)
            {
                for (auto &move : moves)
                    move->attempt(config, rng);
                for (auto &est : estimators)
                    est->measure(config);
            }
        }

    private:
        System system;
        Configuration config;
        RNG rng;
        SimulationParameters params;

        std::vector<Move *> moves;
        std::vector<Estimator *> estimators;
    };

    // ============================================================
    // SECTION 9 — initialize_worldline_from_fock
    // ============================================================

    inline void initialize_worldline_from_fock(
        Worldline &wl,
        const std::vector<int> &fock_state)
    {
        int M = wl.latticeSize();

        for (int site = 0; site < M; site++)
        {
            Kink k;
            k.tau = 0.0;
            k.n = fock_state[site];
            k.src = site;
            k.dest = site;
            k.prev = -1;
            k.next = -1;
            k.src_replica = 0;
            k.dest_replica = 0;
            k.partner = -1;
            k.site = site;

            int idx = wl.addKink(k);
            wl.setHead(site, idx);
        }
    }

    // ============================================================
    // SECTION 10 — Miscellaneous
    // ============================================================

    inline void test()
    {
        std::cout << "PIMC skeleton with continuous-time kinks loaded.\n";
    }

} // namespace pimc

#endif // PIMC_HPP
