#pragma once
#include "pimc.hpp"

namespace pimc
{

    // ============================================================
    // Helper: count physical hops (each hop is one src–dest pair)
    // ============================================================
    inline int countPhysicalHops(const Configuration &C)
    {
        const Worldline &wl = C.replica(0).worldline();
        int M = C.latticeSize();
        int count = 0;

        for (int site = 0; site < M; ++site)
        {
            int idx = wl.firstKink(site);
            while (idx != -1)
            {
                const Kink &k = wl[idx];
                if (k.partner != -1 && idx < k.partner)
                    count++; // one physical hop
                idx = k.next;
            }
        }
        return count;
    }

    // ============================================================
    // Diagonal Energy Estimator
    // E_diag = S_diag / beta
    // ============================================================
    struct DiagonalEnergyEstimator : public Estimator
    {
        double beta;

        explicit DiagonalEnergyEstimator(double beta_) : beta(beta_) {}

        double measure(const Configuration &C) override
        {
            double Sdiag = diagonalAction(C, 0, beta);
            return Sdiag / beta;
        }
    };

    // ============================================================
    // Density Estimator for a single site
    // ============================================================
    struct DensityEstimator : public Estimator
    {
        int site;
        double beta;

        DensityEstimator(int site_, double beta_)
            : site(site_), beta(beta_) {}

        double measure(const Configuration &C) override
        {
            const Worldline &wl = C.replica(0).worldline();

            int idx = wl.firstKink(site);
            if (idx == -1)
                return 0.0;

            double integral = 0.0;
            double prevTau = 0.0;
            int n = wl[idx].n;

            idx = wl[idx].next;
            while (idx != -1)
            {
                const Kink &k = wl[idx];
                double dt = k.tau - prevTau;
                integral += dt * n;

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

            integral += (beta - prevTau) * n;

            return integral / beta;
        }
    };

    // ============================================================
    // Total Occupation Estimator
    // ============================================================
    struct TotalOccupationEstimator : public Estimator
    {
        double beta;

        explicit TotalOccupationEstimator(double beta_) : beta(beta_) {}

        double measure(const Configuration &C) override
        {
            int M = C.latticeSize();
            double total = 0.0;

            for (int i = 0; i < M; i++)
            {
                DensityEstimator dens(i, beta);
                total += dens.measure(C);
            }

            return total;
        }
    };

    // ============================================================
    // Kinetic Energy Estimator
    // E_kin = - t * N_hops / beta
    // ============================================================
    struct KineticEnergyEstimator : public Estimator
    {
        double beta;

        explicit KineticEnergyEstimator(double beta_) : beta(beta_) {}

        double measure(const Configuration &C) override
        {
            int N_hops = countPhysicalHops(C);
            double t = C.system().t();
            return -t * static_cast<double>(N_hops) / beta;
        }
    };

    // ============================================================
    // Total Energy Estimator
    // E_tot = E_diag + E_kin
    // ============================================================
    struct TotalEnergyEstimator : public Estimator
    {
        double beta;

        explicit TotalEnergyEstimator(double beta_) : beta(beta_) {}

        double measure(const Configuration &C) override
        {
            DiagonalEnergyEstimator Ed(beta);
            KineticEnergyEstimator Ek(beta);
            return Ed.measure(C) + Ek.measure(C);
        }
    };

    // ============================================================
    // Mid-time Energy Estimator (consistent kinetic definition)
    // ============================================================
    class MidTimeEnergyEstimator
    {
    public:
        explicit MidTimeEnergyEstimator(double beta)
            : beta_(beta) {}

        double measure(const Configuration &C) const
        {
            const System &sys = C.system();
            const Lattice &lat = C.lattice();
            const Worldline &wl = C.replica(0).worldline();

            double tau_mid = 0.5 * beta_;
            int M = lat.size();

            // Diagonal energy at mid-time
            double E_diag = 0.0;
            for (int i = 0; i < M; ++i)
            {
                int n = wl.occupationAt(i, tau_mid);

                double U = sys.U();
                double mu = sys.mu();

                double e_i = 0.0;
                e_i += 0.5 * U * n * (n - 1);
                e_i -= mu * n;

                E_diag += e_i;
            }

            // Kinetic energy from physical hop count
            int N_hops = countPhysicalHops(C);
            double E_kin = -sys.t() * static_cast<double>(N_hops) / beta_;

            return E_diag + E_kin;
        }

    private:
        double beta_;
    };

} // namespace pimc
