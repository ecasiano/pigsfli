#pragma once
#include "pimc.hpp"

namespace pimc
{

    // ============================================================
    // Diagonal Energy Estimator
    // ============================================================
    struct DiagonalEnergyEstimator : public Estimator
    {
        double beta;

        explicit DiagonalEnergyEstimator(double beta_) : beta(beta_) {}

        double measure(const Configuration &C) override
        {
            // Always use replica 0 for now
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
    // E_kin = - N_hops / beta
    // ============================================================
    struct KineticEnergyEstimator : public Estimator
    {
        double beta;

        explicit KineticEnergyEstimator(double beta_) : beta(beta_) {}

        double measure(const Configuration &C) override
        {
            const Worldline &wl = C.replica(0).worldline();

            int count = 0;

            // Count all kinks (each kink is one hop)
            for (int site = 0; site < C.latticeSize(); ++site)
            {
                int idx = wl.firstKink(site);
                while (idx != -1)
                {
                    count++;
                    idx = wl[idx].next;
                }
            }

            return -double(count) / beta;
        }
    };

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

} // namespace pimc
