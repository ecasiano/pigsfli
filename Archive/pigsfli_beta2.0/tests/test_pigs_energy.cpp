#include "../src/pimc.hpp"
#include "../src/moves_mc.hpp"
#include "../src/moves_timeshift.hpp"
#include "../src/moves_boundary.hpp"
#include "../src/estimators.hpp"

#include <iostream>
#include <vector>
#include <cmath>

using namespace pimc;

// ------------------------------------------------------------
// Move statistics tracker
// ------------------------------------------------------------
struct MoveStats
{
    long proposed = 0;
    long accepted = 0;
    long rejected = 0;

    void record(bool acc)
    {
        proposed++;
        if (acc)
            accepted++;
        else
            rejected++;
    }

    void print(const std::string &name) const
    {
        double accRate = proposed > 0 ? double(accepted) / proposed : 0.0;
        std::cout << name << ": proposed=" << proposed
                  << " accepted=" << accepted
                  << " rejected=" << rejected
                  << " accRate=" << accRate << "\n";
    }
};

// ------------------------------------------------------------
// Mid-time diagonal energy
// ------------------------------------------------------------
double diagonal_energy_mid(const Configuration &C, double beta)
{
    const System &sys = C.system();
    const Worldline &wl = C.replica(0).worldline();
    double tau_mid = 0.5 * beta;

    double E = 0.0;
    for (int site = 0; site < C.latticeSize(); ++site)
    {
        int n = wl.occupationAt(site, tau_mid);
        E += 0.5 * sys.U() * n * (n - 1) - sys.mu() * n;
    }
    return E;
}

// ------------------------------------------------------------
// Global kinetic estimator
// ------------------------------------------------------------
double kinetic_energy_global(const Configuration &C, double beta)
{
    const Worldline &wl = C.replica(0).worldline();
    int M = C.latticeSize();
    int N_pairs = 0;

    for (int site = 0; site < M; ++site)
    {
        int idx = wl.firstKink(site);
        while (idx != -1)
        {
            const Kink &k = wl[idx];
            if (k.partner != -1 && idx < k.partner)
                N_pairs++;
            idx = k.next;
        }
    }

    return -static_cast<double>(N_pairs) / beta;
}

// ------------------------------------------------------------
// Two-window symmetric kinetic estimator
// ------------------------------------------------------------
double kinetic_energy_two_window(const Configuration &C, double beta)
{
    const Worldline &wl = C.replica(0).worldline();

    double width = 0.10 * beta;
    double half_w = 0.5 * width;
    double centers[2] = {0.4 * beta, 0.6 * beta};

    int M = C.latticeSize();
    double E_kin_sum = 0.0;

    for (int w = 0; w < 2; ++w)
    {
        double tau_c = centers[w];
        int N_pairs_window = 0;

        for (int site = 0; site < M; ++site)
        {
            int idx = wl.firstKink(site);
            while (idx != -1)
            {
                const Kink &k = wl[idx];
                if (k.partner != -1 && idx < k.partner)
                {
                    double tau = k.tau;
                    if (tau > 0.0 && tau < beta &&
                        std::fabs(tau - tau_c) <= half_w)
                    {
                        N_pairs_window++;
                    }
                }
                idx = k.next;
            }
        }

        double E_kin_w = -static_cast<double>(N_pairs_window) / width;
        E_kin_sum += E_kin_w;
    }

    return 0.5 * E_kin_sum;
}

// ------------------------------------------------------------
int main()
{
    std::cout << "Running test_pigs_energy...\n\n";

    int N = 3;
    int L = 2;
    double t = 0.5;
    double U = 0.1; // <-- you can change this freely
    double mu = 0.0;

    std::cout << "beta\tE_diag_mid\tE_kin_mid\tE_kin_global\tE_tot_mid\tE_tot_global\n";

    std::vector<double> betas = {1.0, 2.0, 4.0, 6.0, 8.0, 16.0, 32.0};

    for (double beta : betas)
    {
        SimulationParameters params(beta);

        System sys(L, 1, t, U, mu);
        Lattice lat(L, 1);
        Configuration config(sys, lat, 1);

        BoseHubbardHamiltonian H(sys, lat);
        config.setHamiltonian(&H);

        std::vector<int> fock = {N, 0};
        config.initialize(fock);

        RNG rng(12345 + int(beta * 100));

        // Moves
        KinkAntikinkInsertionMC insertBulk(params);
        KinkAntikinkRemovalMC removeBulk(params);
        HopPairTimeShiftMC timeshift(params);
        BoundaryKinkPairInsertion0MC insertB0(params);
        BoundaryKinkPairRemoval0MC removeB0(params);
        BoundaryKinkPairInsertionBetaMC insertBbeta(params);
        BoundaryKinkPairRemovalBetaMC removeBbeta(params);

        // Move statistics
        MoveStats stats_insertBulk, stats_removeBulk, stats_timeshift;
        MoveStats stats_insertB0, stats_removeB0, stats_insertBbeta, stats_removeBbeta;

        // Long-run parameters
        int nSteps = 500000;
        int therm = 5000;
        int measureEvery = 20;

        double accum_diag = 0.0;
        double accum_kin_mid = 0.0;
        double accum_kin_global = 0.0;
        int count = 0;

        for (int step = 0; step < nSteps; ++step)
        {
            double r = rng.uniform();

            if (r < 0.15)
                stats_insertBulk.record(insertBulk.attempt(config, rng));
            else if (r < 0.30)
                stats_removeBulk.record(removeBulk.attempt(config, rng));
            else if (r < 0.45)
                stats_timeshift.record(timeshift.attempt(config, rng));
            else if (r < 0.60)
                stats_insertB0.record(insertB0.attempt(config, rng));
            else if (r < 0.75)
                stats_removeB0.record(removeB0.attempt(config, rng));
            else if (r < 0.875)
                stats_insertBbeta.record(insertBbeta.attempt(config, rng));
            else
                stats_removeBbeta.record(removeBbeta.attempt(config, rng));

            if (step >= therm && (step % measureEvery == 0))
            {
                accum_diag += diagonal_energy_mid(config, beta);
                accum_kin_mid += kinetic_energy_two_window(config, beta);
                accum_kin_global += kinetic_energy_global(config, beta);
                count++;
            }
        }

        double E_diag = accum_diag / count;
        double E_kin_mid = accum_kin_mid / count;
        double E_kin_global = accum_kin_global / count;

        double E_tot_mid = E_diag + E_kin_mid;
        double E_tot_global = E_diag + E_kin_global;

        std::cout << beta << "\t"
                  << E_diag << "\t"
                  << E_kin_mid << "\t"
                  << E_kin_global << "\t"
                  << E_tot_mid << "\t"
                  << E_tot_global << "\n";

        // Print move statistics
        std::cout << "\nMove statistics for beta=" << beta << "\n";
        stats_insertBulk.print("insertBulk");
        stats_removeBulk.print("removeBulk");
        stats_timeshift.print("timeshift");
        stats_insertB0.print("insertB0");
        stats_removeB0.print("removeB0");
        stats_insertBbeta.print("insertBbeta");
        stats_removeBbeta.print("removeBbeta");
        std::cout << "\n";
    }

    std::cout << "Done.\n";
    return 0;
}
