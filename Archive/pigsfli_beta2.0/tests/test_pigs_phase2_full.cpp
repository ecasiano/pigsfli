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
double kinetic_energy_mid_local(const Configuration &C, double beta)
{
    const System &sys = C.system();
    const Worldline &wl = C.replica(0).worldline();

    double tau_mid = 0.5 * beta;
    double window = 0.5;
    double half_w = 0.5 * window;

    int M = C.latticeSize();
    int N_window = 0;

    for (int site = 0; site < M; ++site)
    {
        int idx = wl.firstKink(site);
        while (idx != -1)
        {
            double tau = wl[idx].tau;

            if (tau == 0.0 || tau == beta)
            {
                idx = wl[idx].next;
                continue;
            }

            double dt = std::fabs(tau - tau_mid);
            if (dt <= half_w)
                N_window++;

            idx = wl[idx].next;
        }
    }

    return -sys.t() * (double)N_window / window;
}

// ------------------------------------------------------------
double average_kink_count(const Configuration &C)
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
                count++;
            idx = k.next;
        }
    }
    return (double)count;
}

// ------------------------------------------------------------
std::pair<double, double> midtime_occupations(const Configuration &C, double beta)
{
    const Worldline &wl = C.replica(0).worldline();
    double tau_mid = 0.5 * beta;
    return {
        (double)wl.occupationAt(0, tau_mid),
        (double)wl.occupationAt(1, tau_mid)};
}

// ------------------------------------------------------------
int main()
{
    std::cout << "Running test_pigs_phase2_full...\n";

    int N = 2;
    int L = 2;
    double t = 0.5;
    double U = 1.0;
    double mu = 0.0;

    double E_exact = -0.618034;
    std::cout << "Exact ground state energy = " << E_exact << "\n\n";

    std::cout << "beta\tE_total\tE_diag_mid\tE_kin_mid\tE_mid_total\t"
              << "Nk_avg\tNk_density\tn0_mid\tn1_mid\tE_kin_global\n";

    std::vector<double> betas = {0.5, 1.0, 2.0, 4.0, 8.0};

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

        TotalEnergyEstimator Etot(beta);

        int nSteps = 2000000;
        int therm = 100000;

        double accumE_tot = 0.0;
        double accumE_diag_mid = 0.0;
        double accumE_kin_mid = 0.0;
        double accumNk = 0.0;
        double accum_n0_mid = 0.0;
        double accum_n1_mid = 0.0;
        int count = 0;

        for (int step = 0; step < nSteps; ++step)
        {
            double r = rng.uniform();

            if (r < 0.15)
                insertBulk.attempt(config, rng);
            else if (r < 0.30)
                removeBulk.attempt(config, rng);
            else if (r < 0.45)
                timeshift.attempt(config, rng);
            else if (r < 0.60)
                insertB0.attempt(config, rng);
            else if (r < 0.75)
                removeB0.attempt(config, rng);
            else if (r < 0.875)
                insertBbeta.attempt(config, rng);
            else
                removeBbeta.attempt(config, rng);

            if (step >= therm)
            {
                accumE_tot += Etot.measure(config);
                accumE_diag_mid += diagonal_energy_mid(config, beta);
                accumE_kin_mid += kinetic_energy_mid_local(config, beta);

                accumNk += average_kink_count(config);
                auto [n0_mid, n1_mid] = midtime_occupations(config, beta);
                accum_n0_mid += n0_mid;
                accum_n1_mid += n1_mid;

                count++;
            }
        }

        double E_mc_tot = accumE_tot / count;
        double E_mc_diag_mid = accumE_diag_mid / count;
        double E_mc_kin_mid = accumE_kin_mid / count;
        double E_mc_mid_tot = E_mc_diag_mid + E_mc_kin_mid;

        double Nk_avg = accumNk / count;
        double Nk_density = Nk_avg / beta;
        double n0_mid_avg = accum_n0_mid / count;
        double n1_mid_avg = accum_n1_mid / count;

        double E_kin_global = -sys.t() * Nk_avg / beta;

        std::cout << beta << "\t"
                  << E_mc_tot << "\t"
                  << E_mc_diag_mid << "\t"
                  << E_mc_kin_mid << "\t"
                  << E_mc_mid_tot << "\t"
                  << Nk_avg << "\t"
                  << Nk_density << "\t"
                  << n0_mid_avg << "\t"
                  << n1_mid_avg << "\t"
                  << E_kin_global << "\n";
    }

    std::cout << "\nDone.\n";
    return 0;
}
