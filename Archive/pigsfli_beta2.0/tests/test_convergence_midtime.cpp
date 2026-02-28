#include "../src/pimc.hpp"
#include "../src/moves_mc.hpp"
#include "../src/moves_timeshift.hpp"
#include "../src/estimators.hpp"
#include "../src/ed_bose_hubbard_2site.hpp"
#include "../src/moves_boundary.hpp"

#include <iostream>
#include <vector>
#include <cmath>

using namespace pimc;

// ------------------------------------------------------------
// Diagonal energy at tau = beta/2
// ------------------------------------------------------------
double diagonal_energy_mid(const Configuration &C, double beta)
{
    const System &sys = C.system();
    const Lattice &lat = C.lattice();
    const Worldline &wl = C.replica(0).worldline();

    double tau_mid = 0.5 * beta;
    int M = lat.size();

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

    return E_diag;
}

// ------------------------------------------------------------
// Local-window kinetic energy around tau = beta/2
// ------------------------------------------------------------
double kinetic_energy_mid_local(const Configuration &C, double beta)
{
    const System &sys = C.system();
    const Lattice &lat = C.lattice();
    const Worldline &wl = C.replica(0).worldline();

    double tau_mid = 0.5 * beta;
    double window = 0.5;
    if (window > beta)
        window = beta;
    double half_w = 0.5 * window;

    int M = lat.size();
    int N_window = 0;

    for (int site = 0; site < M; ++site)
    {
        int idx = wl.firstKink(site);
        while (idx != -1)
        {
            double tau_k = wl[idx].tau;

            double dt = std::fabs(tau_k - tau_mid);
            if (dt > 0.5 * beta)
                dt = beta - dt;

            if (dt <= half_w)
                N_window++;

            idx = wl[idx].next;
        }
    }

    double E_kin = -sys.t() * static_cast<double>(N_window) / window;
    return E_kin;
}

// ------------------------------------------------------------
// Count physical hops
// ------------------------------------------------------------
int count_hops(const Configuration &C)
{
    const Worldline &wl = C.replica(0).worldline();
    int M = C.lattice().size();
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
    return count;
}

// ------------------------------------------------------------
// Mid-time occupations
// ------------------------------------------------------------
std::pair<double, double> midtime_occupations(const Configuration &C, double beta)
{
    const Worldline &wl = C.replica(0).worldline();
    double tau_mid = 0.5 * beta;
    int n0 = wl.occupationAt(0, tau_mid);
    int n1 = wl.occupationAt(1, tau_mid);
    return {static_cast<double>(n0), static_cast<double>(n1)};
}

// ------------------------------------------------------------
// Main test
// ------------------------------------------------------------
int main()
{
    std::cout << "Running test_convergence_midtime...\n";

    int N = 2;
    int L = 2;
    double t = 0.5;
    double U = 1.0;
    double mu = 0.0;

    EDResult2Site ed = exactDiagonalization2Site(N, t, U, mu);
    double E_exact = ed.E0;

    std::cout << "Exact ground state energy = " << E_exact << "\n\n";
    std::cout << "beta\tE_total\tE_diag_mid\tE_kin_mid\tE_mid_total\t"
              << "Nk_avg\tNk_density\tn0_mid\tn1_mid\tE_kin_global\n";

    std::vector<double> betas = {0.5, 1.0, 2.0, 4.0, 8.0, 16.0};

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

        KinkAntikinkInsertionMC insertMove(params);
        KinkAntikinkRemovalMC removeMove(params);
        HopPairTimeShiftMC timeshiftMove(params);
        BoundaryInsertionMC boundaryInsert(params);
        BoundaryRemovalMC boundaryRemove(params);

        TotalEnergyEstimator Etot(beta);

        int nSteps = 2000000;
        int therm = 200000;

        double accumE_tot = 0.0;
        double accumE_diag_mid = 0.0;
        double accumE_kin_mid = 0.0;
        double accumNk = 0.0;
        double accum_n0_mid = 0.0;
        double accum_n1_mid = 0.0;
        double accumNhop = 0.0;
        int count = 0;

        for (int step = 0; step < nSteps; ++step)
        {
            double r = rng.uniform();

            double r = rng.uniform();

            if (r < 0.3)
                insertMove.attempt(config, rng);
            else if (r < 0.6)
                removeMove.attempt(config, rng);
            else if (r < 0.75)
                timeshiftMove.attempt(config, rng);
            else if (r < 0.875)
                boundaryInsert.attempt(config, rng);
            else
                boundaryRemove.attempt(config, rng);

            if (step >= therm)
            {
                accumE_tot += Etot.measure(config);
                accumE_diag_mid += diagonal_energy_mid(config, beta);
                accumE_kin_mid += kinetic_energy_mid_local(config, beta);

                accumNk += average_kink_count(config);
                auto [n0_mid, n1_mid] = midtime_occupations(config, beta);
                accum_n0_mid += n0_mid;
                accum_n1_mid += n1_mid;

                accumNhop += count_hops(config);

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

        double Nhop_avg = accumNhop / count;
        double E_kin_global = -sys.t() * Nhop_avg / beta;

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
