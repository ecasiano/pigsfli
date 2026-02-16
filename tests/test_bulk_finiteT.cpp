#include "../src/pimc.hpp"
#include "../src/moves_mc.hpp"
#include "../src/moves_timeshift.hpp"
#include "../src/estimators.hpp"

#include <iostream>
#include <vector>
#include <cmath>

using namespace pimc;

// ------------------------------------------------------------
// Finite-T diagonal energy (exact for 2-site Bose-Hubbard)
// ------------------------------------------------------------
double exact_energy_finiteT(double beta, double t, double U, double mu)
{
    // Hilbert space for 2 bosons on 2 sites:
    // |2,0>, |1,1>, |0,2>
    double E20 = 0.5 * U * 2 * (2 - 1) - mu * 2; // = U - 2mu
    double E11 = 0.0 - mu * 2;                   // = -2mu
    double E02 = E20;

    // Hopping couples |2,0> <-> |1,1> <-> |0,2>
    // Exact diagonalization for 2-site BH with N=2:
    // Energies: E0 = -sqrt(2)*t + U/2, E1 = U, E2 = +sqrt(2)*t + U/2
    double E0 = 0.5 * U - std::sqrt(2.0) * t;
    double E1 = U;
    double E2 = 0.5 * U + std::sqrt(2.0) * t;

    double Z = std::exp(-beta * E0) + std::exp(-beta * E1) + std::exp(-beta * E2);
    double E = (E0 * std::exp(-beta * E0) +
                E1 * std::exp(-beta * E1) +
                E2 * std::exp(-beta * E2)) /
               Z;

    return E;
}

// ------------------------------------------------------------
// Main test
// ------------------------------------------------------------
int main()
{
    std::cout << "Running test_bulk_finiteT...\n";

    int N = 2;
    int L = 2;
    double t = 0.5;
    double U = 1.0;
    double mu = 0.0;

    double beta = 1.0; // finite-T test
    double E_exact = exact_energy_finiteT(beta, t, U, mu);

    std::cout << "Exact finite-T energy at beta=" << beta << " is " << E_exact << "\n\n";

    // ------------------------------------------------------------
    // Build system
    // ------------------------------------------------------------
    SimulationParameters params(beta);

    System sys(L, 1, t, U, mu);
    Lattice lat(L, 1);
    Configuration config(sys, lat, 1);

    BoseHubbardHamiltonian H(sys, lat);
    config.setHamiltonian(&H);

    // Initial Fock state: |2,0>
    std::vector<int> fock = {N, 0};
    config.initialize(fock);

    RNG rng(12345);

    // ------------------------------------------------------------
    // Moves: BULK ONLY
    // ------------------------------------------------------------
    KinkAntikinkInsertionMC insertMove(params);
    KinkAntikinkRemovalMC removeMove(params);
    HopPairTimeShiftMC timeshiftMove(params); // fixed version

    TotalEnergyEstimator Etot(beta);

    int nSteps = 200000;
    int therm = 20000;

    double accumE = 0.0;
    int count = 0;

    for (int step = 0; step < nSteps; ++step)
    {
        double r = rng.uniform();

        if (r < 0.33)
            insertMove.attempt(config, rng);
        else if (r < 0.66)
            removeMove.attempt(config, rng);
        else
            timeshiftMove.attempt(config, rng);

        if (step >= therm)
        {
            accumE += Etot.measure(config);
            count++;
        }
    }

    double E_mc = accumE / count;

    std::cout << "MC energy = " << E_mc << "\n";
    std::cout << "Exact     = " << E_exact << "\n";

    std::cout << "\nDone.\n";
    return 0;
}
