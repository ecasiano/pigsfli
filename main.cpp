#include "pimc.hpp"
#include "export.hpp"
#include <iostream>

using namespace pimc;

/**
 * Entry point for basic initialization and worldline inspection.
 *
 * This program:
 *   1. Builds the static System and Lattice
 *   2. Constructs a Configuration with one replica
 *   3. Generates a random Fock state with N bosons
 *   4. Initializes the worldline from that Fock state
 *   5. Prints and exports the resulting worldline
 *
 * This is a structural test, not a Monte Carlo simulation.
 */

int main()
{
    int L = 2;
    int D = 2;
    int num_replicas = 1;

    pimc::RNG rng(1234);

    // Build system + lattice
    pimc::System sys(L, D, 1.0, 1.0, 0.0);
    pimc::Lattice lat(L, D);

    // Build configuration
    pimc::Configuration config(sys, lat, num_replicas);

    // --- NEW: attach Hamiltonian ---
    BoseHubbardHamiltonian H(sys, lat);
    config.setHamiltonian(&H);

    int M = config.latticeSize();
    int N = 3;

    // Random initial Fock state
    auto fock = random_fock_state(M, N, rng);

    // Initialize all replicas
    config.initialize(fock);

    // Print worldline for replica 0
    print_worldline(config.replica(0).worldline(), M);

    // Export for Python
    export_worldline(config.replica(0).worldline(), M, 1.0, "worldline.txt");

    std::cout << "Exported worldline to worldline.txt\n";
    return 0;
}
