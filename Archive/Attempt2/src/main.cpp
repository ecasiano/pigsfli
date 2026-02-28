#include <iostream>

#include "../include/worldline.hpp"
#include "../include/models/bose_hubbard_model.hpp"
#include "../include/trial_wavefunctions/flat.hpp"
#include "../include/simulation_parameters.hpp"

int main()
{
    // 1. Build model
    Model model;
    // fill model.bonds() with one bond (0,1)
    // set t, U, mu, n_max, etc.

    TrialWavefunction psi_trial; // can be trivial for now

    SimulationParams params;
    params.n_sites = 2;
    params.beta = 4.0;
    params.p_insert_kink_antikink_pair = 0.5;
    params.p_delete_kink_antikink_pair = 0.5;

    Worldline wl(model, psi_trial, params, /*max_kinks=*/1024);

    // set boundary occupations, e.g. (1,0)
    // wl.config_.left_boundary.occupations = {1,0};
    // wl.config_.right_boundary.occupations = {1,0}; // or same for now

    Estimator est; // implement diagonal energy estimator using cfg_view_

    // 2. Thermalization
    for (int i = 0; i < 100000; ++i)
        wl.attempt_move();

    // 3. Measurement
    const int N_meas = 100000;
    double E_acc = 0.0;
    double E2_acc = 0.0;

    for (int i = 0; i < N_meas; ++i)
    {
        wl.attempt_move();

        double E = wl.diagonal_energy_from_segments();
        E_acc += E;
        E2_acc += E * E;

        // optional: check invariants occasionally
        // if ((i % 1000) == 0 && !wl.check_invariants()) { ... }
    }

    double E_mean = E_acc / N_meas;
    double E_var = E2_acc / N_meas - E_mean * E_mean;

    std::cout << "E_mean = " << E_mean << "\n";
    std::cout << "E_std  = " << std::sqrt(E_var / N_meas) << "\n";

    return 0;
}
