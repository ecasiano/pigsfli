#pragma once

struct SimulationParams
{
    // ------------------------------------------------------------
    // System size
    // ------------------------------------------------------------
    int n_sites = 0;
    int n_replicas = 1;

    // ------------------------------------------------------------
    // Imaginary-time extent
    // ------------------------------------------------------------
    double beta = 0.0;

    // ------------------------------------------------------------
    // Monte Carlo control
    // ------------------------------------------------------------
    int n_thermalization = 0;
    int n_measurements = 0;

    // ------------------------------------------------------------
    // Move probabilities (explicit, descriptive names)
    // ------------------------------------------------------------
    double p_insert_kink_antikink_pair = 0.5;
    double p_delete_kink_antikink_pair = 0.5;

    // ------------------------------------------------------------
    // Time window for kink–antikink insertion
    // (You will later replace this with per-segment sampling)
    // ------------------------------------------------------------
    double tau_min = 0.0;
    double tau_max = 0.0;
};
