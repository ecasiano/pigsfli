#include "worldline.hpp"
#include "bose_hubbard_model.hpp"
#include "trial_wavefunction.hpp"
#include "estimator.hpp"
#include "simulation_parameters.hpp"

#include <iostream>
#include <cmath>

class DummyTrialWavefunction : public TrialWavefunction
{
public:
    double trial_wavefunction_ratio_one_site(
        const BoundaryFockState &, const BoundaryFockState &, int) const override
    {
        return 1.0;
    }

    double trial_wavefunction_ratio_two_sites(
        const BoundaryFockState &, const BoundaryFockState &, int, int) const override
    {
        return 1.0;
    }
};

class DummyEstimator : public Estimator
{
public:
    void accumulate(const ConfigurationView &, double) override {}
    void normalize() override {}
    void write_output(std::ostream &) const override {}
};

int main()
{
    int n_sites = 2;
    double t = 1.0;
    double U = 4.0;
    double mu = 0.0;

    std::vector<std::pair<int, int>> bonds = {{0, 1}};
    BoseHubbardModel model(n_sites, t, U, mu, bonds);

    DummyTrialWavefunction psi_trial;
    DummyEstimator est;

    SimulationParams params;
    params.n_sites = n_sites;
    params.beta = 4.0;
    params.p_insert_kink_antikink_pair = 0.5;
    params.p_delete_kink_antikink_pair = 0.5;

    Worldline wl(model, psi_trial, params, 4096);
    // Assume params.n_sites == 2
    wl.config().left_boundary.occupations = {1, 1};
    wl.config().right_boundary.occupations = {1, 1};

    std::cout << "start therm\n";
    int N_therm = 10; // was 1000
    for (int i = 0; i < N_therm; ++i)
        wl.attempt_move();
    std::cout << "done therm\n";

    std::cout << "start meas\n";
    int N_meas = 5000;
    double E_acc = 0.0;
    double E2_acc = 0.0;

    long ins_attempt = 0, ins_accept = 0;
    long del_attempt = 0, del_accept = 0;

    for (int i = 0; i < N_meas; ++i)
    {
        wl.attempt_move_with_stats(ins_attempt, ins_accept,
                                   del_attempt, del_accept);

        double E = wl.energy();
        E_acc += E;
        E2_acc += E * E;
    }
    std::cout << "done meas\n";

    double E_mean = E_acc / N_meas;
    double E_var = E2_acc / N_meas - E_mean * E_mean;
    double E_err = std::sqrt(E_var / N_meas);

    std::cout << "=== Bose-Hubbard Test ===\n";
    std::cout << "E_mean = " << E_mean << "\n";
    std::cout << "E_err  = " << E_err << "\n\n";

    std::cout << "Insertion attempts: " << ins_attempt
              << " accepted: " << ins_accept << "\n";
    std::cout << "Deletion attempts:  " << del_attempt
              << " accepted: " << del_accept << "\n";

    return 0;
}
