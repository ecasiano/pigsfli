#pragma once

#include <utility>
#include <vector>

#include "model.hpp"
#include "trial_wavefunction.hpp"
#include "estimator.hpp"
#include "simulation_parameters.hpp"
#include "configuration_view.hpp"

#include "worldline_internal/kink_pool.hpp"
#include "worldline_internal/worldline_topology.hpp"

struct KapiLocalData
{
    double tau_min = 0.0;
    double tau_max = 0.0;

    int occ_dep_before = 0;
    int occ_arr_before = 0;
    int occ_dep_after = 0;
    int occ_arr_after = 0;

    double eps_alpha = 0.0;
    double eps_beta = 0.0;
    double h = 0.0;

    bool allowed = false;
};

class Worldline
{
public:
    Worldline(const Model &model,
              const TrialWavefunction &psi_trial,
              const SimulationParams &params,
              std::size_t max_kinks);

    void attempt_move();
    void sweep(int n_moves);

    void measure(Estimator &est) const;

    const ConfigurationView &view() const { return cfg_view_; }

    bool check_invariants() const;

    // ------------------------------------------------------------
    // Boundary setter (canonical N via occupations)
    // ------------------------------------------------------------
    void set_boundary_occupations(const std::vector<int> &occ)
    {
        config_.left_boundary.occupations = occ;
        config_.right_boundary.occupations = occ;
    }

    // ------------------------------------------------------------
    // Move wrapper with stats
    // ------------------------------------------------------------
    void attempt_move_with_stats(long &ins_attempt, long &ins_accept,
                                 long &del_attempt, long &del_accept);

    // ------------------------------------------------------------
    // Public energy wrapper (diagonal energy)
    // ------------------------------------------------------------
    double energy() const { return diagonal_energy_from_segments(); }

    // ------------------------------------------------------------
    // Total particle number diagnostic (from segments at tau = 0+)
    // ------------------------------------------------------------
    int total_particle_number_from_segments() const;

private:
    void choose_and_execute_move();

    // Non-worm moves (bulk)
    void move_insert_kink_antikink_pair();
    void move_delete_kink_antikink_pair();

    void move_timeshift_kink();

    // Edge, worm, swap moves (stubs)
    void move_insert_kink_after_zero_edge();
    void move_delete_kink_after_zero_edge();
    void move_insert_kink_before_beta_edge();
    void move_delete_kink_before_beta_edge();

    void move_insert_worm();
    void move_delete_worm();
    void move_insert_antiworm();
    void move_delete_antiworm();
    void move_timeshift_head();
    void move_timeshift_tail();
    void move_insert_kink_before_head();
    void move_delete_kink_before_head();
    void move_insert_kink_before_tail();
    void move_delete_kink_before_tail();
    void move_insert_kink_after_head();
    void move_delete_kink_after_head();
    void move_insert_kink_after_tail();
    void move_delete_kink_after_tail();

    void move_insert_worm_zero_edge();
    void move_delete_worm_zero_edge();
    void move_insert_antiworm_zero_edge();
    void move_delete_antiworm_zero_edge();
    void move_insert_worm_beta_edge();
    void move_delete_worm_beta_edge();
    void move_insert_antiworm_beta_edge();
    void move_delete_antiworm_beta_edge();

    void move_insert_swap_kink();
    void move_delete_swap_kink();
    void move_timeshift_head_across_swap_kink();
    void move_timeshift_tail_across_swap_kink();

    // Helpers
    double random_uniform() const;
    double sample_tau(double beta) const;
    std::pair<double, double> sample_tau_pair(double beta) const;

    int random_active_kink() const;
    int random_bond_index() const;

    void delete_kink_pair(int k1, int k2);

    double diagonal_energy_from_segments() const;

    KapiLocalData local_kapi_window_and_energies(int bond_idx,
                                                 int direction,
                                                 double tau_probe) const;

    int count_pairs_on_bond(int bond_idx) const;
    int count_windows_on_bond(int bond_idx, int direction) const;

    struct Segment
    {
        double tau_start;
        double tau_end;
        int occ;
    };

    void build_segments_for_site(int site, std::vector<Segment> &segs) const;

private:
    KinkPool pool_;
    WorldlineTopology topo_;

    const Model &model_;
    const TrialWavefunction &psi_trial_;
    SimulationParams params_;

    ConfigurationView cfg_view_;

    struct
    {
        BoundaryFockState left_boundary;
        BoundaryFockState right_boundary;
    } config_;
};
