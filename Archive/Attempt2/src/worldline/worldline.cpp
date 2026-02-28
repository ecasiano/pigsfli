#include "../../include/worldline.hpp"
#include "../../include/trial_wavefunction.hpp"

#include <random>
#include <cassert>
#include <cmath>
#include <utility>
#include <algorithm>

// RNG
namespace
{
    std::mt19937_64 &rng_engine()
    {
        static std::mt19937_64 eng{std::random_device{}()};
        return eng;
    }
}

// Sampling helpers

double Worldline::random_uniform() const
{
    static std::uniform_real_distribution<double> dist(0.0, 1.0);
    return dist(rng_engine());
}

double Worldline::sample_tau(double beta) const
{
    return random_uniform() * beta;
}

std::pair<double, double> Worldline::sample_tau_pair(double beta) const
{
    double tau1 = random_uniform() * beta;
    double tau2 = random_uniform() * beta;
    if (tau2 < tau1)
        std::swap(tau1, tau2);
    return {tau1, tau2};
}

int Worldline::random_active_kink() const
{
    if (pool_.n_active() == 0)
        return -1;
    int idx = static_cast<int>(random_uniform() * pool_.n_active());
    if (idx == static_cast<int>(pool_.n_active()))
        idx--;
    return idx;
}

int Worldline::random_bond_index() const
{
    const auto &bonds = model_.bonds();
    if (bonds.empty())
        return -1;
    int idx = static_cast<int>(random_uniform() * bonds.size());
    if (idx == static_cast<int>(bonds.size()))
        idx--;
    return idx;
}

// Segment-based diagonal energy

double Worldline::diagonal_energy_from_segments() const
{
    double E = 0.0;

    for (int site = 0; site < params_.n_sites; ++site)
    {
        int n = config_.left_boundary.occupations[site];
        double last_tau = 0.0;

        int k = topo_.head(site);
        while (k != -1)
        {
            double tau = pool_[k].tau;
            double dt = tau - last_tau;
            if (dt < 0.0)
                dt += params_.beta;

            double e_site = model_.onsite_energy(site, n);
            E += e_site * dt;

            n += pool_[k].type;
            last_tau = tau;

            k = pool_[k].next;
        }

        double dt = params_.beta - last_tau;
        if (dt > 0.0)
        {
            double e_site = model_.onsite_energy(site, n);
            E += e_site * dt;
        }
    }

    return E / params_.beta;
}

// Constructor

Worldline::Worldline(const Model &model,
                     const TrialWavefunction &psi_trial,
                     const SimulationParams &params,
                     std::size_t max_kinks)
    : pool_(max_kinks),
      topo_(params.n_sites),
      model_(model),
      psi_trial_(psi_trial),
      params_(params)
{
    config_.left_boundary.occupations.assign(params_.n_sites, 0);
    config_.right_boundary.occupations.assign(params_.n_sites, 0);
}

// Main MC interface

void Worldline::attempt_move()
{
    choose_and_execute_move();
}

void Worldline::sweep(int n_moves)
{
    for (int i = 0; i < n_moves; ++i)
        attempt_move();
}

void Worldline::measure(Estimator &est) const
{
    est.accumulate(cfg_view_, 1.0);
}

// Invariants (unchanged)

bool Worldline::check_invariants() const
{
    const std::size_t n_active = pool_.n_active();

    auto in_range_or_minus1 = [&](int idx)
    {
        return idx == -1 || (idx >= 0 && static_cast<std::size_t>(idx) < n_active);
    };

    for (std::size_t k = 0; k < n_active; ++k)
    {
        const auto &K = pool_[k];

        if (!in_range_or_minus1(K.partner))
            return false;
        if (!in_range_or_minus1(K.pair_partner))
            return false;

        if (K.bond_id < 0 || K.bond_id >= (int)model_.bonds().size())
            return false;
    }

    for (std::size_t k = 0; k < n_active; ++k)
    {
        const auto &K = pool_[k];
        int p = K.partner;
        if (p != -1)
        {
            if (!in_range_or_minus1(p))
                return false;
            if (pool_[p].partner != (int)k)
                return false;
        }
    }

    for (std::size_t k = 0; k < n_active; ++k)
    {
        const auto &K = pool_[k];
        int pp = K.pair_partner;
        if (pp != -1)
        {
            if (!in_range_or_minus1(pp))
                return false;
            if (pool_[pp].pair_partner != (int)k)
                return false;
        }
    }

    std::vector<int> seen(n_active, 0);
    for (int site = 0; site < params_.n_sites; ++site)
    {
        int k = topo_.head(site);
        int prev = -1;
        while (k != -1)
        {
            if (k < 0 || (std::size_t)k >= n_active)
                return false;
            if (seen[k] != 0)
                return false;
            seen[k] = 1;

            if (prev != -1 && pool_[k].tau < pool_[prev].tau)
                return false;

            prev = k;
            k = pool_[k].next;
        }
    }

    for (std::size_t k = 0; k < n_active; ++k)
        if (seen[k] == 0)
            return false;

    for (int site = 0; site < params_.n_sites; ++site)
    {
        std::vector<Segment> segs;
        build_segments_for_site(site, segs);
        for (const auto &s : segs)
            if (s.occ < 0)
                return false;
    }

    return true;
}

// Dispatcher

void Worldline::choose_and_execute_move()
{
    double r = random_uniform();

    double p_insert = params_.p_insert_kink_antikink_pair;
    double p_delete = params_.p_delete_kink_antikink_pair;
    double norm = p_insert + p_delete;
    p_insert /= norm;
    p_delete /= norm;

    if (r < p_insert)
        move_insert_kink_antikink_pair();
    else
        move_delete_kink_antikink_pair();
}

// Segment helper

void Worldline::build_segments_for_site(int site, std::vector<Segment> &segs) const
{
    segs.clear();

    int n = config_.left_boundary.occupations[site];
    double last_tau = 0.0;

    int k = topo_.head(site);
    while (k != -1)
    {
        double tau = pool_[k].tau;

        if (tau > last_tau)
            segs.push_back({last_tau, tau, n});

        n += pool_[k].type;
        last_tau = tau;

        k = pool_[k].next;
    }

    if (last_tau < params_.beta)
        segs.push_back({last_tau, params_.beta, n});
}

// Insert move (unchanged)

void Worldline::move_insert_kink_antikink_pair()
{
    int bond_idx = random_bond_index();
    if (bond_idx == -1)
        return;

    const auto &bond = model_.bonds()[bond_idx];

    int direction = (random_uniform() < 0.5) ? +1 : -1;
    int site_dep = (direction > 0) ? bond.site_i : bond.site_j;
    int site_arr = (direction > 0) ? bond.site_j : bond.site_i;

    std::vector<Segment> seg_dep, seg_arr;
    build_segments_for_site(site_dep, seg_dep);
    build_segments_for_site(site_arr, seg_arr);

    struct Window
    {
        double tau_min, tau_max;
        int occ_dep, occ_arr;
    };
    std::vector<Window> windows;

    std::size_t i = 0, j = 0;
    while (i < seg_dep.size() && j < seg_arr.size())
    {
        double tau_min = std::max(seg_dep[i].tau_start, seg_arr[j].tau_start);
        double tau_max = std::min(seg_dep[i].tau_end, seg_arr[j].tau_end);

        if (tau_max > tau_min)
        {
            int occ_dep = seg_dep[i].occ;
            int occ_arr = seg_arr[j].occ;

            int occ_dep_after = occ_dep - 1;
            int occ_arr_after = occ_arr + 1;

            if (occ_dep_after >= 0 && occ_arr_after >= 0)
                windows.push_back({tau_min, tau_max, occ_dep, occ_arr});
        }

        if (seg_dep[i].tau_end < seg_arr[j].tau_end)
            ++i;
        else if (seg_dep[i].tau_end > seg_arr[j].tau_end)
            ++j;
        else
        {
            ++i;
            ++j;
        }
    }

    if (windows.empty())
        return;

    int N_bonds = (int)model_.bonds().size();
    int N_win = (int)windows.size();

    int w_idx = (int)(random_uniform() * N_win);
    if (w_idx == N_win)
        w_idx--;

    const auto &win = windows[w_idx];
    double delta_tau = win.tau_max - win.tau_min;

    double u1 = random_uniform();
    double u2 = random_uniform();
    double tau_a = win.tau_min + u1 * delta_tau;
    double tau_b = win.tau_min + u2 * delta_tau;
    double tau1 = std::min(tau_a, tau_b);
    double tau2 = std::max(tau_a, tau_b);

    int occ_dep_before = win.occ_dep;
    int occ_arr_before = win.occ_arr;
    int occ_dep_after = occ_dep_before - 1;
    int occ_arr_after = occ_arr_before + 1;

    double eps_alpha =
        model_.onsite_energy(site_dep, occ_dep_before) +
        model_.onsite_energy(site_arr, occ_arr_before);

    double eps_beta =
        model_.onsite_energy(site_dep, occ_dep_after) +
        model_.onsite_energy(site_arr, occ_arr_after);

    double delta_eps = eps_beta - eps_alpha;

    BoundaryFockState alpha_state = config_.left_boundary;
    BoundaryFockState beta_state = config_.left_boundary;

    alpha_state.occupations[site_dep] = occ_dep_before;
    alpha_state.occupations[site_arr] = occ_arr_before;

    beta_state.occupations[site_dep] = occ_dep_after;
    beta_state.occupations[site_arr] = occ_arr_after;

    double h = model_.offdiagonal_element(beta_state, alpha_state);
    double W_ratio = h * h * std::exp(-delta_eps * (tau2 - tau1));

    int N_pairs_prime = count_pairs_on_bond(bond_idx) + 1;

    double p_ins = params_.p_insert_kink_antikink_pair;
    double p_del = params_.p_delete_kink_antikink_pair;
    double norm = p_ins + p_del;
    p_ins /= norm;
    p_del /= norm;

    double T_forward = p_ins * (1.0 / N_bonds) * 0.5 * (1.0 / N_win) * (2.0 / (delta_tau * delta_tau));
    double T_backward = p_del * (1.0 / N_pairs_prime);
    double T_ratio = T_backward / T_forward;

    double A = std::min(1.0, W_ratio * T_ratio);
    if (random_uniform() >= A)
        return;

    int k1 = pool_.allocate();
    int k2 = pool_.allocate();
    int k3 = pool_.allocate();
    int k4 = pool_.allocate();

    pool_[k1].tau = tau1;
    pool_[k1].site = site_dep;
    pool_[k1].type = -1;
    pool_[k1].bond_id = bond_idx;
    pool_[k2].tau = tau1;
    pool_[k2].site = site_arr;
    pool_[k2].type = +1;
    pool_[k2].bond_id = bond_idx;
    pool_[k3].tau = tau2;
    pool_[k3].site = site_arr;
    pool_[k3].type = -1;
    pool_[k3].bond_id = bond_idx;
    pool_[k4].tau = tau2;
    pool_[k4].site = site_dep;
    pool_[k4].type = +1;
    pool_[k4].bond_id = bond_idx;

    pool_[k1].partner = k2;
    pool_[k2].partner = k1;
    pool_[k3].partner = k4;
    pool_[k4].partner = k3;

    pool_[k1].pair_partner = k4;
    pool_[k4].pair_partner = k1;
    pool_[k2].pair_partner = k3;
    pool_[k3].pair_partner = k2;

    topo_.insert_kink_on_site(site_dep, k1, pool_);
    topo_.insert_kink_on_site(site_arr, k2, pool_);
    topo_.insert_kink_on_site(site_arr, k3, pool_);
    topo_.insert_kink_on_site(site_dep, k4, pool_);
}

// Delete move with compaction + relabel

void Worldline::move_delete_kink_antikink_pair()
{
    int k = random_active_kink();
    if (k == -1)
        return;

    int p = pool_[k].partner;
    int pp = pool_[k].pair_partner;
    if (p == -1 || pp == -1)
        return;
    if (p < 0 || p >= (int)pool_.n_active())
        return;
    if (pp < 0 || pp >= (int)pool_.n_active())
        return;

    int q = pool_[pp].partner;
    if (q == -1)
        return;
    if (q < 0 || q >= (int)pool_.n_active())
        return;

    double tau_a = pool_[k].tau;
    double tau_b = pool_[pp].tau;

    int k1 = k, k2 = p, k3 = pp, k4 = q;
    if (tau_b < tau_a)
    {
        std::swap(tau_a, tau_b);
        std::swap(k1, k3);
        std::swap(k2, k4);
    }

    double tau1 = tau_a;
    double tau2 = tau_b;

    int bond_idx = pool_[k1].bond_id;
    const auto &bond = model_.bonds()[bond_idx];

    int site_dep = pool_[k1].site;
    int site_arr = pool_[k2].site;

    auto occ_at_tau = [&](int site, double tau) -> int
    {
        std::vector<Segment> segs;
        build_segments_for_site(site, segs);

        for (const auto &s : segs)
            if (tau >= s.tau_start && tau < s.tau_end)
                return s.occ;

        if (!segs.empty() && std::abs(tau - params_.beta) < 1e-12)
            return segs.back().occ;

        return config_.left_boundary.occupations[site];
    };

    int occ_dep_before = occ_at_tau(site_dep, tau1);
    int occ_arr_before = occ_at_tau(site_arr, tau1);

    int occ_dep_after = occ_dep_before - 1;
    int occ_arr_after = occ_arr_before + 1;
    if (occ_dep_after < 0 || occ_arr_after < 0)
        return;

    double eps_alpha =
        model_.onsite_energy(site_dep, occ_dep_before) +
        model_.onsite_energy(site_arr, occ_arr_before);

    double eps_beta =
        model_.onsite_energy(site_dep, occ_dep_after) +
        model_.onsite_energy(site_arr, occ_arr_after);

    double delta_eps = eps_beta - eps_alpha;

    BoundaryFockState alpha_state = config_.left_boundary;
    BoundaryFockState beta_state = config_.left_boundary;

    alpha_state.occupations[site_dep] = occ_dep_before;
    alpha_state.occupations[site_arr] = occ_arr_before;

    beta_state.occupations[site_dep] = occ_dep_after;
    beta_state.occupations[site_arr] = occ_arr_after;

    double h = model_.offdiagonal_element(beta_state, alpha_state);
    if (h == 0.0)
        return;

    double W_ratio = (1.0 / (h * h)) * std::exp(delta_eps * (tau2 - tau1));

    int N_bonds = (int)model_.bonds().size();
    int N_pairs = count_pairs_on_bond(bond_idx);
    if (N_pairs == 0)
        return;

    int direction = +1;
    if (site_dep == bond.site_j && site_arr == bond.site_i)
        direction = -1;

    int N_win_prime = count_windows_on_bond(bond_idx, direction);
    if (N_win_prime == 0)
        return;

    double p_ins = params_.p_insert_kink_antikink_pair;
    double p_del = params_.p_delete_kink_antikink_pair;
    double norm = p_ins + p_del;
    p_ins /= norm;
    p_del /= norm;

    double delta_tau = tau2 - tau1;

    double T_forward = p_del * (1.0 / N_pairs);
    double T_backward = p_ins * (1.0 / N_bonds) * 0.5 * (1.0 / N_win_prime) * (2.0 / (delta_tau * delta_tau));
    double T_ratio = T_backward / T_forward;

    double A = std::min(1.0, W_ratio * T_ratio);
    if (random_uniform() >= A)
        return;

    // unlink from topology
    topo_.remove_kink(k1, pool_);
    topo_.remove_kink(k2, pool_);
    topo_.remove_kink(k3, pool_);
    topo_.remove_kink(k4, pool_);

    // compact pool with relabeling
    int idxs[4] = {k1, k2, k3, k4};
    std::sort(std::begin(idxs), std::end(idxs));

    for (int t = 3; t >= 0; --t)
    {
        int old_idx = idxs[t];
        int last = (int)pool_.n_active() - 1;

        if (old_idx != last)
        {
            int site = pool_[last].site;
            pool_.move_last_into(old_idx);
            topo_.relabel_kink_index(site, last, old_idx, pool_);
        }

        pool_.decrement_active();
    }
}

// Stubs

void Worldline::move_timeshift_kink() {}

void Worldline::move_insert_kink_after_zero_edge() {}
void Worldline::move_delete_kink_after_zero_edge() {}
void Worldline::move_insert_kink_before_beta_edge() {}
void Worldline::move_delete_kink_before_beta_edge() {}

void Worldline::move_insert_worm() {}
void Worldline::move_delete_worm() {}
void Worldline::move_insert_antiworm() {}
void Worldline::move_delete_antiworm() {}
void Worldline::move_timeshift_head() {}
void Worldline::move_timeshift_tail() {}
void Worldline::move_insert_kink_before_head() {}
void Worldline::move_delete_kink_before_head() {}
void Worldline::move_insert_kink_before_tail() {}
void Worldline::move_delete_kink_before_tail() {}
void Worldline::move_insert_kink_after_head() {}
void Worldline::move_delete_kink_after_head() {}
void Worldline::move_insert_kink_after_tail() {}
void Worldline::move_delete_kink_after_tail() {}

void Worldline::move_insert_worm_zero_edge() {}
void Worldline::move_delete_worm_zero_edge() {}
void Worldline::move_insert_antiworm_zero_edge() {}
void Worldline::move_delete_antiworm_zero_edge() {}
void Worldline::move_insert_worm_beta_edge() {}
void Worldline::move_delete_worm_beta_edge() {}
void Worldline::move_insert_antiworm_beta_edge() {}
void Worldline::move_delete_antiworm_beta_edge() {}

void Worldline::move_insert_swap_kink() {}
void Worldline::move_delete_swap_kink() {}
void Worldline::move_timeshift_head_across_swap_kink() {}
void Worldline::move_timeshift_tail_across_swap_kink() {}

KapiLocalData Worldline::local_kapi_window_and_energies(int bond_idx,
                                                        int direction,
                                                        double tau_probe) const
{
    KapiLocalData out;

    const auto &bond = model_.bonds()[bond_idx];

    int site_dep = (direction > 0) ? bond.site_i : bond.site_j;
    int site_arr = (direction > 0) ? bond.site_j : bond.site_i;

    auto find_segment = [&](int site, double tau,
                            double &tau_prev, double &tau_next, int &occ_before)
    {
        int k = topo_.head(site);
        tau_prev = 0.0;
        tau_next = params_.beta;
        occ_before = config_.left_boundary.occupations[site];

        while (k != -1)
        {
            double t = pool_[k].tau;
            if (t > tau)
            {
                tau_next = t;
                return;
            }

            tau_prev = t;
            occ_before += pool_[k].type;
            k = pool_[k].next;
        }
    };

    double tau_prev_dep, tau_next_dep;
    double tau_prev_arr, tau_next_arr;
    int occ_dep_before, occ_arr_before;

    find_segment(site_dep, tau_probe, tau_prev_dep, tau_next_dep, occ_dep_before);
    find_segment(site_arr, tau_probe, tau_prev_arr, tau_next_arr, occ_arr_before);

    out.tau_min = std::max(tau_prev_dep, tau_prev_arr);
    out.tau_max = std::min(tau_next_dep, tau_next_arr);

    if (out.tau_max <= out.tau_min)
        return out;

    out.occ_dep_before = occ_dep_before;
    out.occ_arr_before = occ_arr_before;

    out.occ_dep_after = occ_dep_before - direction;
    out.occ_arr_after = occ_arr_before + direction;

    if (out.occ_dep_after < 0 || out.occ_arr_after < 0)
        return out;

    out.eps_alpha =
        model_.onsite_energy(site_dep, out.occ_dep_before) +
        model_.onsite_energy(site_arr, out.occ_arr_before);

    out.eps_beta =
        model_.onsite_energy(site_dep, out.occ_dep_after) +
        model_.onsite_energy(site_arr, out.occ_arr_after);

    BoundaryFockState before = config_.left_boundary;
    BoundaryFockState after = config_.left_boundary;

    before.occupations[site_dep] = out.occ_dep_before;
    before.occupations[site_arr] = out.occ_arr_before;

    after.occupations[site_dep] = out.occ_dep_after;
    after.occupations[site_arr] = out.occ_arr_after;

    out.h = model_.offdiagonal_element(before, after);
    out.allowed = true;
    return out;
}

int Worldline::count_pairs_on_bond(int bond_idx) const
{
    int count = 0;
    std::size_t n_active = pool_.n_active();

    for (std::size_t k = 0; k < n_active; ++k)
    {
        const auto &K = pool_[k];
        if (K.bond_id != bond_idx)
            continue;

        int p = K.partner;
        int pp = K.pair_partner;
        if (p == -1 || pp == -1)
            continue;

        int q = pool_[pp].partner;
        if (q == -1)
            continue;

        int k1 = (int)k;
        int k2 = p;
        int k3 = pp;
        int k4 = q;

        int min_idx = std::min(std::min(k1, k2), std::min(k3, k4));
        if (k1 == min_idx)
            ++count;
    }

    return count;
}

int Worldline::count_windows_on_bond(int bond_idx, int direction) const
{
    const auto &bond = model_.bonds()[bond_idx];

    int site_dep = (direction > 0) ? bond.site_i : bond.site_j;
    int site_arr = (direction > 0) ? bond.site_j : bond.site_i;

    std::vector<Segment> seg_dep, seg_arr;
    build_segments_for_site(site_dep, seg_dep);
    build_segments_for_site(site_arr, seg_arr);

    int N_win = 0;
    std::size_t i = 0, j = 0;

    while (i < seg_dep.size() && j < seg_arr.size())
    {
        double tau_min = std::max(seg_dep[i].tau_start, seg_arr[j].tau_start);
        double tau_max = std::min(seg_dep[i].tau_end, seg_arr[j].tau_end);

        if (tau_max > tau_min)
        {
            int occ_dep = seg_dep[i].occ;
            int occ_arr = seg_arr[j].occ;

            int occ_dep_after = occ_dep - direction;
            int occ_arr_after = occ_arr + direction;

            if (occ_dep_after >= 0 && occ_arr_after >= 0)
                ++N_win;
        }

        if (seg_dep[i].tau_end < seg_arr[j].tau_end)
            ++i;
        else if (seg_dep[i].tau_end > seg_arr[j].tau_end)
            ++j;
        else
        {
            ++i;
            ++j;
        }
    }

    return N_win;
}

void Worldline::attempt_move_with_stats(long &ins_attempt, long &ins_accept,
                                        long &del_attempt, long &del_accept)
{
    double r = random_uniform();

    double p_insert = params_.p_insert_kink_antikink_pair;
    double p_delete = params_.p_delete_kink_antikink_pair;
    double norm = p_insert + p_delete;
    p_insert /= norm;
    p_delete /= norm;

    if (r < p_insert)
    {
        ++ins_attempt;
        std::size_t before = pool_.n_active();
        move_insert_kink_antikink_pair();
        if (pool_.n_active() > before)
            ++ins_accept;
    }
    else
    {
        ++del_attempt;
        std::size_t before = pool_.n_active();
        move_delete_kink_antikink_pair();
        if (pool_.n_active() < before)
            ++del_accept;
    }
}
