#pragma once
#include "../model.hpp"
#include <vector>
#include <utility>
#include <cmath>

class BoseHubbardModel : public Model
{
public:
    BoseHubbardModel(int n_sites,
                     double t,
                     double U,
                     double mu,
                     const std::vector<std::pair<int, int>> &bonds)
        : n_sites_(n_sites),
          t_(t),
          U_(U),
          mu_(mu)
    {
        // Fill base-class bonds_
        for (auto &b : bonds)
            bonds_.push_back({b.first, b.second, t_});
    }

    // ------------------------------------------------------------
    // Onsite energy U/2 * n(n-1) - mu * n
    // ------------------------------------------------------------
    double onsite_energy(int site, int n) const override
    {
        return 0.5 * U_ * n * (n - 1) - mu_ * n;
    }

    // ------------------------------------------------------------
    // Full diagonal energy (optional override)
    // ------------------------------------------------------------
    double diagonal_energy(const BoundaryFockState &state) const override
    {
        double E = 0.0;
        for (int i = 0; i < n_sites_; ++i)
        {
            int n = state.occupations[i];
            E += onsite_energy(i, n);
        }
        return E;
    }

    // ------------------------------------------------------------
    // Off-diagonal matrix element for a hop j → i
    // ------------------------------------------------------------
    double offdiagonal_element(const BoundaryFockState &before,
                               const BoundaryFockState &after) const override
    {
        int from = -1;
        int to = -1;

        for (int i = 0; i < n_sites_; ++i)
        {
            int db = after.occupations[i] - before.occupations[i];
            if (db == +1)
                to = i;
            if (db == -1)
                from = i;
        }

        if (from < 0 || to < 0)
            return 0.0;

        // Check if (from,to) is a valid bond
        bool connected = false;
        for (auto &b : bonds_)
        {
            if ((b.site_i == from && b.site_j == to) ||
                (b.site_i == to && b.site_j == from))
            {
                connected = true;
                break;
            }
        }

        if (!connected)
            return 0.0;

        int n_from_before = before.occupations[from];
        int n_to_after = after.occupations[to];

        return -t_ * std::sqrt(n_from_before * n_to_after);
    }

private:
    int n_sites_;
    double t_;
    double U_;
    double mu_;
};
