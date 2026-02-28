#pragma once

#include <vector>
#include "kink_pool.hpp"

class WorldlineTopology
{
public:
    explicit WorldlineTopology(int n_sites);

    // head of per-site list
    int head(int site) const { return heads_[site]; }

    // insert kink into site's time-ordered list
    void insert_kink_on_site(int site, int kink_index, KinkPool &pool);

    // remove kink from its site's list
    void remove_kink(int kink_index, KinkPool &pool);

    // after pool_.move_last_into(old_idx) on given site,
    // relabel all references to old_idx -> new_idx in that site's list
    void relabel_kink_index(int site, int old_idx, int new_idx, KinkPool &pool);

private:
    std::vector<int> heads_; // per-site head indices
};
