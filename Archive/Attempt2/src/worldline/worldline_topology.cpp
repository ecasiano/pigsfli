#include "../../include/worldline_internal/worldline_topology.hpp"
#include <cassert>

WorldlineTopology::WorldlineTopology(int n_sites)
    : heads_(n_sites, -1)
{
}

void WorldlineTopology::insert_kink_on_site(int site, int kink_index, KinkPool &pool)
{
    assert(site >= 0 && site < (int)heads_.size());

    int &head = heads_[site];
    int prev = -1;
    int curr = head;

    double tau_new = pool[kink_index].tau;

    while (curr != -1 && pool[curr].tau < tau_new)
    {
        prev = curr;
        curr = pool[curr].next;
    }

    pool[kink_index].site = site;
    pool[kink_index].prev = prev;
    pool[kink_index].next = curr;

    if (prev != -1)
        pool[prev].next = kink_index;
    else
        head = kink_index;

    if (curr != -1)
        pool[curr].prev = kink_index;
}

void WorldlineTopology::remove_kink(int kink_index, KinkPool &pool)
{
    int site = pool[kink_index].site;
    assert(site >= 0 && site < (int)heads_.size());

    int prev = pool[kink_index].prev;
    int next = pool[kink_index].next;

    if (prev != -1)
        pool[prev].next = next;
    else
        heads_[site] = next;

    if (next != -1)
        pool[next].prev = prev;

    pool[kink_index].prev = -1;
    pool[kink_index].next = -1;
}

void WorldlineTopology::relabel_kink_index(int site, int old_idx, int new_idx, KinkPool &pool)
{
    assert(site >= 0 && site < (int)heads_.size());

    // fix head if needed
    if (heads_[site] == old_idx)
        heads_[site] = new_idx;

    int k = heads_[site];
    while (k != -1)
    {
        Kink &K = pool[k];
        if (K.next == old_idx)
            K.next = new_idx;
        if (K.prev == old_idx)
            K.prev = new_idx;
        k = K.next;
    }
}
