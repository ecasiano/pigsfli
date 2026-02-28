#include "../../include/worldline_internal/kink_pool.hpp"
#include <cassert>

KinkPool::KinkPool(std::size_t max_kinks)
    : kinks_(max_kinks), n_active_(0)
{
}

int KinkPool::allocate()
{
    assert(n_active_ < kinks_.size());
    int idx = static_cast<int>(n_active_);
    ++n_active_;
    return idx;
}

void KinkPool::move_last_into(int dst)
{
    assert(n_active_ > 0);
    std::size_t last = n_active_ - 1;
    std::size_t d = static_cast<std::size_t>(dst);
    assert(d < n_active_);
    if (d != last)
        kinks_[d] = kinks_[last];
}

void KinkPool::decrement_active()
{
    assert(n_active_ > 0);
    --n_active_;
}
