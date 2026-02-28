#pragma once

#include <vector>
#include <cstddef>
#include "kink.hpp"

class KinkPool
{
public:
    explicit KinkPool(std::size_t max_kinks);

    // allocate new active kink at the end, return its index
    int allocate();

    // move the last active kink into dst (does NOT change n_active_)
    void move_last_into(int dst);

    // decrement active count (drops last slot)
    void decrement_active();

    std::size_t n_active() const { return n_active_; }
    std::size_t capacity() const { return kinks_.size(); }

    Kink &operator[](int idx) { return kinks_[static_cast<std::size_t>(idx)]; }
    const Kink &operator[](int idx) const { return kinks_[static_cast<std::size_t>(idx)]; }

private:
    std::vector<Kink> kinks_;
    std::size_t n_active_ = 0;
};
