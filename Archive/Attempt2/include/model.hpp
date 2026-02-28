#pragma once

#include <vector>

// Forward declaration; full type is in configuration_view.hpp
struct BoundaryFockState;

// ------------------------------------------------------------
// Bond structure for model-dependent hopping
// ------------------------------------------------------------
struct Bond
{
    int site_i;
    int site_j;
    double t;
};

// ------------------------------------------------------------
// Abstract base class for all models
// ------------------------------------------------------------
class Model
{
public:
    virtual ~Model() = default;

    // Bonds (nearest neighbor, next-nearest, etc.)
    const std::vector<Bond> &bonds() const { return bonds_; }

    // Onsite diagonal energy for a given site and occupation.
    virtual double onsite_energy(int site, int n) const = 0;

    // Full diagonal energy for a boundary state
    virtual double diagonal_energy(const BoundaryFockState &state) const = 0;

    // Off-diagonal matrix element
    virtual double offdiagonal_element(const BoundaryFockState &before,
                                       const BoundaryFockState &after) const = 0;

protected:
    std::vector<Bond> bonds_;
};
