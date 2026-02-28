#pragma once
#include "../trial_wavefunction.hpp"

// -----------------------------------------------------------------------------
// FlatTrialWavefunction
// Psi_T(alpha) = constant for all Fock states.
// All ratios = 1.
// -----------------------------------------------------------------------------
class FlatTrialWavefunction : public TrialWavefunction
{
public:
    double trial_wavefunction_ratio_one_site(
        const BoundaryFockState &before,
        const BoundaryFockState &after,
        int site) const override
    {
        (void)before;
        (void)after;
        (void)site;
        return 1.0;
    }

    double trial_wavefunction_ratio_two_sites(
        const BoundaryFockState &before,
        const BoundaryFockState &after,
        int site_i,
        int site_j) const override
    {
        (void)before;
        (void)after;
        (void)site_i;
        (void)site_j;
        return 1.0;
    }
};
