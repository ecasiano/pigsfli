// ============================================================================
//  trial_wavefunction.hpp
//  ----------------------
//  Defines the abstract TrialWavefunction interface.
//
//  In PIGS, the boundary states |α₀⟩ and |α_β⟩ come from a trial distribution.
//  Boundary moves (zero-edge and beta-edge) require evaluating:
//
//      • log Ψ_T(α)
//      • log-amplitude ratios log[Ψ_T(β)/Ψ_T(α)]
//
//  The Worldline engine NEVER hard-codes trial-state logic. It only calls this
//  interface when computing Metropolis ratios for boundary updates.
//
//  This allows:
//
//      • Gutzwiller trial states
//      • Jastrow factors
//      • neural-network trial states
//      • DMRG boundary states
//
//  …all without modifying the worldline engine.
// ============================================================================

#pragma once
#include <vector>

// -----------------------------------------------------------------------------
// BoundaryFockState
// Represents the full boundary Fock state at tau = 0 or tau = beta.
// The user never constructs full wavefunctions; only ratios are computed.
// -----------------------------------------------------------------------------
struct BoundaryFockState
{
    std::vector<int> occupations; // one integer per site
};

// -----------------------------------------------------------------------------
// TrialWavefunction (interface)
// The user provides analytic expressions for ratios of trial-state coefficients.
// No full wavefunction is ever constructed.
//
// Two ratio functions are required:
//   • trial_wavefunction_ratio_one_site(before, after, site)
//   • trial_wavefunction_ratio_two_sites(before, after, site_i, site_j)
//
// These correspond to boundary updates that modify one or two sites.
// -----------------------------------------------------------------------------
class TrialWavefunction
{
public:
    virtual ~TrialWavefunction() = default;

    // Ratio for boundary updates that modify ONE site.
    virtual double trial_wavefunction_ratio_one_site(
        const BoundaryFockState &before,
        const BoundaryFockState &after,
        int site) const = 0;

    // Ratio for boundary updates that modify TWO sites.
    virtual double trial_wavefunction_ratio_two_sites(
        const BoundaryFockState &before,
        const BoundaryFockState &after,
        int site_i,
        int site_j) const = 0;
};
