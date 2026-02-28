// ============================================================================
//  configuration_view.hpp
//  ----------------------
//  Defines read-only views of the worldline configuration.
//
//  Estimators receive a ConfigurationView, NOT direct access to kinks.
//  This enforces:
//
//      • safety (estimators cannot mutate the configuration)
//      • modularity (estimators depend only on the view interface)
//      • future-proofing (internal representation can change)
//
//  ReplicaView and ConfigurationView can be expanded as needed to expose
//  segments, occupations, or other derived quantities.
// ============================================================================

#pragma once
#include <vector>

class ReplicaView
{
public:
    // placeholder; extend as needed
    const std::vector<int> &sites() const { return sites_; }

private:
    std::vector<int> sites_;
};

class ConfigurationView
{
public:
    const std::vector<ReplicaView> &replicas() const { return replicas_; }
    const ReplicaView &replica(int r) const { return replicas_[r]; }

private:
    std::vector<ReplicaView> replicas_;
};
