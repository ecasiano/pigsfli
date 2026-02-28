// ============================================================================
//  estimator.hpp
//  -------------
//  Defines the abstract Estimator interface.
//
//  Estimators read the configuration (via ConfigurationView) and accumulate
//  measurements. They NEVER modify the worldline configuration.
//
//  Examples:
//
//      • density profiles
//      • Green’s functions
//      • structure factors
//      • energy estimators
//      • entanglement estimators (via SWAP)
//
//  Estimators are plug-ins: you can add new ones without touching the engine.
// ============================================================================

#pragma once
#include <ostream>

class ConfigurationView; // fwd

class Estimator
{
public:
    virtual ~Estimator() = default;

    virtual void accumulate(const ConfigurationView &cfg, double weight) = 0;
    virtual void normalize() = 0;
    virtual void write_output(std::ostream &out) const = 0;
};
