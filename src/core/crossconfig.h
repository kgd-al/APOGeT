#ifndef _CROSSOVER_CONFIG_HPP_
#define _CROSSOVER_CONFIG_HPP_

#include "kgd/settings/configfile.h"
#include "kgd/settings/mutationbounds.hpp"

/*!
 * \file crossconfig.h
 *
 * Definitions for controlling the crossover process
 */

/// \enum CrossoverDataMutations
/// Fields that should be mutated in the crossover data structure
DEFINE_NAMESPACE_PRETTY_ENUMERATION(
  config, CrossoverDataMutations,
    DISTANCE, INBREED, OUTBREED
)

namespace genotype {
struct BOCData;
}

namespace config {

#define CFILE Crossover

/// Config file for the crossover algorithm
struct CFILE : public ConfigFile<CFILE> {

  /// Handy alias for the bounds type
  using B = MutationSettings::Bounds<float, genotype::BOCData>;

  /// Handy alias for the mutation rates
  using MutationRates = MutationSettings::MutationRates<CrossoverDataMutations>;

  /// Probability of mutating a child after crossover
  DECLARE_PARAMETER(float, mutateChild)

  /// Mutation bounds for the optimal genetic distance
  /// \see genotype::BOCData::optimalDistance
  DECLARE_PARAMETER(B, optimalDistance)

  /// Mutation bounds for the inbreed tolerance
  /// \see genotye::BOCData::inbreedTolerance
  DECLARE_PARAMETER(B, inbreedTolerance)

  /// Mutation bounds for the outbreed tolerance
  /// \see genotype::BOCData::outbreedTolerance
  DECLARE_PARAMETER(B, outbreedTolerance)

  /// Mutation rates for the BOCData fields
  DECLARE_PARAMETER(MutationRates, cdMutations)
};
#undef CFILE

} // end of namespace config

#endif // _CROSSOVER_CONFIG_HPP_
