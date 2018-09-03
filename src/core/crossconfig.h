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
  using B = MutationSettings::Bounds<float, genotype::BOCData>;
  template <typename E> using MutationRates = MutationSettings::MutationRates<E>;

  /// Probability of mutating a child after crossover
  DECLARE_PARAMETER(float, mutateChild)

  DECLARE_PARAMETER(B, optimalDistance)
  DECLARE_PARAMETER(B, inbreedTolerance)
  DECLARE_PARAMETER(B, outbreedTolerance)

  DECLARE_PARAMETER(MutationRates<CrossoverDataMutations>, cdMutations)
};
#undef CFILE

} // end of namespace config

#endif // _CROSSOVER_CONFIG_HPP_
