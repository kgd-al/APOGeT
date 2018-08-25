#ifndef _CROSSOVER_CONFIG_HPP_
#define _CROSSOVER_CONFIG_HPP_

#include "kgd/settings/configfile.h"
#include "kgd/settings/mutationbounds.hpp"

DEFINE_PRETTY_ENUMERATION(
  CrossoverDataMutations,
    DISTANCE, INBREED, OUTBREED
)

#define CFILE CrossoverConfig
struct CFILE : public ConfigFile<CFILE> {
  template <typename T> using B = MutationSettings::B<T>;
  template <typename E> using MutationRates = MutationSettings::MutationRates<E>;

  DECLARE_PARAMETER(float, mutateChild)

  DECLARE_PARAMETER(B<float>, optimalDistance)
  DECLARE_PARAMETER(B<float>, inbreedTolerance)
  DECLARE_PARAMETER(B<float>, outbreedTolerance)

  DECLARE_PARAMETER(MutationRates<CrossoverDataMutations>, cdMutations)
};
#undef CFILE

#endif // _CROSSOVER_CONFIG_HPP_
