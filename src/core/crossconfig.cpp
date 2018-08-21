#include "crossconfig.h"

#define CFILE CrossoverConfig
template <typename T> using B = typename CFILE::B<T>;

DEFINE_PARAMETER(float, mutateChild, .5)

DEFINE_PARAMETER(B<float>, optimalDistance, 0.f, true, 4.f, 100.f)
DEFINE_PARAMETER(B<float>, inbreedTolerance, 0.f, 2.f, 2.f, 10.f)
DEFINE_PARAMETER(B<float>, outbreedTolerance, 0.f, 2.f, 2.f, 10.f)

template <typename ENUM>
using MutationRates = CFILE::MutationRates<ENUM>;

using CDM = CrossoverDataMutations;
DEFINE_MAP_PARAMETER(MutationRates<CrossoverDataMutations>, cdMutations, {
  { CDM::DISTANCE, 1.f },
  {  CDM::INBREED, 1.f },
  { CDM::OUTBREED, 1.f },
})

#undef CFILE
