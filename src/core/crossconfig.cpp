#include "crossconfig.h"
#include "crossover.hpp"

namespace config {

#define CFILE Crossover
using B = typename CFILE::B;
using D = genotype::BOCData;

DEFINE_PARAMETER(float, mutateChild, .5)

DEFINE_PARAMETER(B, optimalDistance, &D::optimalDistance, 0.f, true, 4.f, 100.f)
DEFINE_PARAMETER(B, inbreedTolerance, &D::inbreedTolerance, 0.f, 2.f, 2.f, 10.f)
DEFINE_PARAMETER(B, outbreedTolerance, &D::outbreedTolerance, 0.f, 2.f, 2.f, 10.f)

template <typename ENUM>
using MutationRates = CFILE::MutationRates<ENUM>;

using CDM = CrossoverDataMutations;
DEFINE_MAP_PARAMETER(MutationRates<CrossoverDataMutations>, cdMutations, {
  { CDM::DISTANCE, 1.f },
  {  CDM::INBREED, 1.f },
  { CDM::OUTBREED, 1.f },
})

} // end of namespace config

#undef CFILE
