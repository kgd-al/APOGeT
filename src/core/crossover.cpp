#include "crossover.h"

#define GENOME BOCData

DEFINE_GENOME_FIELD_WITH_BOUNDS(float, optimalDistance, 0.f, true, 4.f, 100.f)
DEFINE_GENOME_FIELD_WITH_BOUNDS(float, inbreedTolerance, 0.f, 2.f, 2.f, 10.f)
DEFINE_GENOME_FIELD_WITH_BOUNDS(float, outbreedTolerance, 0.f, 2.f, 2.f, 10.f)

DEFINE_GENOME_MUTATION_RATES({
  MUTATION_RATE(optimalDistance, 1.f),
  MUTATION_RATE(inbreedTolerance, 1.f),
  MUTATION_RATE(outbreedTolerance, 1.f),
})

#undef GENOME

#define CFILE config::SAGConfigFile<genotype::BOCData>

DEFINE_PARAMETER(float, mutateChild, .5)

#undef CFILE
