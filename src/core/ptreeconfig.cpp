#include "ptreeconfig.h"

#define CFILE PTreeConfig

DEFINE_PARAMETER(double, compatibilityThreshold, .1)
DEFINE_PARAMETER(double, similarityThreshold, .5)
DEFINE_PARAMETER(uint, enveloppeSize, 5)
DEFINE_PARAMETER(bool, simpleNewSpecies, true)
DEFINE_PARAMETER(bool, ignoreHybrids, true)
DEFINE_PARAMETER(bool, winningPathOnly, false)

#undef CFILE
