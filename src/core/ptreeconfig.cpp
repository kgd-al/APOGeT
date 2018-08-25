#include "ptreeconfig.h"

#define CFILE PTreeConfig

DEFINE_PARAMETER(double, compatibilityThreshold, .1)
DEFINE_PARAMETER(double, similarityThreshold, .5)
DEFINE_PARAMETER(double, outperformanceThreshold, .75)
DEFINE_PARAMETER(uint, enveloppeSize, 5)
DEFINE_PARAMETER(bool, simpleNewSpecies, true)
DEFINE_PARAMETER(bool, ignoreHybrids, true)
DEFINE_PARAMETER(bool, winningPathOnly, false)

DEFINE_PARAMETER(int, DEBUG, 0)

#undef CFILE
