#include "ptreeconfig.h"

namespace config {

#define CFILE PTree

DEFINE_PARAMETER(double, compatibilityThreshold, .1)
DEFINE_PARAMETER(double, similarityThreshold, .5)
DEFINE_PARAMETER(uint, enveloppeSize, 5)
DEFINE_PARAMETER(bool, simpleNewSpecies, true)
DEFINE_PARAMETER(bool, ignoreHybrids, true)
DEFINE_PARAMETER(bool, winningPathOnly, false)

DEFINE_DEBUG_PARAMETER(int, ENV_CRIT, 0)
DEFINE_DEBUG_PARAMETER(int, DEBUG, 0)

#undef CFILE

} // end of namespace config

