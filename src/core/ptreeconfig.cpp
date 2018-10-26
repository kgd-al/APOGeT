#include "ptreeconfig.h"

namespace config {

#define CFILE PTree

DEFINE_PARAMETER(double, compatibilityThreshold, .1)
DEFINE_PARAMETER(double, similarityThreshold, .5)
DEFINE_PARAMETER(double, avgCompatibilityThreshold, .25)
DEFINE_PARAMETER(uint, enveloppeSize, 5)
DEFINE_PARAMETER(bool, simpleNewSpecies, true)
DEFINE_PARAMETER(bool, winningPathOnly, false)

DEFINE_DEBUG_PARAMETER(bool, FULL_CONTINUOUS, true)
DEFINE_DEBUG_PARAMETER(int, ENV_CRIT, 1)
DEFINE_DEBUG_PARAMETER(uint, DEBUG_LEVEL, 0)
DEFINE_DEBUG_PARAMETER(bool, DEBUG_PTREE, 1)
DEFINE_DEBUG_PARAMETER(bool, DEBUG_ENVELOPPE, 0)
DEFINE_DEBUG_PARAMETER(bool, DEBUG_CONTRIBUTORS, 1)

#undef CFILE

} // end of namespace config

