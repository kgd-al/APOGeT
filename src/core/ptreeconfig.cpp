#include "ptreeconfig.h"

namespace config {

#define CFILE PTree

DEFINE_PARAMETER(double, compatibilityThreshold, .1)
DEFINE_PARAMETER(double, similarityThreshold, .5)
DEFINE_PARAMETER(double, avgCompatibilityThreshold, .25)
DEFINE_PARAMETER(uint, enveloppeSize, 5)
DEFINE_PARAMETER(bool, simpleNewSpecies, true)
DEFINE_PARAMETER(bool, winningPathOnly, false)

DEFINE_PARAMETER(uint, stillbornTrimmingPeriod, 100)
DEFINE_PARAMETER(float, stillbornTrimmingThreshold, .25)
DEFINE_PARAMETER(float, stillbornTrimmingDelay, 4)
DEFINE_PARAMETER(uint, stillbornTrimmingMinDelay, 200)

DEFINE_DEBUG_PARAMETER(bool, DEBUG_FULL_CONTINUOUS, true)
DEFINE_DEBUG_PARAMETER(int, DEBUG_ENV_CRIT, 1)

DEFINE_DEBUG_PARAMETER(uint, DEBUG_LEVEL, 0)
DEFINE_DEBUG_PARAMETER(bool, DEBUG_PTREE, 0)
DEFINE_DEBUG_PARAMETER(bool, DEBUG_ENVELOPPE, 0)
DEFINE_DEBUG_PARAMETER(bool, DEBUG_CONTRIBUTORS, 0)
DEFINE_DEBUG_PARAMETER(bool, DEBUG_ID2SPECIES, 0)
DEFINE_DEBUG_PARAMETER(bool, DEBUG_STILLBORNS, 0)

#undef CFILE

} // end of namespace config

