#include "pviewerconfig.h"
#include "../core/ptreeconfig.h"

namespace config {

#define CFILE PViewer

DEFINE_PARAMETER(bool, showNodeNames, true)
DEFINE_PARAMETER(uint, minNodeSurvival, 0)
DEFINE_PARAMETER(uint, minNodeEnveloppe, 0)
DEFINE_PARAMETER(bool, survivorNodesOnly, false)
DEFINE_PARAMETER(uint, speciesDetailVerbosity, 1)

DEFINE_SUBCONFIG(PTree, coreConfig)

#undef CFILE

} // end of namespace config

