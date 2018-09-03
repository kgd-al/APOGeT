#ifndef _P_TREE_CONFIG_H_
#define _P_TREE_CONFIG_H_

#include "kgd/settings/configfile.h"

namespace config {

#define CFILE PTree
struct CFILE : public ConfigFile<CFILE> {

  /// Threshold for being considered a "viable mate"
  DECLARE_PARAMETER(double, compatibilityThreshold)

  /// Number of votes required to being considered as belonging
  DECLARE_PARAMETER(double, similarityThreshold)

  /// Number of enveloppe points to out-perform to be considered better
  DECLARE_PARAMETER(double, outperformanceThreshold)

  /// Number of genome stored as enveloppe points
  DECLARE_PARAMETER(uint, enveloppeSize)

  DECLARE_PARAMETER(bool, simpleNewSpecies)
  DECLARE_PARAMETER(bool, ignoreHybrids)
  DECLARE_PARAMETER(bool, winningPathOnly)

  DECLARE_DEBUG_PARAMETER(int, DEBUG, 0)
};
#undef CFILE

} // end of namespace config

#endif // _P_TREE_CONFIG_H_
