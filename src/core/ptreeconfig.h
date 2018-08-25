#ifndef _P_TREE_CONFIG_H_
#define _P_TREE_CONFIG_H_

#include "kgd/settings/configfile.h"

#define CFILE PTreeConfig
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

  DECLARE_PARAMETER(int, DEBUG)
};
#undef CFILE

#endif // _P_TREE_CONFIG_H_
