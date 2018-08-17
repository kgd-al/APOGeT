#ifndef _P_TREE_CONFIG_H_
#define _P_TREE_CONFIG_H_

#include "settings/configfile.h"

#define CFILE PTreeConfig
struct CFILE : public ConfigFile<CFILE> {
  DECLARE_PARAMETER(double, compatibilityThreshold)
  DECLARE_PARAMETER(double, similarityThreshold)
  DECLARE_PARAMETER(uint, enveloppeSize)
  DECLARE_PARAMETER(bool, simpleNewSpecies)
  DECLARE_PARAMETER(bool, ignoreHybrids)
  DECLARE_PARAMETER(bool, winningPathOnly)
};
#undef CFILE

#endif // _P_TREE_CONFIG_H_
