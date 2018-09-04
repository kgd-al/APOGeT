#ifndef _P_TREE_CONFIG_H_
#define _P_TREE_CONFIG_H_

#include "kgd/settings/configfile.h"

/*!
 * \file ptreeconfig.h
 *
 * Definitions for controlling the phylogenic process
 */

namespace config {

#define CFILE PTree

/// Config file for the phylogenic algorithms
struct CFILE : public ConfigFile<CFILE> {

  /// Threshold for being considered a "viable mate"
  DECLARE_PARAMETER(double, compatibilityThreshold)

  /// Number of votes required to being considered as belonging
  DECLARE_PARAMETER(double, similarityThreshold)

  /// Number of enveloppe points to out-perform to be considered better
  DECLARE_PARAMETER(double, outperformanceThreshold)

  /// Number of genome stored as enveloppe points
  DECLARE_PARAMETER(uint, enveloppeSize)

  /// Whether to put extra effort in creating new species or plain singletons
  DECLARE_PARAMETER(bool, simpleNewSpecies)

  /// What to do with hybrid individuals
  DECLARE_PARAMETER(bool, ignoreHybrids)

  /// Whether to bother drawing extinct paths
  DECLARE_PARAMETER(bool, winningPathOnly)

  /// Debug constant (compile-constant equal to false in all but debug builds)
  DECLARE_DEBUG_PARAMETER(int, DEBUG, 0)
};
#undef CFILE

} // end of namespace config

#endif // _P_TREE_CONFIG_H_
