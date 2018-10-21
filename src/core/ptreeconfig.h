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

  /// Number of genome stored as enveloppe points
  DECLARE_PARAMETER(uint, enveloppeSize)

  /// Whether to put extra effort in creating new species or plain singletons
  DECLARE_PARAMETER(bool, simpleNewSpecies)

  /// Whether to bother drawing extinct paths
  DECLARE_PARAMETER(bool, winningPathOnly)

  /// (Debug) selector for the enveloppe criteria
  DECLARE_DEBUG_PARAMETER(int, ENV_CRIT, 0)

  /// Debug constant (compile-constant equal to false in all but debug builds)
  DECLARE_DEBUG_PARAMETER(int, DEBUG, 0)
};
#undef CFILE

} // end of namespace config

#endif // _P_TREE_CONFIG_H_
