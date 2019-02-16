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

  /// Average compatibility required to being considered as member of species
  DECLARE_PARAMETER(double, avgCompatibilityThreshold)

  /// Number of genome stored as enveloppe points
  DECLARE_PARAMETER(uint, enveloppeSize)

  /// Whether to put extra effort in creating new species or plain singletons
  DECLARE_PARAMETER(bool, simpleNewSpecies)

  /// Whether to bother drawing extinct paths
  DECLARE_PARAMETER(bool, winningPathOnly)

  /// (Debug) selector for the species matching score computing type
  DECLARE_DEBUG_PARAMETER(bool, DEBUG_FULL_CONTINUOUS, true)

  /// (Debug) selector for the enveloppe criteria
  DECLARE_DEBUG_PARAMETER(int, DEBUG_ENV_CRIT, 0)

  /// How much debug information should be printed out
  DECLARE_DEBUG_PARAMETER(uint, DEBUG_LEVEL, 0)

  /// Should debug info about the ptree be printed out ?
  DECLARE_DEBUG_PARAMETER(bool, DEBUG_PTREE, 0)

  /// Should debug info about the enveloppes be printed out ?
  DECLARE_DEBUG_PARAMETER(bool, DEBUG_ENVELOPPE, 0)

  /// Should debug info about the contributors be printed out ?
  DECLARE_DEBUG_PARAMETER(bool, DEBUG_CONTRIBUTORS, 0)

  /// Should debug info about the idToSpecies map be printed out ?
  DECLARE_DEBUG_PARAMETER(bool, DEBUG_ID2SPECIES, 0)
};
#undef CFILE

} // end of namespace config

#endif // _P_TREE_CONFIG_H_
