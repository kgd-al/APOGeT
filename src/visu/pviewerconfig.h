#ifndef _P_VIEWER_CONFIG_H_
#define _P_VIEWER_CONFIG_H_

#include "kgd/settings/configfile.h"

/*!
 * \file pviewerconfig.h
 *
 * Definitions for displaying the phylogenic process
 */

namespace config {

struct PTree;

#define CFILE PViewer

/// Config file for the phylogenic viewer
struct CFILE : public ConfigFile<CFILE> {
  /// Whether to draw species nodes id
  DECLARE_PARAMETER(bool, showNodeNames)

  /// Do not draw species having lived less than that amount
  DECLARE_PARAMETER(uint, minNodeSurvival)

  /// Do not draw species having less than this ratio of fullness
  DECLARE_PARAMETER(uint, minNodeEnveloppe)

  /// Whether to bother drawing extinct paths
  DECLARE_PARAMETER(bool, survivorNodesOnly)

  /// How much detail to print for species aggregation
  DECLARE_PARAMETER(uint, speciesDetailVerbosity)

  /// Core config
  DECLARE_SUBCONFIG(PTree, coreConfig)
};
#undef CFILE

} // end of namespace config

#endif // _P_TREE_CONFIG_H_
