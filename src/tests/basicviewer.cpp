#include "../core/crossover.hpp"
#include "../visu/standaloneviewer.hpp"

/*!
 * \file basicviewer.cpp
 *
 * Contains the &nbsp; \copydoc main
 */

/// Decoy genome with no internal structure
struct Genome {
  /// Decoy alignment for decoy genome
  struct Alignment {};

  /// Should convert the genome to json but, in fact, does nothing
  friend void to_json (nlohmann::json&, const Genome&) {}

  /// Should convert the json into the genome but, in fact, does nothing
  friend void from_json (const nlohmann::json&, Genome&) {}
};

namespace config {
#define CFILE Basic
/// Bare-bones configuration file. Placeholder for coreconfig
struct Basic : public ConfigFile<CFILE> {
  /// Reference to the parameters for the crossover algorithms
  DECLARE_SUBCONFIG(Crossover, crossoverConfig)

  /// Reference to the parameters for the phylogenic algorithms
  DECLARE_SUBCONFIG(PTree, phenotypicTreeConfig)
};
DEFINE_SUBCONFIG(Crossover, crossoverConfig)
DEFINE_SUBCONFIG(PTree, phenotypicTreeConfig)
#undef CFILE
} // end of namespace config

/// main for a straightforward executable that can display any PTree
/// by throwing genetic information away
int main(int argc, char *argv[]) {
  return run<Genome, config::Basic>(argc, argv);
}
