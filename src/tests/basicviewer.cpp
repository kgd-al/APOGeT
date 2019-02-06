#include "../core/crossover.h"
#include "../visu/standaloneviewer.hpp"

/*!
 * \file basicviewer.cpp
 *
 * Contains the &nbsp; \copydoc main
 */

/// Decoy genome with no internal structure
struct Genome {
  /// Should convert the genome to json but, in fact, does nothing
  friend void to_json (nlohmann::json&, const Genome&) {}

  /// Should convert the json into the genome but, in fact, does nothing
  friend void from_json (const nlohmann::json&, Genome&) {}
};

/// main for a straightforward executable that can display any PTree
/// by throwing genetic information away
int main(int argc, char *argv[]) {
  return run<Genome, config::PTree>(argc, argv);
}
