#ifndef SPECIESDATA_HPP
#define SPECIESDATA_HPP

/*!
 * \file speciesdata.hpp
 *
 * Contains the definition for the species data storage structure
 */

#include "treetypes.h"

namespace phylogeny {

/// Contains various data about a specific species
struct SpeciesData {

  /// Time at which this species was first observed
  uint firstAppearance;

  /// Time at which this species was last observed
  uint lastAppearance;

  /// Number of individuals of the species observed
  uint count;

  /// Serialize to json
  friend void to_json (json &j, const SpeciesData &d) {
    j = {d.firstAppearance, d.lastAppearance, d.count};
  }

  /// Deserialize from json
  friend void from_json (const json &j, SpeciesData &d) {
    uint i=0;
    d.firstAppearance = j[i++];
    d.lastAppearance = j[i++];
    d.count = j[i++];
  }
};

} // end of namespace phylogeny

#endif // SPECIESDATA_HPP
