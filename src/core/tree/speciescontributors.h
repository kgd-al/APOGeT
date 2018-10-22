#ifndef SPECIESCONTRIBUTORS_H
#define SPECIESCONTRIBUTORS_H

/*!
 * \file speciescontributors.h
 *
 * Contains the definition for species hybridism watch mechanism
 */

#include "treetypes.h"

namespace phylogeny {

/// Contributor field for a species node
class NodeContributor {
  SID _speciesID;  ///< Reference to the contributor
  uint _count; ///< Number of contributions

public:
  /// No-argument contructor
  NodeContributor(void) : NodeContributor(SID::INVALID, -1) {}

  /// Constructor
  NodeContributor(SID sid, uint initialCount)
    : _speciesID(sid), _count(initialCount) {}

  /// Accessor to the contributor reference
  SID speciesID (void) const {
    return _speciesID;
  }

  /// Increment the number of contributions
  NodeContributor& operator+=(uint k) {
    _count += k;
    return *this;
  }

  /// Compare according to the respective number of contributions
  friend bool operator< (const NodeContributor &lhs, const NodeContributor &rhs) {
    // bigger contributions go first in the array
    return !(lhs._count < rhs._count);
  }

  /// Serialize Contributor \p c into a json
  friend void to_json (json &j, const NodeContributor &c) {
    j = {c._speciesID, c._count};
  }

  /// Deserialize Contributor \p c from json \p j
  friend void from_json (const json &j, NodeContributor &c) {
    uint i=0;
    c._speciesID = j[i++];
    c._count = j[i++];
  }
};

/// Sorted collection of contributors for a species node
class Contributors {
  /// The buffer containing the individual contributions
  std::vector<NodeContributor> vec;

  /// The associated node identificator
  SID nodeID;

public:
  /// Alias for the data structure containing the contributing SIDs
  using Contribution = std::multiset<SID>;

  /// No-argument constructor. Leaves the class in an invalid state
  Contributors (void) : Contributors(SID::INVALID) {}

  /// Constructor. Registers the node whose contributor collection it manages
  Contributors (SID id) : nodeID(id) {}

  /// \return The id of the monitored node
  SID getNodeID (void) const {
    return nodeID;
  }

  /// Register new contributions, updates internal data and returns the new
  /// main contributor
  SID update (Contribution sids);

  /// Update the main contributor cached variable so that if nodeID = X
  /// [] -> INVALID
  /// [X] -> INVALID
  /// [X,Y,...] -> Y
  /// [Y,...,X,...] -> Y
  SID currentMain (void);

  /// Serialize Contributors \p c into a json
  friend void to_json (json &j, const Contributors &c);

  /// Deserialize Contributors \p c from json \p j
  friend void from_json (const json &j, Contributors &c);
};

} // end of namespace phylogeny

#endif // SPECIESCONTRIBUTORS_H
