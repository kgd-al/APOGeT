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
class Contributor {
  SID _speciesID;  ///< Reference to the contributor
  uint _count;  ///< Number of contributions
  bool _elligible;  ///< Valid candidate for being the major contributor ?

public:
  /// No-argument contructor
  Contributor(void) : Contributor(SID::INVALID, -1, false) {}

  /// Constructor
  Contributor(SID sid, uint initialCount, bool elligible)
    : _speciesID(sid), _count(initialCount), _elligible(elligible) {}

  /// Accessor to the contributor reference
  SID speciesID (void) const {
    return _speciesID;
  }

  /// Accessor to the contribution count
  uint count (void) const {
    return _count;
  }

  /// \returns the validity of this contributor
  /// \see _elligible
  bool elligible (void) const { return _elligible;  }

  /// Updates the validity of this contributor
  /// \see _elligible
  void setElligible (bool e) {  _elligible = e;  }

  /// Increment the number of contributions
  Contributor& operator+=(uint k) {
    _count += k;
    return *this;
  }

  /// Contributor comparison functor
  ///
  /// Compare according to the respective number of contributions
  struct CMP {

    /// Contributor comparison function
    ///
    /// \copydetails CMP
    bool operator() (const Contributor &lhs, const Contributor &rhs) {
      // bigger contributions go first in the array
      return lhs._count > rhs._count;
    }
  };

  /// Serialize Contributor \p c into a json
  friend void to_json (json &j, const Contributor &c) {
    j = {c._speciesID, c._count};
  }

  /// Deserialize Contributor \p c from json \p j
  friend void from_json (const json &j, Contributor &c) {
    uint i=0;
    c._speciesID = j[i++];
    c._count = j[i++];
  }
};

/// Describes a contribution update
struct Contribution {
  SID species;  ///< The contributing species
  uint count;   ///< The contribution amount

  /// Constructor
  Contribution (SID s, uint c) : species(s), count(c) {}

  /// Allow comparing a contribution's species id directly to another sid
  friend bool operator== (const Contribution &c, SID sid) {
    return c.species == sid;
  }
};

/// Sorted collection of contributors for a species node.
///
/// Use both for maintaining phylogenic data and determining the major
/// contributor.
///
/// The major contributor is defined as follow:
/// Given A, an instance of contributors, return the first element in this (sorted)
/// container that is flagged as elligible
///
/// A species B is elligible as another species A's parent iff:
///   - B is not a node in A's subtree (including itself)
///   - B is not younger than A (in terms of first appearance)
class Contributors {
  /// The buffer containing the individual contributions
  std::vector<Contributor> vec;

  /// The associated node identificator
  SID nodeID;

public:
  /// Alias for the data structure containing the contributing SIDs
  using Contributions = std::vector<Contribution>;

  /// Alias for the function type used to check is a species is elligible as a
  /// major contributor
  using ValidityEvaluator = std::function<bool(SID,SID)>;

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
  SID update (Contributions ctbs, const ValidityEvaluator &elligible);

  /// \return the id of the node's main contributor or SID::INVALID if none is found
  SID currentMain (void);

  /// Updates, for each contributions, whether it is coming from a valid
  /// candidate to being a major contributor or not
  ///
  /// \return the updated parent
  SID updateElligibilities (const ValidityEvaluator &elligible);

  /// Allow const iteration of the underlying container
  const auto begin (void) const {
    return vec.begin();
  }

  /// \copydoc begin
  const auto end (void) const {
    return vec.end();
  }

  /// Serialize Contributors \p c into a json
  friend void to_json (json &j, const Contributors &c);

  /// Deserialize Contributors \p c from json \p j
  friend void from_json (const json &j, Contributors &c);

  /// Stream Contributors \p c to \p os
  friend std::ostream& operator<< (std::ostream &os, const Contributors &c);

  /// \returns Whether or not \p rhs is a valid candidate for being the major
  /// contributor of species \p lhs
  template <typename T>
  static bool elligibile (SID lhs, SID rhs, const T &nodes) {
    auto n = nodes.at(lhs).get(),
         p = nodes.at(rhs).get();

    // Do not allow younger species to serve as parent (would be quite ugly and
    // is probably wrong anyway)
    if (n->data.firstAppearance <= p->data.firstAppearance)
      return false;

    // Assert that candidate is not in n's subtree
    while (p && p != n)
      p = p->parent();
    return p != n;
  }
};

} // end of namespace phylogeny

#endif // SPECIESCONTRIBUTORS_H
