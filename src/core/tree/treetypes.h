#ifndef _PTREE_TYPES_HPP_
#define _PTREE_TYPES_HPP_

/*!
 * \file treetypes.h
 *
 * Contains the definition of the generic types used across the whole phylogenic
 * algorithms
 */

#include <type_traits>
#include <set>

#include "kgd/external/json.hpp"

#include "../crossover.h"

namespace phylogeny {

/// Helper alias to the json type used for (de)serialization
using json = nlohmann::json;

/// Helper alias for the genome identificator
using GID = genotype::BOCData::GID;

/// Alias for the species identificator
enum class SID : uint {
  INVALID = uint(-1)  ///< Value indicating an unspecified species
};

/// Auto convert outstream operator
std::ostream& operator<< (std::ostream &os, SID sid);

/// Collections of still-alive species identificators
using LivingSet = std::set<SID>;


namespace _details {

/// Distance & compatibilities cache
struct DCCache {
  /// Cache collection of distances
  std::vector<float> distances;

  /// Cache collection of compatibilities
  std::vector<float> compatibilities;

  /// Remove all contents
  void clear (void) {
    distances.clear(),  compatibilities.clear();
  }

  /// Prepare exactly \p n units of storage space
  void reserve (uint n) {
    distances.reserve(n), compatibilities.reserve(n);
  }

  /// Append values
  void push_back (float d, float c) {
    distances.push_back(d), compatibilities.push_back(c);
  }

  /// \returns the size of the cache
  size_t size (void) const {
    assert(distances.size() == compatibilities.size());
    return distances.size();
  }
};


/// Helper structure for ensuring that the pair values are ordered
///
/// \tparam T stored type
template <typename T>
struct ordered_pair {
  T first;  ///< First value (lower or equal to second)
  T second; ///< Second value (greater or equal to first)

  /// Create an ordered pair from un-ordered pair of values
  ordered_pair(T first, T second)
    : first(std::min(first, second)), second(std::max(first, second)) {}

  /// Compare two ordered_pairs in lexicographic order
  friend bool operator< (const ordered_pair &lhs, const ordered_pair &rhs) {
    if (lhs.first != rhs.first) return lhs.first < rhs.first;
    return lhs.second < rhs.second;
  }
};

/// Helper alias to a map whose keys are ordered so that
/// \f$\forall i,j \in M i<j\f$
using DistanceMap = std::map<ordered_pair<uint>, float>;

/// Description of the contribution of a genome to a species enveloppe
struct EnveloppeContribution {
  bool better;  ///< Should an enveloppe point be replaced
  uint than;    ///< Index of the enveloppe point to replace
  float value;  ///< Confidence of the replacement pertinence
};

/// Helper alias to a genome identificator
using GID = genotype::BOCData::GID;

/// Computes whether or not the considered species would be better described by
/// replacing a point from the current enveloppe (with distance map \p edist)
/// by an incoming genome (with distances \p gdist)
EnveloppeContribution computeContribution(const DistanceMap &edist,
                                          const std::vector<float> &gdist,
                                          GID gid, const std::vector<GID> &ids);

} // end of namespace _details

} // end of namespace phylogeny

#endif // _PTREE_TYPES_HPP_
