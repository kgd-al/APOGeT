#ifndef KGD_PTREE_TYPES_HPP
#define KGD_PTREE_TYPES_HPP

/*!
 * \file treetypes.h
 *
 * Contains the definition of the generic types used across the whole phylogenic
 * algorithms
 */

#include <type_traits>
#include <set>

#include "kgd/external/json.hpp"
#include "kgd/utils/utils.h"
#include "kgd/utils/assertequal.hpp"

namespace phylogeny {

/// Helper alias to the json type used for (de)serialization
using json = nlohmann::json;

/// Defined type for GenomeID
enum class GID : uint {
  INVALID = uint(-1)  ///< Value indicating an unspecified genome
};

/// Auto-convert outstream operator
std::ostream& operator<< (std::ostream &os, GID gid);

/// Manager for generating the succession of genetic identificators
class GIDManager {
  /// Helper alias to the underlying type of the genetic identificators
  using GID_ut = std::underlying_type<GID>::type;

  /// Next genome identificator
  GID_ut next;

public:
  /// Sets the first identificator to 0
  GIDManager (void) {
    setNext(GID(0));
  }

  /// Sets next genome id
  void setNext (GID next) {
    this->next = GID_ut(next) + 1;
  }

  /// Generate next id value
  GID operator() (void) {
    if (next == std::numeric_limits<GID_ut>::max())
      utils::Thrower<std::out_of_range>(
        "Exhausted all possible identifers for the underlying type ",
        utils::className<GID_ut>());
    return GID(next++);
  }

  /// \returns value for next id
  /// \warning This does not modify the current value. Use operator() instead
  explicit operator GID(void) const {
    return GID(next);
  }

  /// Used to assert correct cloning (designed for the EDEnS algorithm)
  friend void assertEqual (const GIDManager &lhs, const GIDManager &rhs,
                           bool deepcopy) {
    using utils::assertEqual;
    assertEqual(lhs.next, rhs.next, deepcopy);
  }
};

/// Alias for the species identificator
enum class SID : uint {
  INVALID = uint(-1)  ///< Value indicating an unspecified species
};

/// Auto convert outstream operator
std::ostream& operator<< (std::ostream &os, SID sid);

/// Collections of still-alive species identificators
using LivingSet = std::set<SID>;

/// Holds the identificators for a given individual
struct PID {
  GID gid; ///< Identificator in the genomic population
  SID sid;  ///< Identificator of the associated species

  /// Creates an invalid identificator
  explicit PID (void) : PID(GID::INVALID) {}

  /// Creates an identificator from a genomic id
  explicit PID (GID gid) : gid(gid), sid(SID::INVALID) {}

  /// \returns Whether this ID belong to an existing genome
  bool valid (void) const { return gid != GID::INVALID; }

  /// Compares two IDs for equality
  friend bool operator== (const PID &lhs, const PID &rhs) {
    return lhs.gid == rhs.gid && lhs.sid == rhs.sid;
  }

  /// Serializes an ID into json form
  friend void to_json (json &j, const PID &pid) {
    j["g"] = pid.gid;
    j["s"] = pid.sid;
  }

  /// Deserializes an ID from a json
  friend void from_json (const json &j, PID &pid) {
    pid.gid = j["g"];
    pid.sid = j["s"];
  }

  /// Writes an ID into the provided outstream
  friend std::ostream& operator<< (std::ostream &os, const PID &pid) {
    return os << "{G: " << pid.gid << ", S: " << pid.sid << "}";
  }
};

/// Holds the identificators of a given genome and its parents (if any)
struct Genealogy {
  PID mother;  ///< Identificators of the genome's primary parent
  PID father;  ///< Identificators of the genome's secondary parent
  PID self;  ///< Identificators of the given genome

  uint generation;  ///< Generation of the given geno;e

  /// Sets the species ID
  void setSID (SID sid) { self.sid = sid; }

  /// Updates internal data to reflect the clone status of the associated genome
  void updateAfterCloning (GIDManager &m) {
    mother = self;
    father = PID();
    self = PID(m());
    generation++;
  }

  /// Updates internal data to reflect the child status of the associated genome
  void updateAfterCrossing(const Genealogy &mother, const Genealogy &father,
                           GIDManager &m) {
    this->mother = mother.self;
    this->father = father.self;
    self = PID(m());
    generation = std::max(mother.generation, father.generation) + 1;
  }

  /// Sets up value to indicate a primordial genome
  void setAsPrimordial(GIDManager &m) {
    self = PID(m());
    mother = PID();
    father = PID();
    generation = 0;
  }

  /// Serializes genealogy into json form
  friend void to_json (json &j, const Genealogy &g) {
    j["m"] = g.mother;
    j["f"] = g.father;
    j["s"] = g.self;
  }

  /// Deserializes genealogy from a json
  friend void from_json (const json &j, Genealogy &g) {
    g.mother = j["m"];
    g.father = j["f"];
    g.self = j["s"];
  }

  /// Compare two genealogies for equality
  friend bool operator== (const Genealogy &lhs, const Genealogy &rhs) {
    return lhs.self == rhs.self
        && lhs.mother == rhs.mother
        && lhs.father == rhs.father
        && lhs.generation == rhs.generation;
  }

  /// Writes a genealogy into an outstream
  friend std::ostream& operator<< (std::ostream &os, const Genealogy &g) {
    os << "{ ";
    if (g.mother.valid()) os << "M: " << g.mother << ", ";
    if (g.father.valid()) os << "F: " << g.father << ", ";
    if (g.self.valid())   os << "S: " << g.self;
    else                  os << "PRIMORDIAL";
    return os << " }";
  }
};

/// Contains the result from an insertion into the tree
template <typename UserData>
struct InsertionResult {
  SID sid;  ///< Associated species
  UserData *udata;  ///< User data (if part of the rset)
};

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

  /// Asserts that two ordered pairs are equal
  friend void assertEqual (const ordered_pair &lhs, const ordered_pair &rhs,
                           bool deepcopy) {
    using utils::assertEqual;
    assertEqual(lhs.first, rhs.first, deepcopy);
    assertEqual(lhs.second, rhs.second, deepcopy);
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

/// Computes whether or not the considered species would be better described by
/// replacing a point from the current enveloppe (with distance map \p edist)
/// by an incoming genome (with distances \p gdist)
EnveloppeContribution computeContribution(const DistanceMap &edist,
                                          const std::vector<float> &gdist,
                                          GID gid, const std::vector<GID> &ids);

} // end of namespace _details

} // end of namespace phylogeny

#endif // KGD_PTREE_TYPES_HPP
