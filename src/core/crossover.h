#ifndef _CROSSOVER_HPP_
#define _CROSSOVER_HPP_

#include <cmath>
#include <functional>
#include <atomic>

#include "kgd/external/json.hpp"

#include "kgd/genotype/selfawaregenome.hpp"

/*!
 * \file crossover.h
 *
 * Definition of structures and algorithms for implementing the Bail-Out Crossover
 */

namespace genotype {

/// Common crossover control data
class SELF_AWARE_GENOME(BOCData) {
  /// The compatibility function. Half of an unnormalized gaussian.
  double gaussoid (double x, double mu, double sigma) const {
    return exp(-((x-mu)*(x-mu)/(2.*sigma*sigma)));
  }

  // ========================================================================
  // == Compatibility function

  /// Genetic distance that maximises reproduction compatibility
  DECLARE_GENOME_FIELD(float, optimalDistance)

  /// Standard deviation for distances below optimal
  DECLARE_GENOME_FIELD(float, inbreedTolerance)

  /// Standard deviation for distances above optimal
  DECLARE_GENOME_FIELD(float, outbreedTolerance)

  /// Needs privileged access
  friend struct config::SAGConfigFile<BOCData>;

public:
  /// Helper alias to the source of randomness
  using Dice = SelfAwareGenome<BOCData>::Dice;

  /// The possible sexes
  enum Sex { FEMALE, MALE };

  /// Which sex the associated genome codes for
  DECLARE_GENOME_FIELD(Sex, sex)

  // ========================================================================
  // == Classification data

  /// Defined type for GenomeID
  using GID = uint;

  /// Definition of invalid ID
  static constexpr GID NO_ID = GID(-1);

  /// The id of the associated genome
  GID id;

  /// Generate next id value (thread-safe)
  static GID nextID (void) {
    static std::atomic<GID> NEXT_ID (0);
    return NEXT_ID++;
  }

  /// The set of parents
  enum Parent : uint { MOTHER = 0, FATHER = 1 };

  /// The id of the associated genome's parents, if any. NO_ID otherwise
  GID parents [2];

  /// The generation at which the associated genome appeared
  uint generation;

  // ========================================================================
  // == Member function(s)

  /// \return Whether or not the genome has registered a parent \p p
  bool hasParent (Parent p) const {
    return parents[p] != NO_ID;
  }

  /// \return The id of parent \p p
  GID parent (Parent p) const {
    return parents[p];
  }

  /// \return the optimal genetic distance
  /// \see optimalDistance
  float getOptimalDistance (void) const { return optimalDistance; }

  /// \return the inbreed tolerance
  /// \see inbreedTolerance
  float getInbreedTolerance (void) const { return inbreedTolerance; }

  /// \return the outbreed tolerance
  /// \see outbreedTolerance
  float getOutbreedTolerance (void) const { return outbreedTolerance; }

  /// Evaluates the compatibility between this genome and another at \p distance
  /// According to the formula:
  ///
  /// \f{equation}{ e^{\frac{-(d-\mu)^2}{2\sigma^2}} \f}
  ///
  /// with
  ///   - d, the distance
  ///   - \f$ \mu \f$, the optimal distance
  ///   - \f$ \sigma = \f$
  ///     - \f$ \sigma_i,\f$ if \f$ d < \mu \f$ (inbreed tolerance)
  ///     - \f$ \sigma_o, \f$ otherwise (outbreed tolerance)
  ///
  /// \see optimalDistance, inbreedTolerance, outbreedTolerance, gaussoid
  double operator() (double distance) const {
    return gaussoid(distance, optimalDistance,
                    distance < optimalDistance ? inbreedTolerance : outbreedTolerance);
  }

  // ========================================================================
  // == Specific genetic operators

  /// Updates internal data to reflect the clone status of the associated genome
  void updateCloneLineage (GID parent, Dice &dice) {
    sex = dice.toss(Sex::FEMALE, Sex::MALE);
    id = nextID();
    parents[MOTHER] = parent;
    parents[FATHER] = -1;
    generation++;
  }

  /// Low-level crossing of manually managed fields (id, parents, generation)
  void crossExtension(const BOCData &lhs, const BOCData &rhs, Dice&) override {
    id = nextID();
    parents[MOTHER] = lhs.id;
    parents[FATHER] = rhs.id;
    generation = std::max(lhs.generation, rhs.generation) + 1;
  }

  /// Sets up manually managed field (id, parents, generation)
  void randomExtension(Dice&) override {
    id = nextID();
    parents[MOTHER] = -1;
    parents[FATHER] = -1;
    generation = 0;
  }

  // ========================================================================
  // == Serialization

  /// Overridden to include manual fields in the comparison
  void equalExtension (const BOCData &that, bool &eq) const override {
    eq &= id == that.id
       && parents[MOTHER] == that.parents[MOTHER]
       && parents[FATHER] == that.parents[FATHER]
       && generation == that.generation;
  }

  /// Overridden to include manual fields in the serialization
  void to_jsonExtension (nlohmann::json &j) const override {
    j["id"] = id;
    j["p0"] = parents[MOTHER];
    j["p1"] = parents[FATHER];
    j["G"] = generation;
  }

  /// Overridden to include manual fields in the deserialization
  void from_jsonExtension (nlohmann::json &j) override {
    id = j["id"];               j.erase("id");
    parents[MOTHER] = j["p0"];  j.erase("p0");
    parents[FATHER] = j["p1"];  j.erase("p1");
    generation = j["G"];        j.erase("G");
  }
};

/// Pretty prints a sex enumeration value
std::ostream& operator<< (std::ostream &os, BOCData::Sex s);

/// Converts a pretty printed sex enumeration value to its integer representation
std::istream& operator>> (std::istream &is, BOCData::Sex &s);

} // end of namespace genotype

namespace config {

/// Config file for the crossover algorithms
template <> struct SAG_CONFIG_FILE(BOCData) {
  /// Helper alias to bounds object for floating point fields
  using Bf = Bounds<float>;

  /// Probability of mutating a child after crossover
  DECLARE_PARAMETER(float, mutateChild)

  /// Mutation bounds for the optimal genetic distance
  /// \see genotype::BOCData::optimalDistance
  DECLARE_PARAMETER(Bf, optimalDistanceBounds)

  /// Mutation bounds for the inbreed tolerance
  /// \see genotye::BOCData::inbreedTolerance
  DECLARE_PARAMETER(Bf, inbreedToleranceBounds)

  /// Mutation bounds for the outbreed tolerance
  /// \see genotype::BOCData::outbreedTolerance
  DECLARE_PARAMETER(Bf, outbreedToleranceBounds)

  /// Mutation bounds for the sex
  /// \see genotype::BOCData::sex
  DECLARE_PARAMETER(Bounds<genotype::BOCData::Sex>, sexBounds)

  /// Mutation rates for the BOCData fields
  DECLARE_PARAMETER(MutationRates, mutationRates)
};

} // end of namespace config

namespace genotype {

/// Crossing of \p mother and \p father.
/// The algorithm is:
///   - Compute the distance based on the genomes and alignment
///   - Request a compatiblity rating \p r from the mother based on this distance
///   - Toss coin with success probability \p r
///     - if unsuccessfull, bail-out
///     - otherwise generate a child through the appropriate crossover algorithm
/// and potentially mutate it a bit.
///
/// \tparam GENOME The genome structure to cross
///
/// \param mother Female genome, requested for compatibility data
/// \param father Male genome
/// \param child Container for child genome in case of successfull mating.
/// Untouched otherwise
/// \param dice Source of randomness.
template <typename GENOME>
bool
bailOutCrossver(const GENOME &mother, const GENOME &father,
                     GENOME &child, rng::AbstractDice &dice) {

  double dist = distance(mother, father);
  assert(0 <= dist);

  double compat = mother.compatibility(dist);
  assert(0 <= compat && compat <= 1);

  if (dice(compat)) {
    child = cross(mother, father, dice);
    if (dice(config::SAGConfigFile<BOCData>::mutateChild())) child.mutate(dice);
    return true;
  }

  return false;
}

/// Crossing of \p mother and \p father.
/// The algorithm is:
///   - Compute the alignment between \p mother and \p father
///  (cache data, can be omitted for simple genomes)
///   - Compute the distance based on the genomes and alignment
///   - Request a compatiblity rating \p r from the mother based on this distance
///   - Toss coin with success probability \p r
///     - if unsuccessfull, bail-out
///     - otherwise generate a child through the appropriate crossover algorithm
/// and potentially mutate it a bit.
///
/// \tparam GENOME The genome structure to cross
///
/// \param mother Female genome, requested for compatibility data
/// \param father Male genome
/// \param child Container for child genome in case of successfull mating.
/// Untouched otherwise
/// \param dice Source of randomness.
template <typename GENOME>
std::enable_if_t<std::is_same<std::void_t<typename GENOME::Alignment>, void>::value, bool>
bailOutCrossver(const GENOME &mother, const GENOME &father,
                     GENOME &child, rng::AbstractDice &dice) {

  typename GENOME::Alignment alg = align(mother, father);

  double dist = distance(mother, father, alg);
  assert(0 <= dist);

  double compat = mother.compatibility(dist);
  assert(0 <= compat && compat <= 1);

  if (dice(compat)) {
    child = cross(mother, father, dice, alg);
    if (dice(config::SAGConfigFile<BOCData>::mutateChild())) child.mutate(dice);
    return true;
  }

  return false;
}

} // end of namespace genotype (again)

#endif // _CROSSOVER_HPP_
