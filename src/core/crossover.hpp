#ifndef _CROSSOVER_HPP_
#define _CROSSOVER_HPP_

#include <cmath>
#include <functional>
#include <atomic>

#include "kgd/external/json.hpp"

#include "crossconfig.h"

/*!
 * \file crossover.hpp
 *
 * Definition of structures and algorithms for implementing the Bail-Out Crossover
 */

namespace genotype {

/// Common crossover control data
class BOCData {
  /// Helper alias to the enumeration of possible mutations
  using CDM = config::CrossoverDataMutations;

  /// Helper alias to the configuration values
  using Config = config::Crossover;

  /// The configuration need privileged access to this structure's fields
  friend Config;

  /// The compatibility function. Half of an unnormalized gaussian.
  double gaussoid (double x, double mu, double sigma) const {
    return exp(-((x-mu)*(x-mu)/(2.*sigma*sigma)));
  }

  // ========================================================================
  // == Compatibility function

  /// Genetic distance that maximises reproduction compatibility
  float optimalDistance;

  /// Standard deviation for distances below optimal
  float inbreedTolerance;

  /// Standard deviation for distances above optimal
  float outbreedTolerance;

public:
  /// The possible sexs
  enum Sex { FEMALE, MALE };

  /// Which sex the associated genome codes for
  Sex sex;

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
  // == Genetic operators

  /// Mutate a single field of this object base on the mutation rates and bounds
  /// in #config::Crossover
  void mutate (rng::AbstractDice &dice) {
    switch (dice.pickOne(Config::cdMutations())) {
    case CDM::DISTANCE: Config::optimalDistance().mutate(*this, dice); break;
    case CDM::INBREED:  Config::inbreedTolerance().mutate(*this, dice); break;
    case CDM::OUTBREED: Config::outbreedTolerance().mutate(*this, dice); break;
    }
  }

  /// Updates internal data to reflect the clone status of the associated genome
  void updateCloneLineage (GID parent, rng::AbstractDice &dice) {
    sex = dice.toss(Sex::FEMALE, Sex::MALE);
    id = nextID();
    parents[MOTHER] = parent;
    parents[FATHER] = -1;
    generation++;
  }

  /// Compute the genetic distance for this structure
  friend double distance (const BOCData &lhs, const BOCData &rhs) {
    double d = 0;
    d += Config::optimalDistance().distance(lhs, rhs);
    d += Config::inbreedTolerance().distance(lhs, rhs);
    d += Config::outbreedTolerance().distance(lhs, rhs);
    return d;
  }

  /// Cross genetic material of both parent to create a new one
  /// Also update internal data (id, parents, generation)
  friend BOCData crossover (const BOCData &lhs, const BOCData &rhs,
                            rng::AbstractDice &dice) {
    BOCData child;
    dice.toss(lhs, rhs, child, &BOCData::optimalDistance);
    dice.toss(lhs, rhs, child, &BOCData::inbreedTolerance);
    dice.toss(lhs, rhs, child, &BOCData::outbreedTolerance);
    dice.toss(lhs, rhs, child, &BOCData::sex);

    child.id = nextID();
    child.parents[MOTHER] = lhs.id;
    child.parents[FATHER] = rhs.id;
    child.generation = std::max(lhs.generation, rhs.generation) + 1;

    return child;
  }

  /// \return a randomly generated structure according to the bounds in
  ///  #config::Crossover
  static BOCData random (rng::AbstractDice &dice) {
    BOCData d;
    d.optimalDistance = Config::optimalDistance().rand(dice);
    d.inbreedTolerance = Config::inbreedTolerance().rand(dice);
    d.outbreedTolerance = Config::outbreedTolerance().rand(dice);
    d.sex = Sex(dice(.5));

    d.id = nextID();
    d.parents[MOTHER] = -1;
    d.parents[FATHER] = -1;
    d.generation = 0;

    return d;
  }

  // ========================================================================
  // == Serialization

  /// Equality comparator. Compares each field for equality.
  friend bool operator== (const BOCData &lhs, const BOCData &rhs) {
    return lhs.optimalDistance == rhs.optimalDistance
        && lhs.inbreedTolerance == rhs.inbreedTolerance
        && lhs.outbreedTolerance == rhs.outbreedTolerance
        && lhs.sex == rhs.sex
        && lhs.id == rhs.id
        && lhs.parents[MOTHER] == rhs.parents[MOTHER]
        && lhs.parents[FATHER] == rhs.parents[FATHER]
        && lhs.generation == rhs.generation;
  }

  /// Converts a #genotype::BOCData structure into a field-named json
  friend void to_json (nlohmann::json &j, const BOCData &d) {
    j["mu"] = d.optimalDistance;
    j["si"] = d.inbreedTolerance;
    j["so"] = d.outbreedTolerance;
    j["S"] = d.sex;

    j["id"] = d.id;
    j["p0"] = d.parents[MOTHER];
    j["p1"] = d.parents[FATHER];
    j["G"] = d.generation;
  }

  /// Reads a #genotype::BOCData structure from a field-named json
  friend void from_json (const nlohmann::json &j, BOCData &d) {
    d.optimalDistance = j["mu"];
    d.inbreedTolerance = j["si"];
    d.outbreedTolerance = j["so"];
    d.sex = j["S"];

    d.id = j["id"];
    d.parents[MOTHER] = j["p0"];
    d.parents[FATHER] = j["p1"];
    d.generation = j["G"];
  }

  // ========================================================================
  // == Debugging

  /// Stream operator. Mostly for debugging purposes
  friend std::ostream& operator<< (std::ostream &os, const BOCData &d) {
    return os << "("
                << d.optimalDistance
                << ", " << d.inbreedTolerance
                << ", " << d.outbreedTolerance
                << ", " << (d.sex == Sex::MALE ? "M" : "F")
              << ")";
  }
};

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
/// \tparam Alignment Either the alignment structure provided by the genome or
/// one specified at call location
///
/// \param mother Female genome, requested for compatibility data
/// \param father Male genome
/// \param child Container for child genome in case of successfull mating.
/// Untouched otherwise
/// \param dice Source of randomness.
template <typename GENOME, typename Alignment = typename GENOME::Alignment>
bool bailOutCrossver (const GENOME &mother, const GENOME &father,
                      GENOME &child, rng::AbstractDice &dice) {

  Alignment alg = align(mother, father);

  double dist = distance(mother, father, alg);
  assert(0 <= dist);

  double compat = mother.compatibility(dist);
  assert(0 <= compat && compat <= 1);

  if (dice(compat)) {
    child = crossover(mother, father, dice, alg);
    if (dice(config::Crossover::mutateChild())) child.mutate(dice);
    return true;
  }

  return false;
}

} // end of namespace genotype

#endif // _CROSSOVER_HPP_
