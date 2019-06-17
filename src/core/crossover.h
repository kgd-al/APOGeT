#ifndef _CROSSOVER_HPP_
#define _CROSSOVER_HPP_

#include <cmath>
#include <functional>
#include <atomic>

#include "kgd/external/json.hpp"

#include "kgd/utils/functions.h"

#include "kgd/genotype/selfawaregenome.hpp"

#include "tree/treetypes.h"

/*!
 * \file crossover.h
 *
 * Definition of structures and algorithms for implementing the Bail-Out Crossover
 */

namespace genotype {

/// Common crossover control data
class BOCData : public EDNA<BOCData> {
  APT_EDNA()

  // ========================================================================
  // == Compatibility function

  /// The compatibility function. Two halves of an unnormalized gaussian.
  static double gaussoid (double d, double mu, double sigma) {
    return utils::gauss(d, mu, sigma);
  }

  /// The inverse compatibility function. Returns the distances for this
  ///  compatibility value
  static double gaussoid_inverse (double c, double mu, double sigma, int sign) {
    return std::max(0., utils::gauss_inverse(c, mu, sigma, sign));
  }

  /// Genetic distance that maximises reproduction compatibility
  float optimalDistance;

  /// Standard deviation for distances below optimal
  float inbreedTolerance;

  /// Standard deviation for distances above optimal
  float outbreedTolerance;

public:
  /// Helper alias to the source of randomness
  using Dice = EDNA<BOCData>::Dice;

  /// The possible sexes
  enum Sex { FEMALE, MALE };

  /// Which sex the associated genome codes for
  Sex sex;

  // ========================================================================
  // == Member functions

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

  /// Evalutes the distance between this genome and another that produced the
  /// compatibility value \p compat
  ///
  /// This is the inverse operation of operator()(double)
  ///
  /// \param compat the compatibility
  /// \param d_inbreed the distance for inbreeding
  /// \param d_outbreed the distance for outbreeding
  void operator() (double compat, double &d_inbreed, double &d_outbreed) const {
    assert(0 <= compat && compat <= 1);
    d_inbreed = gaussoid_inverse(compat, optimalDistance, inbreedTolerance, -1);
    d_outbreed = gaussoid_inverse(compat, optimalDistance, outbreedTolerance, 1);
  }
};

DECLARE_GENOME_FIELD(BOCData, float, optimalDistance)
DECLARE_GENOME_FIELD(BOCData, float, inbreedTolerance)
DECLARE_GENOME_FIELD(BOCData, float, outbreedTolerance)
DECLARE_GENOME_FIELD(BOCData, BOCData::Sex, sex)

/// Pretty prints a sex enumeration value
std::ostream& operator<< (std::ostream &os, BOCData::Sex s);

/// Converts a pretty printed sex enumeration value to its integer representation
std::istream& operator>> (std::istream &is, BOCData::Sex &s);

} // end of namespace genotype

namespace config {

/// Config file for the crossover algorithms
template <> struct EDNA_CONFIG_FILE(BOCData) {
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

  /// Distance weights for the BOCData fields
  DECLARE_PARAMETER(DistanceWeights, distanceWeights)
};

/// Specialization for the sex field which is not included in the distances
template <>
struct MutationSettings::BoundsOperators<genotype::BOCData::Sex, void> {
  using Sex = genotype::BOCData::Sex; ///< Helper alias to the sex enumeration
  using Dice = MutationSettings::Dice; ///< Helper alias to the rng

  /// \returns either Male or Female
  static Sex rand (const Sex&, const Sex&, Dice &dice) {
    return dice.toss(Sex::MALE, Sex::FEMALE);
  }

  /// \returns 0
  static double distance (const Sex&, const Sex&, const Sex&, const Sex&) {
    return 0;
  }

  /// \returns flip the sex
  static void mutate (Sex &s, const Sex&, const Sex&, Dice&) {
    s = (s == Sex::MALE) ? Sex::FEMALE : Sex::MALE;
  }

  /// \returns whether the sex is either male or female
  static bool check (Sex &s, const Sex&, const Sex&) {
    return (s == Sex::MALE) || (s == Sex::FEMALE);
  }
};

} // end of namespace config

namespace genotype {

namespace _details {
template <typename T, typename = void>
struct requiresAlignment : std::false_type {};

template <typename T>
struct requiresAlignment<T, std::void_t<typename T::Alignment>> : std::true_type {};
}

/// Crossing of \p mother and \p father.
/// The algorithm is:
///   - Compute the genomic distqnce
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
/// \param litter Container for children genome. In case of successfull mating
/// will contains as many offsprings as the input size. Untouched otherwise.
/// \param dice Source of randomness.
/// \param oDistance If not null, will be filled with the genomic distance
/// \param oCompatibility If not null, will be filled with perceived compatibility
template <typename GENOME>
std::enable_if_t<!_details::requiresAlignment<GENOME>::value, bool>
bailOutCrossver(const GENOME &mother, const GENOME &father,
                std::vector<GENOME> &litter, rng::AbstractDice &dice,
                float *oDistance = nullptr, float *oCompatibility = nullptr) {

  double dist = distance(mother, father);
  assert(0 <= dist);
  if (oDistance) *oDistance = dist;

  double compat = mother.compatibility(dist);
  assert(0 <= compat && compat <= 1);
  if (oCompatibility)  *oCompatibility = compat;

  if (dice(compat)) {
    for (GENOME &child: litter) {
      child = cross(mother, father, dice);
      if (dice(config::EDNAConfigFile<BOCData>::mutateChild()))
        child.mutate(dice);
    }
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
/// \todo FIXME Forwarding to auto-field managers is a pain...
///
/// \tparam GENOME The genome structure to cross
///
/// \param mother Female genome, requested for compatibility data
/// \param father Male genome
/// \param litter Container for children genome. In case of successfull mating
/// will contains as many offsprings as the input size. Untouched otherwise.
/// \param dice Source of randomness.
/// \param oDistance If not null, will be filled with the genomic distance
/// \param oCompatibility If not null, will be filled with perceived compatibility
template <typename GENOME>
std::enable_if_t<_details::requiresAlignment<GENOME>::value, bool>
bailOutCrossver(const GENOME &mother, const GENOME &father,
                std::vector<GENOME> &litter, rng::AbstractDice &dice,
                float *oDistance = nullptr, float *oCompatibility = nullptr) {

  typename GENOME::Alignment alg = align(mother, father);

  double dist = distance(mother, father, alg);
  assert(0 <= dist);
  if (oDistance) *oDistance = dist;

  double compat = mother.compatibility(dist);
  assert(0 <= compat && compat <= 1);
  if (oCompatibility)  *oCompatibility = compat;

  if (dice(compat)) {
    for (GENOME &child: litter) {
      child = cross(mother, father, dice, alg);
      if (dice(config::EDNAConfigFile<BOCData>::mutateChild()))
        child.mutate(dice);
    }
    return true;
  }

  return false;
}

} // end of namespace genotype (again)

#endif // _CROSSOVER_HPP_
