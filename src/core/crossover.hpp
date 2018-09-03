#ifndef _CROSSOVER_HPP_
#define _CROSSOVER_HPP_

#include <cmath>
#include <functional>
#include <atomic>

#include "kgd/external/json.hpp"

#include "crossconfig.h"

namespace genotype {

class BOCData {
  using CDM = config::CrossoverDataMutations;
  using Config = config::Crossover;
  friend Config;

  double gaussoid (double x, double mu, double sigma) const {
    return exp(-((x-mu)*(x-mu)/(2.*sigma*sigma)));
  }

  // ========================================================================
  // == Compatibility function

  float optimalDistance;
  float inbreedTolerance;
  float outbreedTolerance;

public:
  enum Sex { FEMALE, MALE };
  Sex sex;

  // ========================================================================
  // == Classification data

  using GID = uint;
  static constexpr GID NO_ID = GID(-1);

  GID id;
  static GID nextID (void) {
    static std::atomic<GID> NEXT_ID (0);
    return NEXT_ID++;
  }

  enum Parent : uint { MOTHER = 0, FATHER = 1 };
  GID parents [2];

  uint generation;

  // ========================================================================
  // == Member function(s)

  bool hasParent (Parent p) const {
    return parents[p] != NO_ID;
  }

  GID parent (Parent p) const {
    return parents[p];
  }

  float getOptimalDistance (void) const { return optimalDistance; }
  float getInbreedTolerance (void) const { return inbreedTolerance; }
  float getOutbreedTolerance (void) const { return outbreedTolerance; }

  double operator() (double distance) const {
    return gaussoid(distance, optimalDistance,
                    distance < optimalDistance ? inbreedTolerance : outbreedTolerance);
  }

  // ========================================================================
  // == Genetic operators

  template <typename RND>
  void mutate (RND &rnd) {
    switch (rnd.pickOne(Config::cdMutations())) {
    case CDM::DISTANCE: Config::optimalDistance().mutate(this, rnd); break;
    case CDM::INBREED:  Config::inbreedTolerance().mutate(this, rnd); break;
    case CDM::OUTBREED: Config::outbreedTolerance().mutate(this, rnd); break;
    }
  }

  /// \todo This function is ill-named
  template <typename RNG>
  void updateLineage (GID parent, RNG &rnd) {
    sex = rnd.toss(Sex::FEMALE, Sex::MALE);
    id = nextID();
    parents[MOTHER] = parent;
    parents[FATHER] = -1;
    generation++;
  }

  friend double distance (const BOCData &lhs, const BOCData &rhs) {
    double d = 0;
    d += Config::optimalDistance().distance(lhs, rhs);
    d += Config::inbreedTolerance().distance(lhs, rhs);
    d += Config::outbreedTolerance().distance(lhs, rhs);
    return d;
  }

  template <typename RND>
  friend BOCData crossover (const BOCData &lhs, const BOCData &rhs, RND &rnd) {
    BOCData child;
    rnd.toss(lhs, rhs, child, &BOCData::optimalDistance);
    rnd.toss(lhs, rhs, child, &BOCData::inbreedTolerance);
    rnd.toss(lhs, rhs, child, &BOCData::outbreedTolerance);
    rnd.toss(lhs, rhs, child, &BOCData::sex);

    child.id = nextID();
    child.parents[MOTHER] = lhs.id;
    child.parents[FATHER] = rhs.id;
    child.generation = std::max(lhs.generation, rhs.generation) + 1;

    return child;
  }

  template <typename RND>
  static BOCData random (RND &rnd) {
    BOCData d;
    d.optimalDistance = Config::optimalDistance().rand(rnd);
    d.inbreedTolerance = Config::inbreedTolerance().rand(rnd);
    d.outbreedTolerance = Config::outbreedTolerance().rand(rnd);
    d.sex = Sex(rnd(.5));

    d.id = nextID();
    d.parents[MOTHER] = -1;
    d.parents[FATHER] = -1;
    d.generation = 0;

    return d;
  }

  // ========================================================================
  // == Serialization

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

  friend std::ostream& operator<< (std::ostream &os, const BOCData &d) {
    return os << "("
                << d.optimalDistance
                << ", " << d.inbreedTolerance
                << ", " << d.outbreedTolerance
                << ", " << (d.sex == Sex::MALE ? "M" : "F")
              << ")";
  }
};


template <typename GENOME, typename Alignment = typename GENOME::Alignment>
struct BailOutCrossover {
  using Data = BOCData;

  template <typename RND>
  bool operator() (const GENOME &mother, const GENOME &father, GENOME &child, RND &rnd) {
    Alignment alg = align(mother, father);

    double dist = distance(mother, father, alg);
    assert(0 <= dist);

    double compat = mother.compatibility(dist);
    assert(0 <= compat && compat <= 1);

    if (rnd(compat)) {
      child = crossover(mother, father, rnd, alg);
      if (rnd(config::Crossover::mutateChild())) child.mutate(rnd);
      return true;
    }

    return false;
  }
};

} // end of namespace genotype

#endif // _CROSSOVER_HPP_
