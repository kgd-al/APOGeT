#ifndef _CROSSOVER_HPP_
#define _CROSSOVER_HPP_

#include <cmath>
#include <functional>
#include <atomic>

#include "kgd/external/json.hpp"

template <typename GENOME, typename CONFIG>
struct BailOutCrossover {
  using Alignment = typename GENOME::Alignment;

  class Data {
    double gaussoid (double x, double mu, double sigma) {
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

    using ID = uint;
    static constexpr ID NO_ID = ID(-1);
    static std::atomic<ID> NEXT_ID;
    ID id;

    enum Parent { MOTHER, FATHER };
    ID parents [2];

    uint generation;

    // ========================================================================
    // == Member function(s)

    double operator() (double distance) {
      return gaussoid(distance, optimalDistance,
                      distance < optimalDistance ? inbreedTolerance : outbreedTolerance);
    }

    // ========================================================================
    // == Genetic operators

    template <typename RND>
    void mutate (RND &rnd) {
      using CDM = CrossoverDataMutations;
      switch (rnd.pickOne(CONFIG::cdMutations())) {
      case CDM::DISTANCE: CONFIG::optimalDistance().mutate(optimalDistance, rnd); break;
      case CDM::INBREED:  CONFIG::inbreedTolerance().mutate(inbreedTolerance, rnd); break;
      case CDM::OUTBREED: CONFIG::outbreedTolerance().mutate(outbreedTolerance, rnd); break;
      }
    }

    /// TODO This function is ill-named
    template <typename RNG>
    void updateLineage (ID parent, RNG &rnd) {
      sex = rnd.toss(Sex::FEMALE, Sex::MALE);
      id = NEXT_ID++;
      parents[MOTHER] = parent;
      parents[FATHER] = -1;
      generation++;
    }

    friend double distance (const Data &lhs, const Data &rhs) {
      double d = 0;
      d += fabs(lhs.optimalDistance - rhs.optimalDistance) / CONFIG::optimalDistance().span();
      d += fabs(lhs.inbreedTolerance - rhs.inbreedTolerance) / CONFIG::inbreedTolerance().span();
      d += fabs(lhs.outbreedTolerance - rhs.outbreedTolerance) / CONFIG::outbreedTolerance().span();
      return d;
    }

    template <typename RND>
    friend Data crossover (const Data &lhs, const Data &rhs, RND &rnd) {
      Data child;
      child.optimalDistance = rnd.toss(lhs.optimalDistance, rhs.optimalDistance);
      child.inbreedTolerance = rnd.toss(lhs.inbreedTolerance, rhs.inbreedTolerance);
      child.outbreedTolerance = rnd.toss(lhs.outbreedTolerance, rhs.outbreedTolerance);
      child.sex = rnd.toss(lhs.sex, rhs.sex);

      child.id = NEXT_ID++;
      child.parents[MOTHER] = lhs.id;
      child.parents[FATHER] = rhs.id;
      child.generation = std::max(lhs.generation, rhs.generation) + 1;

      return child;
    }

    template <typename RND>
    static Data random (RND &rnd) {
      Data d;
      d.optimalDistance = CONFIG::optimalDistance().rand(rnd);
      d.inbreedTolerance = CONFIG::inbreedTolerance().rand(rnd);
      d.outbreedTolerance = CONFIG::outbreedTolerance().rand(rnd);
      d.sex = Sex(rnd(.5));

      d.id = NEXT_ID++;
      d.parents[MOTHER] = -1;
      d.parents[FATHER] = -1;
      d.generation = 0;

      return d;
    }

    // ========================================================================
    // == Serialization

    friend bool operator== (const Data &lhs, const Data &rhs) {
      return lhs.optimalDistance == rhs.optimalDistance
          && lhs.inbreedTolerance == rhs.inbreedTolerance
          && lhs.outbreedTolerance == rhs.outbreedTolerance
          && lhs.sex == rhs.sex
          && lhs.id == rhs.id
          && lhs.parents[MOTHER] == rhs.parents[MOTHER]
          && lhs.parents[FATHER] == rhs.parents[FATHER]
          && lhs.generation == rhs.generation;
    }

    friend void to_json (nlohmann::json &j, const Data &d) {
      j["distance"] = d.optimalDistance;
      j["inbreed"] = d.inbreedTolerance;
      j["outbreed"] = d.outbreedTolerance;
      j["sex"] = d.sex;

      j["id"] = d.id;
      j["mother"] = d.parents[MOTHER];
      j["father"] = d.parents[FATHER];
      j["gen"] = d.generation;
    }

    friend void from_json (const nlohmann::json &j, Data &d) {
      d.optimalDistance = j["distance"];
      d.inbreedTolerance = j["inbreed"];
      d.outbreedTolerance = j["outbreed"];
      d.sex = j["sex"];

      d.id = j["id"];
      d.parents[MOTHER] = j["mother"];
      d.parents[FATHER] = j["father"];
      d.generation = j["gen"];
    }

    // ========================================================================
    // == Debugging

    friend std::ostream& operator<< (std::ostream &os, const Data &d) {
      return os << "("
                  << d.optimalDistance
                  << ", " << d.inbreedTolerance
                  << ", " << d.outbreedTolerance
                  << ", " << (d.sex == Sex::MALE ? "M" : "F")
                << ")";
    }
  };

  template <typename RND>
  bool operator() (const GENOME &mother, const GENOME &father, GENOME &child, RND &rnd) {
    Alignment alg = align(mother, father);

    double dist = distance(mother, father, alg);
    assert(0 <= dist);

    double compat = mother.compatibility(dist);
    assert(0 <= compat && compat <= 1);

    if (rnd(compat)) {
      child = crossover(mother, father, alg, rnd);
      child.mutate(rnd);
      return true;
    }

    return false;
  }
};

template <typename G, typename C>
std::atomic<typename BailOutCrossover<G,C>::Data::ID> BailOutCrossover<G,C>::Data::NEXT_ID (0);

#endif // _CROSSOVER_HPP_
