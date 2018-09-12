#include "crossover.h"

#define GENOME BOCData

DEFINE_GENOME_FIELD_WITH_BOUNDS(float, optimalDistance, "mu", 0.f, true, 4.f, 100.f)
DEFINE_GENOME_FIELD_WITH_BOUNDS(float, inbreedTolerance, "si", 0.f, 2.f, 2.f, 10.f)
DEFINE_GENOME_FIELD_WITH_BOUNDS(float, outbreedTolerance, "so", 0.f, 2.f, 2.f, 10.f)

using Sex = genotype::BOCData::Sex;
DEFINE_GENOME_FIELD_WITH_BOUNDS(Sex, sex, "S", Sex::FEMALE, Sex::MALE)

std::ostream& operator<< (std::ostream &os, Sex s) {
  switch (s) {
  case Sex::FEMALE: return os << "F";
  case Sex::MALE: return os << "M";
  }
  os.setstate(std::ios_base::failbit);
  return os << "ERROR";
}

std::istream& operator>> (std::istream &is, Sex &s) {
  char c;
  is >> c;
  switch (c) {
  case 'F': s = Sex::FEMALE; break;
  case 'M': s = Sex::MALE;   break;
  }
  is.setstate(std::ios_base::failbit);
  return is;
}

DEFINE_GENOME_MUTATION_RATES({
  MUTATION_RATE(optimalDistance, 1.f),
  MUTATION_RATE(inbreedTolerance, 1.f),
  MUTATION_RATE(outbreedTolerance, 1.f),
//  MUTATION_RATE(sex, 1.f),
})

#undef GENOME

#define CFILE config::SAGConfigFile<genotype::BOCData>

DEFINE_PARAMETER(float, mutateChild, .5)

#undef CFILE
