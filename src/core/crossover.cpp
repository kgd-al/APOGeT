#include "crossover.h"

std::atomic<genotype::BOCData::GID_ut> genotype::BOCData::NEXT_ID = 0;

#define GENOME BOCData

DEFINE_GENOME_FIELD_WITH_BOUNDS(float, optimalDistance, "mu", 0.f, true, 4.f, 100.f)
DEFINE_GENOME_FIELD_WITH_BOUNDS(float, inbreedTolerance, "si", 0.f, 2.f, 2.f, 10.f)
DEFINE_GENOME_FIELD_WITH_BOUNDS(float, outbreedTolerance, "so", 0.f, 2.f, 2.f, 10.f)

using Sex = genotype::BOCData::Sex;
DEFINE_GENOME_FIELD_WITH_BOUNDS(Sex, sex, "S", Sex::FEMALE, Sex::MALE)

std::ostream& genotype::operator<< (std::ostream &os, Sex s) {
  switch (s) {
  case Sex::FEMALE: return os << "F";
  case Sex::MALE: return os << "M";
  }
  os.setstate(std::ios_base::failbit);
  return os << "ERROR";
}

std::istream& genotype::operator>> (std::istream &is, Sex &s) {
  char c;
  is >> c;
  switch (c) {
  case 'F': s = Sex::FEMALE; break;
  case 'M': s = Sex::MALE;   break;
  default:  is.setstate(std::ios_base::failbit);
  }
  return is;
}

DEFINE_GENOME_MUTATION_RATES({
  __EDNA_PAIR(  optimalDistance, 1.f),
  __EDNA_PAIR( inbreedTolerance, 1.f),
  __EDNA_PAIR(outbreedTolerance, 1.f),
  __EDNA_PAIR(              sex, 1.f),
})

DEFINE_GENOME_DISTANCE_WEIGHTS({
  __EDNA_PAIR(  optimalDistance, 1.f),
  __EDNA_PAIR( inbreedTolerance, 1.f),
  __EDNA_PAIR(outbreedTolerance, 1.f),
  __EDNA_PAIR(              sex, 1.f),
})

#undef GENOME

#define CFILE config::EDNAConfigFile<genotype::BOCData>

DEFINE_PARAMETER(float, mutateChild, .5)

#undef CFILE
