#include <iomanip>

#include "../ptreeconfig.h"
#include "treetypes.h"

namespace phylogeny {
namespace _details {

using Config = config::PTree;

auto debug = [] {
  return Config::DEBUG_LEVEL() * Config::DEBUG_ENVELOPPE();
};

/// \returns a vector of indices into \p values allowing sorted access
template <typename T>
static std::vector<size_t> ordered(std::vector<T> const& values) {
  std::vector<size_t> indices(values.size());
  std::iota(begin(indices), end(indices), static_cast<size_t>(0));

  std::sort(
    begin(indices), end(indices),
    [&values](size_t a, size_t b) { return values[a] > values[b]; }
  );
  return indices;
}

void computeAvgAndStdDev (const DistanceMap &m, double &avg, double &stdDev) {
  avg = 0;
  for (auto &it: m) avg += it.second;
  avg /= double(m.size());

  stdDev = 0;
  for (auto &it: m) stdDev += std::pow(avg - it.second, 2);
  stdDev = std::sqrt(stdDev / double(m.size()));
}


// Maximize average (has a known pitfall)
EnveloppeContribution maxAverage (const DistanceMap &edist,
                                  const std::vector<float> &gdist,
                                  GID gid, const std::vector<GID> &ids) {

  const uint k = ids.size();
  EnveloppeContribution ec;
  ec.value = -std::numeric_limits<double>::max();
  ec.than = -1;
  ec.better = false;

  const auto pad = [gid] {
    return std::setw(ceil(log10(std::underlying_type<GID>::type(gid))));
  };

  // Compute variance contributions and least contributor
  for (uint i=0; i<k; i++) {
    if (debug() >= 2)
      std::cerr << "\n\t\tc(" << pad() << ids[i]
                << "/" << pad() << gid << ") =";

    double c = 0;
    for (uint j=0; j<k; j++) {
      if (i==j) continue;
      c += - edist.at({i,j}) + gdist.at(j);

      if (debug() >= 2) {
        std::cerr << std::left;
        if (j>0)  std::cerr << "\t\t  " << pad() << " "
                            << " " << pad() << " " << "   ";
        std::cerr << " - " << std::setw(8) << edist.at({i,j})
                  << " + " << std::setw(8) << gdist.at(j);
        if (j<k-2)  std::cerr << "\n";
      }
    }

    if (debug() >= 2) std::cerr << " = " << c << std::endl;
    if (ec.value < c) {
      ec.value = c;
      ec.than = i;
    }
  }

  ec.better = (ec.value > 0);
  return ec;
}

// Just maximize min distance
EnveloppeContribution maxMinDist (const DistanceMap &edist,
                                  const std::vector<float> &gdist,
                                  GID gid, const std::vector<GID> &ids) {

  const uint k = ids.size();

  EnveloppeContribution ec;
  ec.value = -std::numeric_limits<double>::max();
  ec.than = -1;
  ec.better = false;

  const auto pad = [gid] {
    return std::setw(ceil(log10(std::underlying_type<GID>::type(gid))));
  };

  // Compare with each vertex
  for (uint i=0; i<k; i++) {
    if (debug() >= 2)
      std::cerr << "\t\tc(" << pad() << ids[i]
                << "/" << pad() << gid << ") =" << std::left;

    float minBase = std::numeric_limits<float>::max(),
          minNew = std::numeric_limits<float>::max();

    for (uint j=0; j<k; j++) {
      if (i==j) continue;
      minBase = std::min(minBase, edist.at({i,j}));
      minNew = std::min(minNew, gdist.at(j));
    }

    double c = - minBase + minNew;

    if (debug() >= 2)
      std::cerr << " - " << std::setw(8) << minBase
                << " + " << std::setw(8) << minNew
                << " = " << std::setw(8) << c
                << std::endl;

    if (ec.value < c) {
      ec.value = c;
      ec.than = i;
    }
  }

  ec.better = (ec.value > 0);
  return ec;
}

// Maximize mean distance while reducing deviation
EnveloppeContribution maxAvgMinStdDev (const DistanceMap &edist,
                                       const std::vector<float> &gdist,
                                       GID gid, const std::vector<GID> &ids) {

  const uint k = ids.size();

  EnveloppeContribution ec;
  ec.value = -std::numeric_limits<double>::max();
  ec.than = -1;
  ec.better = false;

  const auto pad = [gid] {
    return std::setw(ceil(log10(std::underlying_type<GID>::type(gid))));
  };

  // Average internal distance
  double baseAVG, baseStdDev;
  computeAvgAndStdDev(edist, baseAVG, baseStdDev);

  // Compare with each vertex
  for (uint i=0; i<k; i++) {
    if (debug() >= 2)
      std::cerr << "\t\tc(" << pad() << ids[i]
                << "/" << pad() << gid << ") =" << std::left;

    DistanceMap newMap = edist;
    for (uint j=0; j<k; j++)
      if (i != j)
        newMap.at({i,j}) = gdist[j];

    double newAVG, newStdDev;
    computeAvgAndStdDev(newMap, newAVG, newStdDev);

    double c = - baseAVG + newAVG
               + baseStdDev - newStdDev;

    if (debug() >= 2)
      std::cerr << " - " << std::setw(8) << baseAVG
                << " + " << std::setw(8) << newAVG
                << " + " << std::setw(8) << baseStdDev
                << " - " << std::setw(8) << newStdDev
                << " = " << std::setw(8) << c
                << std::endl;

//      if (debug() >= 2) {
//        std::cerr << std::left;
//        if (j>0)  std::cerr << "\t\t  " << pad() << " "
//                            << " " << pad() << " " << "   ";
//        std::cerr << "\t"
//                  << std::setw(8) << w << " * (";
//        std::cerr << std::setw(9) << nc
//                  << " + " << std::setw(8) << pc
//                  << ")";
//        if (j<k-2)  std::cerr << "\n";
//      }
//    }
//    if (debug() >= 2) std::cerr << " = " << c << std::endl;

    if (ec.value < c) {
      ec.value = c;
      ec.than = i;
    }
  }

  ec.better = (ec.value > 0);
  return ec;
}

// Weighted by distance to mean. Shitty
EnveloppeContribution maxWeightedDist2Avg (const DistanceMap &edist,
                                           const std::vector<float> &gdist,
                                           GID gid, const std::vector<GID> &ids) {

  const uint k = ids.size();
  EnveloppeContribution ec;
  ec.value = -std::numeric_limits<double>::max();
  ec.than = -1;
  ec.better = false;

  const auto pad = [gid] {
    return std::setw(ceil(log10(std::underlying_type<GID>::type(gid))));
  };

  // Average internal distance
  double A = 0;
  for (auto &it: edist)  A += it.second;
  A /= double(edist.size());

  //
  auto weight = [A] (double d) {
    return 1 - exp(- (d-A)*(d-A) / (2. * A * A / 16.));
  };

  // Compute variance contributions and least contributor
  for (uint i=0; i<k; i++) {
    std::vector<double> d_i (k), d_g(k);
    if (debug() >= 2)
      std::cerr << "\n\t\tc(" << pad() << ids[i]
                << "/" << pad() << gid << ") =";

    for (uint j=0; j<k; j++) {
      if (i == j) continue;
      d_i.push_back(edist.at({i,j}));
      d_g.push_back(gdist[j]);
    }

    double c = 0;
    const auto i_i = ordered(d_i), i_g = ordered(d_g);
    for (uint j=0; j<k-1; j++) {
      double nc = - d_i[i_i[j]];
      double pc = + d_g[i_g[j]];
      double w = weight(pc);
      c += w * (nc + pc);

      if (debug() >= 2) {
        std::cerr << std::left;
        if (j>0)  std::cerr << "\t\t  " << pad() << " "
                            << " " << pad() << " " << "   ";
        std::cerr << "\t"
                  << std::setw(8) << w << " * (";
        std::cerr << std::setw(9) << nc
                  << " + " << std::setw(8) << pc
                  << ")";
        if (j<k-2)  std::cerr << "\n";
      }
    }

    if (debug() >= 2) std::cerr << " = " << c << std::endl;
    if (ec.value < c) {
      ec.value = c;
      ec.than = i;
    }
  }

  ec.better = (ec.value > 0);
  return ec;
}


EnveloppeContribution computeContribution (const DistanceMap &edist,
                                           const std::vector<float> &gdist,
                                           GID gid, const std::vector<GID> &ids) {
  auto f = computeContribution;
  switch (Config::DEBUG_ENV_CRIT()) {
  case 0: f = maxAverage;           break;
  case 1: f = maxMinDist;           break;
  case 2: f = maxAvgMinStdDev;      break;
  case 3: f = maxWeightedDist2Avg;  break;
  default:
    throw std::logic_error("No function for this use case");
  }
  assert(f != computeContribution);
  return f(edist, gdist, gid, ids);
}
}
}
