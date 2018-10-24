#include "speciescontributors.h"
#include "../ptreeconfig.h"

namespace phylogeny {

auto debug = [] {
  return config::PTree::DEBUG_LEVEL() * config::PTree::DEBUG_CONTRIBUTORS();
};

bool operator== (const NodeContributor &lhs, const NodeContributor &rhs) {
  return lhs.speciesID() == rhs.speciesID()
      && lhs.count() == rhs.count();
}

SID Contributors::update (Contribution sids,
                          const ValidityEvaluator &elligible) {

  assert(nodeID != SID::INVALID);

  uint i=0;
  uint unprocessed = sids.size();

  if (debug() >= 1)
    std::cerr << "Updating contributions for " << nodeID << std::endl;

  // Ignore invalid(s)
  sids.erase(SID::INVALID);

  // Update already known contributors
  while(i < vec.size() && unprocessed > 0) {
    auto &c = vec[i];
    SID sid = c.speciesID();
    uint k = sids.count(sid);

    if (k > 0) {
      c += k;
      unprocessed -= k;
      sids.erase(sid);

      if (debug() >= 2)
        std::cerr << "\tAdded " << k << " at pos "
                  << i << " (SID=" << sid << ")" << std::endl;
    }

    i++;
  }

  // Register new contributors
  i = 0;
  while (!sids.empty()) {
    SID sid = *sids.begin();
    auto range = sids.equal_range(sid);
    uint k = std::distance(range.first, range.second);
    sids.erase(sid);

    bool e = elligible(nodeID, sid);
    vec.emplace_back(sid, k, e);

    if (debug() >= 2)
      std::cerr << "\tAppend " << k
                << " (SID=" << sid << ", elligible ? " << std::boolalpha << e
                << ")" << std::endl;
  }

  // sort by decreasing contribution
  std::stable_sort(vec.begin(), vec.end());

  return currentMain();
}

SID Contributors::currentMain (void) {
  assert(nodeID != SID::INVALID);

  if (vec.empty())
    return SID::INVALID;

  uint i = 1;
  auto mc = vec[0];
  while (!mc.elligible() && i < vec.size())
    mc = vec[i++];

  if (debug() >= 1)
    std::cerr << "Main contributor for " << nodeID << " is "
              << mc.speciesID() << " based on " << *this << std::endl;

  return mc.elligible() ? mc.speciesID() : SID::INVALID;
}

/// Serialize Contributors \p c into a json
void to_json (json &j, const Contributors &c) {
  j = {c.nodeID, c.vec};
}

/// Deserialize Contributors \p c from json \p j
void from_json (const json &j, Contributors &c) {
  uint i=0;
  c.nodeID = j[i++];
  c.vec = j[i++].get<decltype(Contributors::vec)>();
}

std::ostream& operator<< (std::ostream &os, const Contributors &c) {
  os << "[ ";
  for (const NodeContributor &nc: c.vec)
    os << "{" << nc.speciesID() << "," << nc.count() << "} ";
  return os << "]";
}

} // end of namespace phylogeny
