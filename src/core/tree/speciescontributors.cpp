#include "speciescontributors.h"
#include "../ptreeconfig.h"

namespace phylogeny {

auto debug = [] {
  return config::PTree::DEBUG_LEVEL() * config::PTree::DEBUG_CONTRIBUTORS();
};

bool operator== (const Contributor &lhs, const Contributor &rhs) {
  return lhs.speciesID() == rhs.speciesID()
      && lhs.count() == rhs.count();
}

SID Contributors::update (Contributions ctbs,
                          const ValidityEvaluator &elligible) {

  assert(nodeID != SID::INVALID);

  uint i=0;
  uint unprocessed = ctbs.size();

  if (debug() >= 1)
    std::cerr << "Updating contributions for " << nodeID << std::endl;

  // Ignore invalid(s)
  ctbs.erase(std::remove(ctbs.begin(), ctbs.end(), SID::INVALID),
             ctbs.end());

  // Update already known contributors
  while(i < vec.size() && unprocessed > 0) {
    auto &ctor = vec[i];
    SID sid = ctor.speciesID();
    auto it = std::find(ctbs.begin(), ctbs.end(), sid);

    if (it != ctbs.end()) {
      Contribution &ction = *it;
      ctor += ction.count;
      unprocessed--;
      ctbs.erase(it);

      if (debug() >= 2)
        std::cerr << "\tAdded " << ction.count << " at pos "
                  << i << " (SID=" << sid << ")" << std::endl;
    }

    i++;
  }

  // Register new contributors
  for (Contribution &ction: ctbs) {
    bool e = elligible(nodeID, ction.species);
    vec.emplace_back(ction.species, ction.count, e);

    if (debug() >= 2)
      std::cerr << "\tAppend " << ction.count
                << " (SID=" << ction.species << ", elligible ? "
                << std::boolalpha << e << ")" << std::endl;
  }

  // sort by decreasing contribution
  std::stable_sort(vec.begin(), vec.end(), Contributor::CMP());

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

SID Contributors::updateElligibilities(const ValidityEvaluator &elligible) {
  for (Contributor &c: vec)
    c.setElligible(elligible(nodeID, c.speciesID()));

  return currentMain();
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
  for (const Contributor &nc: c.vec)
    os << "{" << nc.speciesID() << "," << nc.count() << "} ";
  return os << "]";
}

} // end of namespace phylogeny
