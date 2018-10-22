#include "speciescontributors.h"
#include "../ptreeconfig.h"

namespace phylogeny {

using Config = config::PTree;

SID Contributors::update (Contribution sids) {
  assert(nodeID != SID::INVALID);

  uint i=0;
  uint processed = 0;

  // Ignore invalid(s)
  sids.erase(SID::INVALID);

  // Update already known contributors
  while(i < vec.size() && processed < sids.size()) {
    SID sid = vec[i].speciesID();
    uint k = sids.count(sid);

    if (k > 0) {
      vec[i] += k;
      processed += k;
      sids.erase(sid);
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

    vec.emplace_back(sid, k);
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
  SID sid = vec[0].speciesID();
  while (sid == nodeID && i < vec.size())
    sid = vec[i++].speciesID();

  if (Config::DEBUG() >= 2)
    std::cerr << "Main contributor for " << nodeID << " is "
              << sid << " based on " << vec << std::endl;

  if (sid == nodeID)
        return SID::INVALID;
  else  return sid;
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

} // end of namespace phylogeny
