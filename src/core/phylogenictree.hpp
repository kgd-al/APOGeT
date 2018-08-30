#ifndef _PHYLOGENIC_TREE_H_
#define _PHYLOGENIC_TREE_H_

#include <vector>
#include <map>
#include <memory>
#include <fstream>

#include <cassert>
#include <iostream>

#include "ptreeconfig.h"
#include "kgd/external/json.hpp"

namespace phylogeny {

template <typename PT>
struct Callbacks_t {
  using GID = typename PT::GID;
  using SID = typename PT::SID;
  using LivingSet = typename PT::LivingSet;

  void onStepped (uint step, const LivingSet &living);
  void onNewSpecies (SID pid, SID sid);
  void onGenomeEntersEnveloppe (SID sid, GID gid);
  void onGenomeLeavesEnveloppe (SID sid, GID gid);
};

struct SpeciesData {
  uint firstAppearance;
  uint lastAppearance;
  uint count;

  friend void to_json (nlohmann::json &j, const SpeciesData &d) {
    j = {d.firstAppearance, d.lastAppearance, d.count};
  }

  friend void from_json (const nlohmann::json &j, SpeciesData &d) {
    uint i=0;
    d.firstAppearance = j[i++];
    d.lastAppearance = j[i++];
    d.count = j[i++];
  }
};

template <typename GENOME>
class PhylogenicTree {
public:
  using GID = typename GENOME::CData::GID;
  using SID = uint;

  using Parent = typename GENOME::CData::Parent;
  static constexpr SID NoID = SID(-1);

  using LivingSet = std::set<SID>;

  using Callbacks = Callbacks_t<PhylogenicTree<GENOME>>;


  PhylogenicTree(void) {
    _nextNodeID = 0;
    _step = 0;
    _hybrids = 0;
    _root = nullptr;
    _callbacks = nullptr;
  }

  void setCallbacks (Callbacks *c) const { _callbacks = c; }
  Callbacks* callbacks (void) {   return _callbacks; }

  const auto& root (void) const {
    return _root;
  }

  uint width (void) const {
    return _nodes.size();
  }

  const auto& nodeAt (uint i) const {
    return _nodes.at(i);
  }

  uint step (void) const {
    return _step;
  }

  void setStep (uint step) {
    _step = step;
  }

  template <typename IT, typename F>
  void step (uint step, IT begin, IT end, F geneticID) {
    LivingSet aliveSpecies;
    for (IT it = begin; it != end; ++it)
      aliveSpecies.insert(_idToSpecies.at(geneticID(*it)));

    for (SID sid: aliveSpecies)
      _nodes.at(sid)->data.lastAppearance = step;

    _step = step;
    if (_callbacks) _callbacks->onStepped(step, aliveSpecies);

    if (PTreeConfig::DEBUG())
      std::cerr << _idToSpecies.size() << " id>species pairs stored" << std::endl;
  }

  SID addGenome (const GENOME &g) {
    if (!_root)
      _root = makeNode(nullptr);

    if (!g.hasParent(Parent::FATHER) || !g.hasParent(Parent::MOTHER))
      return addGenome(g, _root);

    uint mSID = _idToSpecies.parentSID(g, Parent::MOTHER),
         pSID = _idToSpecies.parentSID(g, Parent::FATHER);

    assert(PTreeConfig::ignoreHybrids() || mSID == pSID);
    if (mSID != pSID) _hybrids++;
    if (mSID == pSID)
      return addGenome(g, _nodes[mSID]);

    else if (PTreeConfig::ignoreHybrids()) {
      if (PTreeConfig::DEBUG() >= 0)
        std::cerr << "Linking hybrid genome " << g.id() << " to mother species" << std::endl;

      return addGenome(g, _nodes[mSID]);

    } else {
      assert(false);
      if (PTreeConfig::DEBUG())
        std::cerr << "Managing hybrid genome " << g.id() << std::endl;
    }

    return NoID;
  }

  void delGenome (const GENOME &g) {
    auto sid = _idToSpecies.remove(g);
    if (PTreeConfig::DEBUG())
      std::cerr << "New last appearance of species " << sid << " is " << _step << std::endl;

    _nodes[sid]->data.lastAppearance = _step;
  }

  friend std::ostream& operator<< (std::ostream &os, const PhylogenicTree &pt) {
    os << pt._hybrids << " Hybrids;\n";
    return os << *pt._root;
  }

  void logTo (const std::string &filename) const {
    std::ofstream ofs (filename);
    ofs << "digraph {\n";
    _root->logTo(ofs);
    ofs << "}\n";
  }

protected:
  SID _nextNodeID;

  template <typename T>
  struct ordered_pair {
      T first, second;
      ordered_pair(T first, T second) : first(std::min(first, second)), second(std::max(first, second)) {}

      friend bool operator< (const ordered_pair &lhs, const ordered_pair &rhs) {
        if (lhs.first != rhs.first) return lhs.first < rhs.first;
        return lhs.second < rhs.second;
      }
  };

  struct Node;
  using Node_ptr = std::shared_ptr<Node>;
  struct Node {
    SID id;
    SpeciesData data;

    Node *parent;
    std::vector<Node_ptr> children;

    std::vector<GENOME> enveloppe;
    std::map<ordered_pair<uint>, float> distances;

    Node (SID id, Node_ptr parent) : id(id), parent(parent.get()) {}

    friend std::ostream& operator<< (std::ostream &os, const Node &n) {
      std::string spacing = "> ";
      const Node *p = &n;
      while ((p = p->parent))  spacing += "  ";

      os << spacing << "[" << n.id << "] ( ";
      for (const GENOME &g: n.enveloppe)    os << g.id() << " ";
      os << ")\n";

      for (const Node_ptr &ss: n.children)  os << *ss.get();

      return os;
    }

    void logTo (std::ostream &os) const {
      os << "\t" << id << ";\n";
      for (const Node_ptr &n: children) {
        os << "\t" << id << " -> " << n->id << ";\n";
        n->logTo(os);
      }
    }
  };

  struct DCCache {
    std::vector<float> distances, compatibilities;

    void reserve (uint n) {
      distances.reserve(n), compatibilities.reserve(n);
    }

    void push_back (float d, float c) {
      distances.push_back(d), compatibilities.push_back(c);
    }
  };

  struct IdToSpeciesMap {
    struct ITSMData {
      SID species;
      uint refCount;
    };
    std::map<GID, ITSMData> map;

    uint size (void) const {
      return map.size();
    }

    SID parentSID (const GENOME &g, Parent p) {
      auto &d = map.at(g.parent(p));
      d.refCount++;
      return d.species;
    }

    SID remove (const GID id) {
      auto it = map.find(id);
      auto sid = it->second.species;
      auto &ref = it->second.refCount;

      assert(ref > 0);
      ref--;

      if (ref == 0)  map.erase(it);

      return sid;
    }

    SID remove (const GENOME &g) {
      auto sid = remove(g.id());
      for (Parent p: {Parent::MOTHER, Parent::FATHER})
        if (g.hasParent(p))
          remove(g.parent(p));
      return sid;
    }

    void insert (GID gid, SID sid) {
      ITSMData d;
      d.refCount = 1;
      d.species = sid;
      map[gid] = d;
    }

    SID at (GID gid) const {
      return map.at(gid).species;
    }
  };

  Node_ptr _root;

  std::vector<Node_ptr> _nodes;

  IdToSpeciesMap _idToSpecies;

  mutable Callbacks *_callbacks;

  uint _hybrids;
  uint _step;

  Node_ptr makeNode (Node_ptr parent) {
    Node_ptr p = std::make_shared<Node>(_nextNodeID++, parent);
    p->data.firstAppearance = _step;
    p->data.lastAppearance = _step;
    p->data.count = 1;

    _nodes.push_back(p);

    if (parent) parent->children.push_back(p);
    if (_callbacks)  _callbacks->onNewSpecies(parent ? parent->id : NoID, p->id);

    return p;
  }

  SID addGenome (const GENOME &g, Node_ptr species) {
    if (PTreeConfig::DEBUG())
      std::cerr << "Adding genome " << g.id() << " to species " << species->id << std::endl;

    DCCache dccache;
    // Compatible enough with current species
    if (matchesSpecies(g, species, dccache)) {
      insertInto(_step, g, species, dccache, _callbacks);
      _idToSpecies.insert(g.id(), species->id);
      return species->id;
    }

    if (PTreeConfig::DEBUG())
      std::cerr << "\tIncompatible with " << species->id << std::endl;

    // Belongs to subspecies ?
    for (Node_ptr &subspecies: species->children) {
      if (matchesSpecies(g, subspecies, dccache)) {
        insertInto(_step, g, subspecies, dccache, _callbacks);
        _idToSpecies.insert(g.id(), subspecies->id);
        return subspecies->id;
      }
    }

    // Need to create new species
    if (PTreeConfig::simpleNewSpecies()) {
      Node_ptr subspecies = makeNode(species);

      insertInto(_step, g, subspecies, dccache, _callbacks);
      _idToSpecies.insert(g.id(), subspecies->id);
      return subspecies->id;

    } else {
      assert(false);
    }

    assert(false);
    return NoID;
  }

  friend bool matchesSpecies (const GENOME &g, Node_ptr species, DCCache &dccache) {
    uint k = species->enveloppe.size();
    dccache.reserve(k);

    uint matable = 0;
    for (const GENOME &e: species->enveloppe) {
      double d = distance(g, e);
      double c = std::min(g.const_cdata()(d), e.const_cdata()(d));

      if (c >= PTreeConfig::compatibilityThreshold()) matable++;
      dccache.push_back(d, c);
    }

    return matable >= PTreeConfig::similarityThreshold() * k;
  }

  friend void insertInto (uint step, const GENOME &g, Node_ptr species,
                          const DCCache &dccache, Callbacks *callbacks) {

    const uint k = species->enveloppe.size();

    if (PTreeConfig::DEBUG())
      std::cerr << "\tCompatible with " << species->id << std::endl;

    auto &dist = species->distances;

    // Populate the enveloppe
    if (species->enveloppe.size() < PTreeConfig::enveloppeSize()) {
      if (PTreeConfig::DEBUG())  std::cerr << "\tAppend to the enveloppe" << std::endl;

      species->enveloppe.push_back(g);
      if (callbacks)  callbacks->onGenomeEntersEnveloppe(species->id, g.id());
      for (uint i=0; i<k; i++)
        dist[{i, k}] = dccache.distances[i];

    // Better enveloppe point ?
    } else {
      assert(k == PTreeConfig::enveloppeSize());

      /// Find most similar current enveloppe point
      double minDistance = dccache.distances[0];
      uint closest = 0;
      for (uint i=1; i<k; i++) {
        double d = dccache.distances[i];
        if (d < minDistance) {
          minDistance = d;
          closest = i;
        }
      }
      if (PTreeConfig::DEBUG() >= 2)
        std::cerr << "\t\tClosest to "
                  << closest
                  << " (id: " << species->enveloppe[closest].id()
                  << ", d = " << minDistance << ")" << std::endl;

      /// Compute number of times 'g' is better (i.e. more distinct) than the one on the ejectable seat
      uint newIsBest = 0;
      for (uint i=0; i<k; i++) {
        if (i != closest) {
          if (PTreeConfig::DEBUG() >= 2)
            std::cerr << "\t\t" << i << "(" << species->enveloppe[i].id()
                      << "): " << dccache.distances[i] << " >? "
                      << dist[{i,closest}] << std::endl;

          newIsBest += (dccache.distances[i] > dist[{i,closest}]);
        }
      }

      // Genome inside the enveloppe. Nothing to do
      if (newIsBest < PTreeConfig::outperformanceThreshold() * (k-1)) {
        if (PTreeConfig::DEBUG())
          std::cerr << "\tGenome deemed unremarkable with "
                    << k - 1 - newIsBest << " to " << newIsBest << std::endl;

      // Replace closest enveloppe point with new one
      } else {
        if (PTreeConfig::DEBUG())
          std::cerr << "\tReplaced enveloppe point " << closest
                    << " with a vote of " << newIsBest << " to " << k - 1 - newIsBest
                    << std::endl;

        if (callbacks) {
          callbacks->onGenomeLeavesEnveloppe(species->id, species->enveloppe[closest].id());
          callbacks->onGenomeEntersEnveloppe(species->id, g.id());
        }
        species->enveloppe[closest] = g;
        for (uint i=0; i<k; i++)
          if (i != closest)
            dist[{i,closest}] = dccache.distances[i];
      }
    }

    species->data.count++;
    species->data.lastAppearance = step;
  }

// =============================================================================
// == Json conversion

private:
  Node_ptr rebuildHierarchy(Node_ptr parent, const nlohmann::json &j) {
    Node_ptr n = makeNode(parent);

    uint i=0;
    n->id = j[i++];
    n->data = j[i++];
    n->enveloppe = j[i++].get<decltype(Node::enveloppe)>();
    const nlohmann::json jd = j[i++];
    const nlohmann::json jc = j[i++];

    for (const auto &d: jd)
      n->distances[{d[0], d[1]}] = d[2];

    for (const auto &c: jc)
      rebuildHierarchy(n, c);

    return n;
  }

public:
  friend void to_json (nlohmann::json &j, const PhylogenicTree &pt) {
    j = {pt._step, *pt._root};
  }

  friend void to_json (nlohmann::json &j, const Node &n) {
    nlohmann::json jd;
    for (const auto &d: n.distances) jd.push_back({d.first.first, d.first.second, d.second});

    nlohmann::json jc;
    for (const auto &c: n.children) jc.push_back(*c);

    j = {n.id, n.data, n.enveloppe, jd, jc};
  }

  friend void from_json (const nlohmann::json &j, PhylogenicTree &pt) {
    uint i=0;
    pt._step = j[i++];
    pt._root = pt.rebuildHierarchy(nullptr, j[i++]);
  }
};

} // end of namespace phylogeny

#endif // _PHYLOGENIC_TREE_H_
