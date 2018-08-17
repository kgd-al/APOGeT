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

namespace genotype {

static constexpr int D_PTREE = 0;

template <typename GENOME>
class PhylogenicTree {
  using NodeID = uint;
  static constexpr NodeID NoID = NodeID(-1);

public:

  PhylogenicTree(void) {
    _nextNodeID = 0;
    _step = 0;
    _hybrids = 0;
    _root = makeNode(nullptr);
  }

  void step (uint step, const std::set<uint> &alivePlants) {
    std::set<NodeID> aliveSpecies;
    for (uint pid: alivePlants)
      aliveSpecies.insert(_idToSpecies.at(pid));

    for (NodeID sid: aliveSpecies)
      _nodes.at(sid)->data.lastAppearance = step;

    _step = step;
  }

  NodeID addGenome (int x, const GENOME &g) {
    if (!g.hasParent(1) || !g.hasParent(0))
      return addGenome(x, g, _root);

    uint mID = _idToSpecies.at(g.parent(0)),
         pID = _idToSpecies.at(g.parent(1));

    assert(PTreeConfig::ignoreHybrids() || mID == pID);
    if (mID != pID) _hybrids++;
    if (mID == pID)
      return addGenome(x, g, _nodes[mID]);

    else if (PTreeConfig::ignoreHybrids()) {
      if (D_PTREE)  std::cerr << "Linking hybrid genome " << g.id() << " to mother species" << std::endl;
      return addGenome(x, g, _nodes[mID]);

    } else {
      assert(false);
      if (D_PTREE)  std::cerr << "Managing hybrid genome " << g.id() << std::endl;
    }

    return NoID;
  }

  void delGenome (uint step, uint id) {
    auto it = _idToSpecies.find(id);
    if (it != _idToSpecies.end()) {
      if (D_PTREE)  std::cerr << "New last appearance of species " << it->second << " is " << step << std::endl;
      _nodes[it->second]->data.lastAppearance = step;
    }
//    _idToSpecies.erase(id);
  }

  NodeID speciesID (uint id) {
    auto it = _idToSpecies.find(id);
    if (it != _idToSpecies.end())
          return it->second;
    else  return NoID;
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
  template <typename _GENOME>
  friend class PTreeIntrospecter;

  NodeID _nextNodeID;

  std::map<uint, NodeID> _idToSpecies;

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
    NodeID id;

    struct Data {
      uint firstAppearance;
      uint lastAppearance;
      uint count;
      int xmin, xmax;
    } data;

    friend void to_json (nlohmann::json &j, const Data &d) {
      j = {d.firstAppearance, d.lastAppearance, d.count, d.xmin, d.xmax};
    }

    friend void from_json (const nlohmann::json &j, Data &d) {
      uint i=0;
      d.firstAppearance = j[i++];
      d.lastAppearance = j[i++];
      d.count = j[i++];
      d.xmin = j[i++];
      d.xmax = j[i++];
    }

    Node *parent;
    std::vector<Node_ptr> children;

    std::vector<GENOME> enveloppe;
    std::map<ordered_pair<uint>, float> distances;

    Node (NodeID id, Node_ptr parent) : id(id), parent(parent.get()) {}

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

  Node_ptr _root;

  std::vector<Node_ptr> _nodes;

  uint _hybrids;
  uint _step;

  Node_ptr makeNode (Node_ptr parent) {
    Node_ptr p = std::make_shared<Node>(_nextNodeID++, parent);
    p->data.firstAppearance = 0;
    p->data.count = 0;
    _nodes.push_back(p);
    return p;
  }

  NodeID addGenome (int x, const GENOME &g, Node_ptr species) {
    if (D_PTREE)  std::cerr << "Adding genome " << g.id() << " to species " << species->id << std::endl;

    std::vector<float> compatibilities;
    // Compatible enough with current species
    if (matchesSpecies(g, species, compatibilities)) {
      insertInto(_step,x, g, species, compatibilities);
      _idToSpecies[g.id()] = species->id;
      return species->id;
    }

    if (D_PTREE)  std::cerr << "\tIncompatible with " << species->id << std::endl;

    // Belongs to subspecies ?
    for (Node_ptr &subspecies: species->children) {
      if (matchesSpecies(g, subspecies, compatibilities)) {
        insertInto(_step,x, g, subspecies, compatibilities);
        _idToSpecies[g.id()] = subspecies->id;
        return subspecies->id;
      }
    }

    // Need to create new species
    if (PTreeConfig::simpleNewSpecies()) {
      Node_ptr subspecies = makeNode(species);
      subspecies->data.firstAppearance = _step;
      subspecies->data.xmin = x;
      subspecies->data.xmax = x;
      species->children.push_back(subspecies);
      insertInto(_step, x, g, subspecies, compatibilities);
      _idToSpecies[g.id()] = subspecies->id;
      return subspecies->id;

    } else {
      assert(false);
    }

    assert(false);
    return NoID;
  }

  friend bool matchesSpecies (const GENOME &g, Node_ptr species, std::vector<float> &compatibilities) {
    uint k = species->enveloppe.size();
    compatibilities.reserve(k);

    uint matable = 0;
    for (const GENOME &e: species->enveloppe) {
      double d = distance(g, e);
      double c = std::min(g.compatibility(d), e.compatibility(d));

      if (c >= PTreeConfig::compatibilityThreshold()) matable++;
      compatibilities.push_back(c);
    }

    return matable >= PTreeConfig::similarityThreshold() * k;
  }

  friend void insertInto (uint step, int x, const GENOME &g, Node_ptr species,
                          const std::vector<float> &compatibilities) {

    const uint k = species->enveloppe.size();

    if (D_PTREE)  std::cerr << "\tCompatible with " << species->id << std::endl;

    auto &dist = species->distances;

    // Populate the enveloppe
    if (species->enveloppe.size() < PTreeConfig::enveloppeSize()) {
      if (D_PTREE)  std::cerr << "\tAppend to the enveloppe" << std::endl;

      species->enveloppe.push_back(g);
      for (uint i=0; i<k; i++)
        dist[{i, k}] = compatibilities[i];

    // Better enveloppe point ?
    } else {
      assert(k == PTreeConfig::enveloppeSize());

      bool further = false;
      double min = compatibilities[0];
      uint closer = 0;

      for (uint i=0; i<k; i++) {
        for (uint j=0; j<i; j++) {
          double c = compatibilities[j];
          further = (dist[{i,j}] < c);

          if (c < min) {
            min = c;
            closer = i;
          }
        }
      }

      assert(!further || closer < PTreeConfig::enveloppeSize());

      // Quiet genome inside the enveloppe. Nothing to good
      if (!further) {
        if (D_PTREE)  std::cerr << "\tUnremarkable genome" << std::endl;

      // Replace closest enveloppe point with new one
      } else {
        if (D_PTREE)  std::cerr << "\tReplaced enveloppe point " << closer << std::endl;
        species->enveloppe[closer] = g;
        for (uint i=0; i<k; i++)
          if (i != closer)
            dist[{i,closer}] = compatibilities[i];
      }
    }

    species->data.count++;
    species->data.lastAppearance = step;
    species->data.xmin = std::min(species->data.xmin, x);
    species->data.xmax = std::max(species->data.xmax, x);
  }

// =============================================================================
// == Json conversion

private:
  Node_ptr rebuildHierarchy(Node_ptr parent, const nlohmann::json &j) {
    Node_ptr n = makeNode(parent);
    if (parent) parent->children.push_back(n);

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

} // end of namespace genotype

#endif // _PHYLOGENIC_TREE_H_
