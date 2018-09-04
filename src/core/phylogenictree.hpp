#ifndef _PHYLOGENIC_TREE_H_
#define _PHYLOGENIC_TREE_H_

#include <vector>
#include <map>
#include <memory>
#include <fstream>

#include <cassert>
#include <iostream>

#include "kgd/external/json.hpp"
#include "ptreeconfig.h"
#include "crossover.hpp"

/*!
 * \file phylogenictree.hpp
 *
 * Contains the core classes for the phylogeny algorithms
 */

namespace phylogeny {

/// \cond internal

/// Contains a set of functions called by the phylogenic tree when each of the
/// corresponding events occur
template <typename PT>
struct Callbacks_t {
  /// Alias to the genetic ID
  using GID = typename PT::GID;

  /// Alias to the species ID
  using SID = typename PT::SID;

  /// Alias to the collection of still-alive species id
  using LivingSet = typename PT::LivingSet;

  /// \brief Called when the PTree has been stepped.
  ///
  /// Provides the current step and the set of still-alive species
  void onStepped (uint step, const LivingSet &living);

  /// \brief Called to notify of a newly created species.
  ///
  /// Provides the identificators of both parent (if any) and new species
  void onNewSpecies (SID pid, SID sid);

  /// \brief Called when a genome is added to an enveloppe.
  ///
  /// Provides the id of the species whose enveloppe just changed and the id
  /// of the newly inserted genome
  void onGenomeEntersEnveloppe (SID sid, GID gid);

  /// \brief Called when a genome is removed from an enveloppe.
  ///
  /// Provides the id of the species whose enveloppe just changed and the id
  /// of the newly removed genome
  void onGenomeLeavesEnveloppe (SID sid, GID gid);
};

/// \endcond

/// Contains data regarding a specific species
struct SpeciesData {

  /// Time at which this species was first observed
  uint firstAppearance;

  /// Time at which this species was last observed
  uint lastAppearance;

  /// Number of individuals of the species observed
  uint count;

  /// Serialize to json
  friend void to_json (nlohmann::json &j, const SpeciesData &d) {
    j = {d.firstAppearance, d.lastAppearance, d.count};
  }

  /// Deserialize from json
  friend void from_json (const nlohmann::json &j, SpeciesData &d) {
    uint i=0;
    d.firstAppearance = j[i++];
    d.lastAppearance = j[i++];
    d.count = j[i++];
  }
};

/// When kept informed about the birth/death and stepping events of a simulation
/// this struct can generate a valid, complete, record of all species event
/// with informations on both the hierarchical and invididual dynamics
///
/// \tparam GENOME the genome of the observed individuals.
template <typename GENOME>
class PhylogenicTree {
public:

  /// Helper alias for the genome identificator
  using GID = genotype::BOCData::GID;

  /// Alias for the species identificator
  using SID = uint;

  /// Value indicating an unspecified species
  static constexpr SID NoID = SID(-1);

  /// Helper alias for the Parent enumeration
  using Parent = genotype::BOCData::Parent;

  /// Collections of still-alive species identificators
  using LivingSet = std::set<SID>;

  /// Specialization used by this tree. Uses CRTP
  using Callbacks = Callbacks_t<PhylogenicTree<GENOME>>;

  /// Helper alias for the configuration data
  using Config = config::PTree;

  /// Create an empty PTree
  PhylogenicTree(void) {
    _nextNodeID = 0;
    _step = 0;
    _hybrids = 0;
    _root = nullptr;
    _callbacks = nullptr;
  }

  /// Sets the callbacks used by this ptree
  void setCallbacks (Callbacks *c) const { _callbacks = c; }

  /// \return the callbacks used by this ptree
  Callbacks* callbacks (void) {   return _callbacks; }

  /// \return a smart pointer to the root (can be null)
  const auto& root (void) const {
    return _root;
  }

  /// \return the number of nodes in this tree
  uint width (void) const {
    return _nodes.size();
  }

  /// \return the node with id/index i
  const auto& nodeAt (uint i) const {
    return _nodes.at(i);
  }

  /// \return the current timestep for this PTree
  uint step (void) const {
    return _step;
  }

  /// Sets the current timestep for this PTree
  void setStep (uint step) {
    _step = step;
  }

  /// Update the set of still-alive species based on the list of still-alive
  /// genomes [\p begin,\p end[ extracted through \p geneticID
  ///
  /// Callbacks:
  ///   - Callbacks_t::onStepped
  ///
  /// \tparam IT Iterator to the begin/end of the population list
  /// \tparam F Functor for extracting the genome id from an iterator
  template <typename IT, typename F>
  void step (uint step, IT begin, IT end, F geneticID) {

    // Determine which species are still alive
    LivingSet aliveSpecies;
    for (IT it = begin; it != end; ++it)
      aliveSpecies.insert(_idToSpecies.at(geneticID(*it)));

    // Update internal data
    for (SID sid: aliveSpecies)
      _nodes.at(sid)->data.lastAppearance = step;
    _step = step;

    // Potentially notify outside world
    if (_callbacks) _callbacks->onStepped(step, aliveSpecies);

    if (Config::DEBUG())
      std::cerr << _idToSpecies.size() << " id>species pairs stored" << std::endl;
  }

  /// Insert \p g into this PTree
  /// \return The species \p g was added to
  SID addGenome (const GENOME &g) {
    // Ensure that the root exists
    if (!_root) _root = makeNode(nullptr);

    // If genome is parent-less, assume that it comes from initialization data
    if (!g.hasParent(Parent::FATHER) || !g.hasParent(Parent::MOTHER))
      return addGenome(g, _root);

    // Retrieve parent's species
    uint mSID = _idToSpecies.parentSID(g, Parent::MOTHER),
         pSID = _idToSpecies.parentSID(g, Parent::FATHER);

    // Manage (poorly) hydrisation
    assert(Config::ignoreHybrids() || mSID == pSID);
    if (mSID != pSID) _hybrids++;

    // All is well go ahead
    if (mSID == pSID)
      return addGenome(g, _nodes[mSID]);

    // Got a hybrid -> use mother species instead
    else if (Config::ignoreHybrids()) {
      if (Config::DEBUG() >= 0)
        std::cerr << "Linking hybrid genome " << g.id() << " to mother species" << std::endl;

      return addGenome(g, _nodes[mSID]);

    /// \todo Find a way to deal with this
    } else {
      assert(false);
      if (Config::DEBUG())
        std::cerr << "Managing hybrid genome " << g.id() << std::endl;
    }

    // Should never reach this point
    assert(false);
    return NoID;
  }

  /// Remove \p g from this PTree (and update relevant internal data)
  void delGenome (const GENOME &g) {
    auto sid = _idToSpecies.remove(g);
    if (Config::DEBUG())
      std::cerr << "New last appearance of species " << sid << " is " << _step << std::endl;

    _nodes[sid]->data.lastAppearance = _step;
  }

  /// Stream \p pt to \p os. Mostly for debugging purpose: output is quickly
  /// unintelligible
  friend std::ostream& operator<< (std::ostream &os, const PhylogenicTree &pt) {
    os << pt._hybrids << " Hybrids;\n";
    return os << *pt._root;
  }

  /// Dump this PTree into dot file \p filename
  void logTo (const std::string &filename) const {
    std::ofstream ofs (filename);
    ofs << "digraph {\n";
    _root->logTo(ofs);
    ofs << "}\n";
  }

protected:
  /// Identificator for the next species
  SID _nextNodeID;

  /// Helper structure for ensuring that the pair values are ordered
  ///
  /// \tparam T stored type
  template <typename T>
  struct ordered_pair {
    T first;  ///< First value (lower or equal to second)
    T second; ///< Second value (greater or equal to first)

    /// Create an ordered pair from un-ordered pair of values
    ordered_pair(T first, T second)
      : first(std::min(first, second)), second(std::max(first, second)) {}

    /// Compare two ordered_pairs in lexicographic order
    friend bool operator< (const ordered_pair &lhs, const ordered_pair &rhs) {
      if (lhs.first != rhs.first) return lhs.first < rhs.first;
      return lhs.second < rhs.second;
    }
  };

  struct Node;

  /// Smart pointer (shared) to a species node
  using Node_ptr = std::shared_ptr<Node>;

  /// Species node
  struct Node {
    SID id; ///< Species identificator
    SpeciesData data; ///< Species additionnal data

    Node *parent; ///< Species parent (if any)

    /// Subspecies of this node
    std::vector<Node_ptr> children;

    /// Collection of borderoids (in opposition to centroids)
    std::vector<GENOME> enveloppe;

    /// Cache map for the intra-enveloppe distances
    std::map<ordered_pair<uint>, float> distances;

    /// Create a node with \p id and \p parent
    Node (SID id, Node_ptr parent) : id(id), parent(parent.get()) {}

    /// Stream this node. Mostly for debugging purpose.
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

    /// Dump this node, in dot format.
    void logTo (std::ostream &os) const {
      os << "\t" << id << ";\n";
      for (const Node_ptr &n: children) {
        os << "\t" << id << " -> " << n->id << ";\n";
        n->logTo(os);
      }
    }
  };

  /// Distance & compatibilities cache
  struct DCCache {
    /// Cache collection of distances
    std::vector<float> distances;

    /// Cache collection of compatibilities
    std::vector<float> compatibilities;

    /// Prepare exactly \p n units of storage space
    void reserve (uint n) {
      distances.reserve(n), compatibilities.reserve(n);
    }

    /// Append values
    void push_back (float d, float c) {
      distances.push_back(d), compatibilities.push_back(c);
    }
  };

  /// Allows mapping a genome ID to its species ID
  struct IdToSpeciesMap {

    /// Efficient data storage that retains the number of uses of a genetic ID
    /// whether by the organism itself or by its direct descendants
    struct ITSMData {
      SID species;    ///< The species
      uint refCount;  ///< The current number of users of this data
    };

    /// Internal container
    std::map<GID, ITSMData> map;

    /// \return the size of the lookup table
    uint size (void) const {
      return map.size();
    }

    /// Queries for the parent species identificator for a genome.
    /// Updates internal reference counter.
    /// \return The species identificator of \p g 's parent \p p
    /// \warning \p g must have a valid parent \p p
    SID parentSID (const GENOME &g, Parent p) {
      auto &d = map.at(g.parent(p));
      d.refCount++;
      return d.species;
    }

    /// Decrement the reference counter for key \p id and removes it upon reaching 0
    SID remove (const GID id) {
      auto it = map.find(id);
      auto sid = it->second.species;
      auto &ref = it->second.refCount;

      assert(ref > 0);
      ref--;

      if (ref == 0)  map.erase(it);

      return sid;
    }

    /// Try to remove \p g (and its eventual parents) from the association map
    SID remove (const GENOME &g) {
      auto sid = remove(g.id());
      for (Parent p: {Parent::MOTHER, Parent::FATHER})
        if (g.hasParent(p))
          remove(g.parent(p));
      return sid;
    }

    /// Insert data for a genome \p gid belonging to species \p sid
    void insert (GID gid, SID sid) {
      ITSMData d;
      d.refCount = 1;
      d.species = sid;
      map[gid] = d;
    }

    /// \return species indentificator for genome \p gid
    /// \throws std::out_of_range if gid &notin; map
    SID at (GID gid) const {
      return map.at(gid).species;
    }
  };

  /// The PTree root. Null until the first genome is inserted
  Node_ptr _root;

  /// Nodes linear collection for direct access
  std::vector<Node_ptr> _nodes;

  /// Genome to species lookup table
  IdToSpeciesMap _idToSpecies;

  /// Pointer to the callbacks object. Null by default
  mutable Callbacks *_callbacks;

  uint _hybrids;  ///< Number of hybrids encountered thus far (\todo DEBUGGING)
  uint _step; ///< Current timestep for this tree

  /// Create a smart pointer to a node created on-the-fly under \p parent
  /// Callbacks:
  ///   - Callbacks_t::onNewSpecies

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

  /// Find the appropriate place for \p g in the subtree rooted at \p species
  SID addGenome (const GENOME &g, Node_ptr species) {
    if (Config::DEBUG())
      std::cerr << "Adding genome " << g.id() << " to species " << species->id << std::endl;

    DCCache dccache;
    // Compatible enough with current species
    if (matchesSpecies(g, species, dccache)) {
      insertInto(_step, g, species, dccache, _callbacks);
      _idToSpecies.insert(g.id(), species->id);
      return species->id;
    }

    if (Config::DEBUG())
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
    if (Config::simpleNewSpecies()) {
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

  /// \return Whether \p g is similar enough to \p species
  friend bool matchesSpecies (const GENOME &g, Node_ptr species, DCCache &dccache) {
    uint k = species->enveloppe.size();
    dccache.reserve(k);

    uint matable = 0;
    for (const GENOME &e: species->enveloppe) {
      double d = distance(g, e);
      double c = std::min(g.const_cdata()(d), e.const_cdata()(d));

      if (c >= Config::compatibilityThreshold()) matable++;
      dccache.push_back(d, c);
    }

    return matable >= Config::similarityThreshold() * k;
  }

  /// Insert \p g into node \p species, possibly changing the enveloppe.
  ///
  /// Callbacks:
  ///   - Callbacks_t::onGenomeEntersEnveloppe
  ///   - Callbacks_t::onGenomeLeavesEnveloppe
  friend void insertInto (uint step, const GENOME &g, Node_ptr species,
                          const DCCache &dccache, Callbacks *callbacks) {

    const uint k = species->enveloppe.size();

    if (Config::DEBUG())
      std::cerr << "\tCompatible with " << species->id << std::endl;

    auto &dist = species->distances;

    // Populate the enveloppe
    if (species->enveloppe.size() < Config::enveloppeSize()) {
      if (Config::DEBUG())  std::cerr << "\tAppend to the enveloppe" << std::endl;

      species->enveloppe.push_back(g);
      if (callbacks)  callbacks->onGenomeEntersEnveloppe(species->id, g.id());
      for (uint i=0; i<k; i++)
        dist[{i, k}] = dccache.distances[i];

    // Better enveloppe point ?
    } else {
      assert(k == Config::enveloppeSize());

      // Find most similar current enveloppe point
      double minDistance = dccache.distances[0];
      uint closest = 0;
      for (uint i=1; i<k; i++) {
        double d = dccache.distances[i];
        if (d < minDistance) {
          minDistance = d;
          closest = i;
        }
      }
      if (Config::DEBUG() >= 2)
        std::cerr << "\t\tClosest to "
                  << closest
                  << " (id: " << species->enveloppe[closest].id()
                  << ", d = " << minDistance << ")" << std::endl;

      // Compute number of times 'g' is better (i.e. more distinct) than the one on the ejectable seat
      uint newIsBest = 0;
      for (uint i=0; i<k; i++) {
        if (i != closest) {
          if (Config::DEBUG() >= 2)
            std::cerr << "\t\t" << i << "(" << species->enveloppe[i].id()
                      << "): " << dccache.distances[i] << " >? "
                      << dist[{i,closest}] << std::endl;

          newIsBest += (dccache.distances[i] > dist[{i,closest}]);
        }
      }

      // Genome inside the enveloppe. Nothing to do
      if (newIsBest < Config::outperformanceThreshold() * (k-1)) {
        if (Config::DEBUG())
          std::cerr << "\tGenome deemed unremarkable with "
                    << k - 1 - newIsBest << " to " << newIsBest << std::endl;

      // Replace closest enveloppe point with new one
      } else {
        if (Config::DEBUG())
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
  /// Rebuilds PTree hierarchy and internal structure based on the contents of
  /// json \p j
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
  /// Serialise PTree \p pt into a json
  friend void to_json (nlohmann::json &j, const PhylogenicTree &pt) {
    j = {pt._step, *pt._root};
  }

  /// Serialize Node \p n into a json
  friend void to_json (nlohmann::json &j, const Node &n) {
    nlohmann::json jd;
    for (const auto &d: n.distances) jd.push_back({d.first.first, d.first.second, d.second});

    nlohmann::json jc;
    for (const auto &c: n.children) jc.push_back(*c);

    j = {n.id, n.data, n.enveloppe, jd, jc};
  }

  /// Deserialise PTree \p pt from json \p j
  friend void from_json (const nlohmann::json &j, PhylogenicTree &pt) {
    uint i=0;
    pt._step = j[i++];
    pt._root = pt.rebuildHierarchy(nullptr, j[i++]);
  }
};

} // end of namespace phylogeny

#endif // _PHYLOGENIC_TREE_H_
