#ifndef _PHYLOGENIC_TREE_H_
#define _PHYLOGENIC_TREE_H_

#include <vector>
#include <map>
#include <memory>
#include <fstream>

#include <cassert>
#include <iostream>

#include "../ptreeconfig.h"

#include "treetypes.h"
#include "node.hpp"
#include "callbacks.hpp"

/*!
 * \file phylogenictree.hpp
 *
 * Contains the core classes for the phylogeny algorithms
 */

namespace phylogeny {

/// When kept informed about the birth/death and stepping events of a simulation
/// this struct can generate a valid, complete, record of all species event
/// with informations on both the hierarchical and invididual dynamics
///
/// \tparam GENOME the genome of the observed individuals.
template <typename GENOME>
class PhylogenicTree {
  /// Helper lambda for debug printing
  static constexpr auto debug = [] {
    return config::PTree::DEBUG_LEVEL() * config::PTree::DEBUG_PTREE();
  };

public:
  /// Helper alias for the Parent enumeration
  using Parent = genotype::BOCData::Parent;

  /// Helper alias to a species node
  using Node = phylogeny::Node<GENOME>;

  /// \copydoc Node::Ptr
  using Node_ptr = typename Node::Ptr;

  /// \copydoc Node::Collection
  using Nodes = typename Node::Collection;

  /// Specialization used by this tree. Uses CRTP
  using Callbacks = Callbacks_t<PhylogenicTree<GENOME>>;

  /// Helper alias for the configuration data
  using Config = config::PTree;

  /// Helper alias to the type used to cache distance/compatibilities values
  using DCCache = _details::DCCache;

  /// \copydoc Contributors::Contribution
  using SpeciesContribution = typename Contributors::Contribution;

  /// Create an empty PTree
  PhylogenicTree(void) {
    _nextNodeID = SID(0);
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

  /// \return the node with \p sid
  const auto& nodeAt (SID i) const {
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
      nodeAt(sid)->data.lastAppearance = step;
    _step = step;

    // Potentially notify outside world
    if (_callbacks) _callbacks->onStepped(step, aliveSpecies);

    if (debug())
      std::cerr << _idToSpecies.size() << " id>species pairs stored" << std::endl;
  }

  /// Insert \p g into this PTree
  /// \return The species \p g was added to
  SID addGenome (const GENOME &g) {
    // Ensure that the root exists
    if (!_root) {
      _root = makeNode(SpeciesContribution{});
      return updateSpeciesContents(g, _root, DCCache{});
    }

    // Retrieve parent's species
    SID mSID = _idToSpecies.parentSID(g, Parent::MOTHER),
        fSID = _idToSpecies.parentSID(g, Parent::FATHER);

    Node_ptr s0 = nullptr, s1 = nullptr;
    if (mSID == SID::INVALID && fSID == SID::INVALID)
      s0 = _root;

    else if (fSID == SID::INVALID || mSID == fSID)
      s0 = _nodes[mSID];

    else {
      s0 = _nodes[mSID];
      s1 = _nodes[fSID];
    }

    SID sid = addGenome(g, s0, s1, {mSID, fSID});

    if (Config::DEBUG_LEVEL())  std::cerr << std::endl;
    return sid;
  }

  /// Remove \p g from this PTree (and update relevant internal data)
  SID delGenome (const GENOME &g) {
    auto sid = _idToSpecies.remove(g);
    if (debug())
      std::cerr << "New last appearance of species " << sid << " is " << _step << std::endl;

    nodeAt(sid)->data.lastAppearance = _step;
    return sid;
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

private:
  /// Identificator for the next species
  SID _nextNodeID;

protected:
  /// Wraps incrementation of the species identificator counter
  SID nextNodeID (void) {
    using SID_t = std::underlying_type<SID>::type;
    SID curr = _nextNodeID;
    _nextNodeID = SID(SID_t(_nextNodeID)+1);
    return curr;
  }

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
    /// \return The species identificator of \p g 's parent \p p or SID::INVALID
    /// if \p g does not have such a parent
    SID parentSID (const GENOME &g, Parent p) {
      if (!g.hasParent(p))  return SID::INVALID;

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
  Nodes _nodes;

  /// Genome to species lookup table
  IdToSpeciesMap _idToSpecies;

  /// Pointer to the callbacks object. Null by default
  mutable Callbacks *_callbacks;

  uint _hybrids;  ///< Number of hybrids encountered thus far \todo WIP
  uint _step; ///< Current timestep for this tree

  /// Create a smart pointer to a node created on-the-fly with contributors
  /// as described in \p initialContrib
  /// Callbacks:
  ///   - Callbacks_t::onNewSpecies
  Node_ptr makeNode (Contributors &&initialContrib,
                     bool withCallbacks = true) {

    Node_ptr p = Node::make_shared(initialContrib, _nodes);
    p->data.firstAppearance = _step;
    p->data.lastAppearance = _step;
    p->data.count = 1;

    assert(p->contributors.getNodeID()
           == SID(std::underlying_type<SID>::type(_nextNodeID)-1));

    _nodes.push_back(p);

    Node *parent = p->parent();
    if (parent) parent->addChild(p);
    if (_callbacks && withCallbacks)
      _callbacks->onNewSpecies(parent ? parent->id() : SID::INVALID, p->id());

    return p;
  }

  /// \overload
  /// The species contribution monitor is build on the fly with the provided
  /// contributions
  Node_ptr makeNode (const SpeciesContribution &contrib) {
    SID id = nextNodeID();
    Contributors c (id);
    c.update(contrib, Node::elligibilityTester(_nodes));
    return makeNode(std::move(c));
  }

  /// \copydoc nodeAt
  auto& nodeAt (SID i) {
    return const_cast<Node_ptr&>(std::as_const(*this).nodeAt(i));
  }

  /// \return Whether \p g is similar enough to \p species
  float speciesMatchingScore (const GENOME &g, Node_ptr species,
                              DCCache &dccache) {
    uint k = species->enveloppe.size();

    dccache.clear();
    dccache.reserve(k);

    uint matable = 0;
    for (const GENOME &e: species->enveloppe) {
      double d = distance(g, e);
      double c = std::min(g.const_cdata()(d), e.const_cdata()(d));

      if (c >= Config::compatibilityThreshold()) matable++;
      dccache.push_back(d, c);
    }

    assert(dccache.size() == k);
    return matable - Config::similarityThreshold() * k;
  }

  /// Find the appropriate place for \p g in the subtree(s) rooted at
  ///  \p species0 (and species1)
  SID addGenome (const GENOME &g, Node_ptr species0, Node_ptr species1,
                 const SpeciesContribution &contrib) {

    if (debug()) {
      std::cerr << "Attempting to add genome " << g.id()
                << " to species ";
      if (!species1)
        std::cerr << species0->id();
      else
        std::cerr << "either " << species0->id() << " or " << species1->id();
      std::cerr << std::endl;
    }

    DCCache dccache, bestSpeciesDCCache;
    Node_ptr bestSpecies = nullptr;
    float bestScore = -std::numeric_limits<float>::max();

    std::vector<Node_ptr> species;
    species.push_back(species0);
    if (species1) species.push_back(species0);

    // Find best top-level species
    for (Node_ptr s: species) {
      float score = speciesMatchingScore(g, s, dccache);
      if (bestScore < score) {
        bestSpecies = s;
        bestScore = score;
        bestSpeciesDCCache = dccache;
      }
    }

    // Compatible enough with current species ?
    if (bestScore > 0)
      return updateSpeciesContents(g, bestSpecies, bestSpeciesDCCache, contrib);

    if (debug()) {
      std::cerr << "\tIncompatible with ";
      if (!species1)  std::cerr << species0->id();
      else  std::cerr << "both " << species0->id() << " and " << species1->id();
      std::cerr << " (score=" << bestScore << ")" << std::endl;
    }

    // Find best derived species
    uint i=0, k = 0;
    size_t max = species0->children().size();
    if (species1) max = std::max(max, species1->children().size());
    while (bestScore <= 0 && i < max) {
      Node_ptr &subspecies = species[k]->child(i);
      float score = speciesMatchingScore(g, subspecies, dccache);
      if (bestScore < score) {
        bestSpecies = subspecies;
        bestScore = score;
        bestSpeciesDCCache = dccache;
      }

      if (!species1 || species[1-k]->children().size() <= i)
        i++;
      else {
        if (k==1) i++;
        k = 1 - k;
      }
    }

    // Belongs to subspecies ?
    if (bestScore > 0) {
      if (debug())
        std::cerr << "Compatible with " << bestSpecies->id()
                  << " (score=" << bestScore << ")" << std::endl;
      return updateSpeciesContents(g, bestSpecies, bestSpeciesDCCache, contrib);
    }

    // Need to create new species
    if (Config::simpleNewSpecies()) {
      Node_ptr subspecies = makeNode(contrib);
      dccache.clear();
      if (debug())
        std::cerr << "Created new species " << subspecies->id() << std::endl;
      return updateSpeciesContents(g, subspecies, dccache);

    } else
      assert(false);

    assert(false);
    return SID::INVALID;
  }

  /// Insert \p g into node \p species, possibly changing the enveloppe.
  ///
  /// Callbacks:
  ///   - Callbacks_t::onGenomeEntersEnveloppe
  ///   - Callbacks_t::onGenomeLeavesEnveloppe
  friend void insertInto (uint step, const GENOME &g, Node_ptr species,
                          const DCCache &dccache, Callbacks *callbacks) {

    const uint k = species->enveloppe.size();

    auto &dist = species->distances;

    // Populate the enveloppe
    if (k < Config::enveloppeSize()) {
      if (debug())  std::cerr << "\tAppend to the enveloppe" << std::endl;

      species->enveloppe.push_back(g);
      if (callbacks)  callbacks->onGenomeEntersEnveloppe(species->id(), g.id());
      for (uint i=0; i<k; i++)
        dist.at({i, k}) = dccache.distances[i];

    // Better enveloppe point ?
    } else {
      assert(k == Config::enveloppeSize());
      std::vector<GID> ids (k);
      for (uint i=0; i<k; i++)  ids[i] = species->enveloppe[i].id();
      _details::EnveloppeContribution ec =
          computeContribution(dist, dccache.distances, g.id(), ids);

      // Genome inside the enveloppe. Nothing to do
      if (!ec.better) {
        if (debug())
          std::cerr << "\t" << g.id() << "'s contribution is too low ("
                    << ec.value << ")" << std::endl;

      // Replace closest enveloppe point with new one
      } else {
        if (debug())
          std::cerr << "\t" << g.id() << "'s contribution is better "
                    << "than enveloppe point " << ec.than << " (id: "
                    << species->enveloppe[ec.than].id()
                    << ", c = " << ec.value << ")" << std::endl;

        if (callbacks) {
          callbacks->onGenomeLeavesEnveloppe(species->id(), species->enveloppe[ec.than].id());
          callbacks->onGenomeEntersEnveloppe(species->id(), g.id());
        }
        species->enveloppe[ec.than] = g;
        for (uint i=0; i<k; i++)
          if (i != ec.than)
            dist.at({i,ec.than}) = dccache.distances[i];
      }
    }

    species->data.count++;
    species->data.lastAppearance = step;
  }

  /// Update species \p s by inserting genome \p g, updating the contributions
  /// and registering the GID>SID association
  SID updateSpeciesContents(const GENOME &g, Node_ptr s,
                            const DCCache &cache,
                            const SpeciesContribution &ctb = {}) {

    insertInto(_step, g, s, cache, _callbacks);
    if (!ctb.empty()) updateContributions(s, ctb);
    _idToSpecies.insert(g.id(), s->id());
    return s->id();
  }

  /// Update species \p s contributions with the provided values
  void updateContributions (Node_ptr s, const SpeciesContribution &contrib) {
    Node *oldMC = s->parent(),
         *newMC = s->update(contrib, _nodes);

    if (oldMC != newMC) {
      /// \todo remove
      if (!oldMC || !newMC)
        s->update({}, _nodes);

      assert(oldMC);
      assert(newMC);

      // Parent changed. Update and notify
      oldMC->delChild(s);
      newMC->addChild(s);
      _root->updateElligibilities(_nodes);

      _callbacks->onMajorContributorChanged(s->id(), oldMC->id(), newMC->id());
    }
  }

// =============================================================================
// == Json conversion

private:
  /// Rebuilds PTree hierarchy and internal structure based on the contents of
  /// json \p j
  Node_ptr rebuildHierarchy(const json &j) {
    Contributors c = j["contribs"];
    Node_ptr n = makeNode(std::move(c), false);

    n->data = j["data"];
    n->enveloppe = j["envlp"].get<decltype(Node::enveloppe)>();
    const json jd = j["dists"];
    const json jc = j["children"];

    for (const auto &d: jd)
      n->distances.at({d[0], d[1]}) = d[2];

    for (const auto &c: jc)
      rebuildHierarchy(c);

    return n;
  }

public:
  /// Serialise PTree \p pt into a json
  friend void to_json (json &j, const PhylogenicTree &pt) {
    j = {pt._step, *pt._root};
  }

  /// Serialize Node \p n into a json
  friend void to_json (json &j, const Node &n) {
    json jd;
    for (const auto &d: n.distances) jd.push_back({d.first.first, d.first.second, d.second});

    json jc;
    for (const auto &c: n.children()) jc.push_back(*c);

    j["id"] = n.id();
    j["data"] = n.data;
    j["envlp"] = n.enveloppe;
    j["contribs"] = n.contributors;
    j["dists"] = jd;
    j["children"] = jc;
  }

  /// Deserialise PTree \p pt from json \p j
  friend void from_json (const json &j, PhylogenicTree &pt) {
    uint i=0;
    pt._step = j[i++];
    pt._root = pt.rebuildHierarchy(j[i++]);
  }
};

} // end of namespace phylogeny

#endif // _PHYLOGENIC_TREE_H_
