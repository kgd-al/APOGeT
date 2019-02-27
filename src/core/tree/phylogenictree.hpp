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

/// Placeholder for end-users not requiering additional phylogenic informations
struct NoUserData {
  NoUserData (void) {}  ///< Default constructor for nlohmann::json
  NoUserData (GID) {}   ///< Discards provided id
  void removedFromEnveloppe (void) const {} ///< Nothing to do
  friend void to_json (json&, const NoUserData&) {}  ///< Nothing to store
  friend void from_json (const json&, NoUserData&) {}  ///< Nothing to read
};

/// When kept informed about the birth/death and stepping events of a simulation
/// this struct can generate a valid, complete, record of all species event
/// with informations on both the hierarchical and invididual dynamics
///
/// \tparam GENOME the genome of the observed individuals.
/// \tparam UDATA user data for collecting sample statistics at the individual
/// level (defaults to nothing)
///
/// FIXME At least once, _idToSpecies.at(foo) throws even though _idToSpecies.contains(foo) == true
template <typename GENOME, typename UDATA>
class PhylogenicTree {
  /// Helper lambda for debug printing
  static constexpr auto debug = [] {
    return config::PTree::DEBUG_LEVEL() * config::PTree::DEBUG_PTREE();
  };

public:
  /// Helper alias for the Parent enumeration
  using Parent = genotype::BOCData::Parent;

  /// Helper alias to a species node
  using Node = phylogeny::Node<GENOME, UDATA>;

  /// \copydoc Node::Ptr
  using Node_ptr = typename Node::Ptr;

  /// \copydoc Node::Collection
  using Nodes = typename Node::Collection;

  /// Specialization used by this tree. Uses CRTP
  using Callbacks = Callbacks_t<PhylogenicTree<GENOME, UDATA>>;

  /// Helper alias for the configuration data
  using Config = config::PTree;

  /// Helper alias to the type used to cache distance/compatibilities values
  using DCCache = _details::DCCache;

  /// \copydoc Contributors::Contributions
  using SpeciesContribution = typename Contributors::Contributions;

  /// Create an empty PTree
  PhylogenicTree(void) {
    _nextNodeID = SID(0);
    _enveloppeSize = Config::enveloppeSize();
    _stillborns = 0;
    _step = 0;
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

  /// \return the user data for enveloppe point \p gid or nullptr if it is a
  /// regular individual
  UDATA* getUserData (GID id) const {
    auto sid = _idToSpecies.at(id);
    const auto &species = _nodes.at(sid);
    for (const auto &ep: species->enveloppe)
      if (ep.genome.crossoverData().id == id)
        return ep.userData.get();
    return nullptr;
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

    static const auto &T = Config::stillbornTrimmingPeriod();
    if ((T > 0) && (_step % T) == 0)  performStillbornTrimming();

    // Potentially notify outside world
    if (_callbacks) _callbacks->onStepped(step, aliveSpecies);

    if (debug())
      std::cerr << _idToSpecies.size() << " id>species pairs stored" << std::endl;
  }

  /// Insert \p g into this PTree
  /// \return The species \p g was added to and, if it was also added to the
  /// enveloppe, a pointer to the associated user data structure
  /// (nullptr otherwise).
  std::pair<SID, UDATA*> addGenome (const GENOME &g) {
    // Ensure that the root exists
    if (!_root) {
      _root = makeNode(SpeciesContribution{});
      return updateSpeciesContents(g, _root, DCCache{}, SpeciesContribution{});
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

    auto ret = addGenome(g, s0, s1, mSID, fSID);

    if (Config::DEBUG_LEVEL())  std::cerr << std::endl;
    return ret;
  }

  /// Remove \p g from this PTree (and update relevant internal data)
  SID delGenome (const GENOME &g) {
    auto sid = _idToSpecies.remove(g);
    if (debug())
      std::cerr << "New last appearance of species " << sid << " is " << _step << std::endl;

    SpeciesData &data = nodeAt(sid)->data;
    data.lastAppearance = _step;
    data.currentlyAlive--;
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

      /// Serializes an item into a json
      friend void to_json (json &j, const ITSMData &d) {
        j = { d.species, d.refCount };
      }

      /// Deserializes an item from a json
      friend void from_json (const json &j, ITSMData &d) {
        d.species = j[0], d.refCount = j[1];
      }
    };

    /// Helper alias to the type of the underlying container
    using map_t = std::map<GID, ITSMData>;

    /// Internal container
    map_t map;

    /// \return the size of the lookup table
    uint size (void) const {
      return map.size();
    }

    /// Queries for the parent species identificator for a genome.
    /// Updates internal reference counter.
    /// \return The species identificator of \p g 's parent \p p or SID::INVALID
    /// if \p g does not have such a parent
    SID parentSID (const GENOME &g, Parent p) {
      if (!g.crossoverData().hasParent(p))  return SID::INVALID;

      auto &d = map.at(g.crossoverData().parent(p));
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

      if (config::PTree::DEBUG_ID2SPECIES())
        std::cerr << "Removed a reference to " << id
                  << " current count is " << ref << std::endl;

      if (ref == 0)  map.erase(it);

      return sid;
    }

    /// Try to remove \p g (and its eventual parents) from the association map
    SID remove (const GENOME &g) {
      auto sid = remove(g.crossoverData().id);
      for (Parent p: {Parent::MOTHER, Parent::FATHER})
        if (g.crossoverData().hasParent(p))
          remove(g.crossoverData().parent(p));
      return sid;
    }

    /// Insert data for a genome \p gid belonging to species \p sid
    void insert (GID gid, SID sid) {
      ITSMData d;
      d.refCount = 1;
      d.species = sid;

      auto res = map.try_emplace(gid, d);
      if (!res.second)
        utils::doThrow<std::logic_error>("Unable to register GID(", gid,
                                         ")>SID(", sid, ") association!");

      if (config::PTree::DEBUG_ID2SPECIES())
        std::cerr << "Inserted a reference to " << gid
                  << " current count is 1" << std::endl;
    }

    /// \return species indentificator for genome \p gid
    /// \throws std::out_of_range if gid &notin; map
    SID at (GID gid) const {
      return map.at(gid).species;
    }

    /// \return Whether or not \p gid is associated to a species
    bool contains (GID gid) const {
      return map.find(gid) != map.end();
    }

    /// Serialize the underlying map into a json
    friend void to_json (json &j, const IdToSpeciesMap &m) {
      j = m.map;
    }

    /// Deserialize the underlying map from a json
    friend void from_json (const json &j, IdToSpeciesMap &m) {
      m.map = j.get<map_t>();
    }
  };

  /// The PTree root. Null until the first genome is inserted
  Node_ptr _root;

  /// Nodes collection for logarithmic access
  Nodes _nodes;

  /// Genome to species lookup table
  IdToSpeciesMap _idToSpecies;

  /// Pointer to the callbacks object. Null by default
  mutable Callbacks *_callbacks;

  uint _enveloppeSize;  ///< Number of enveloppe points
  uint _stillborns; ///< Number of stillborn species removed
  uint _step; ///< Current timestep for this tree

  /// Create a smart pointer to a node created on-the-fly with contributors
  /// as described in \p initialContrib
  /// Callbacks:
  ///   - Callbacks_t::onNewSpecies
  Node_ptr makeNode (const SpeciesContribution &contrib) {

    SID id = nextNodeID();
    Contributors c (id);

    Node_ptr p = Node::make_shared(c);
    p->data.firstAppearance = _step;
    p->data.lastAppearance = _step;
    p->data.count = 0;
    p->data.currentlyAlive = 0;

    assert(p->contributors.getNodeID()
           == SID(std::underlying_type<SID>::type(_nextNodeID)-1));

    _nodes[p->id()] = p;

    // Compute parent
    p->update(contrib, _nodes);

    Node *parent = p->parent();
    if (parent) parent->addChild(p);
    if (_callbacks)
      _callbacks->onNewSpecies(parent ? parent->id() : SID::INVALID, p->id());

    return p;
  }

  /// \copydoc nodeAt
  auto& nodeAt (SID i) {
    return const_cast<Node_ptr&>(std::as_const(*this).nodeAt(i));
  }

  /// \todo remove one
  /// \return Whether \p g is similar enough to \p species
  static float speciesMatchingScoreSimicontinuous (const GENOME &g,
                                                   Node_ptr species,
                                                   DCCache &dccache) {
    uint k = species->enveloppe.size();

    dccache.clear();
    dccache.reserve(k);

    uint matable = 0;
    for (const auto &ep: species->enveloppe) {
      double d = distance(g, ep.genome);
      double c = std::min(g.crossoverData()(d), ep.genome.crossoverData()(d));

      if (c >= Config::compatibilityThreshold()) matable++;
      dccache.push_back(d, c);
    }

    assert(dccache.size() == k);
    return matable - Config::similarityThreshold() * k;
  }

  /// \todo remove one
  /// \return Whether \p g is similar enough to \p species
  static float speciesMatchingScoreContinuous (const GENOME &g,
                                               Node_ptr species,
                                               DCCache &dccache) {
    uint k = species->enveloppe.size();

    dccache.clear();
    dccache.reserve(k);

    float avgCompat = 0;
    for (const auto &ep: species->enveloppe) {
      double d = distance(g, ep.genome);
      double c = std::min(g.crossoverData()(d), ep.genome.crossoverData()(d));

      avgCompat += c;
      dccache.push_back(d, c);
    }

    assert(dccache.size() == k);
    return avgCompat / float(k) - Config::avgCompatibilityThreshold();
  }

  /// \todo remove
  /// Proxy for delegating score computation to the appropriate function
  /// \see Config::FULL_CONTINUOUS
  static float speciesMatchingScore (const GENOME &g, Node_ptr species,
                                     DCCache &dccache) {
    auto f =
      Config::DEBUG_FULL_CONTINUOUS() ?
          speciesMatchingScoreContinuous
        : speciesMatchingScoreSimicontinuous;
    return f(g, species, dccache);
  }

  /// Finds the best derived species amongst the list of parents
  void findBestDerived (const GENOME &g, std::vector<Node_ptr> species,
                        Node_ptr &bestSpecies, float &bestScore,
                        DCCache &bestSpeciesDCCache) {

    DCCache dccache;
    std::vector<size_t> sizes;
    size_t max = 0;

    sizes.reserve(species.size());
    for (Node_ptr &sp: species) {
      size_t s = sp->children().size();
      sizes.push_back(s);
      if (max < s)  max = s;
    }

    if (debug() >= 2) std::cerr << "\tComputing scores:\n";
    for (uint i=0; i<max; i++) {
      for (uint k=0; k<species.size(); k++) {
        if (sizes[k] <= i)  continue;

        Node_ptr &subspecies = species[k]->child(i);
        float score = speciesMatchingScore(g, subspecies, dccache);

        if (debug() >= 2)
          std::cerr << "\t\t" << subspecies->id() << ": " << score << std::endl;

        if (bestScore < score) {
          bestSpecies = subspecies;
          bestScore = score;
          bestSpeciesDCCache = dccache;
        }

        if (bestScore > 0)
          return;
      }
    }
  }

  /// Find the appropriate place for \p g in the subtree(s) rooted at
  ///  \p species0 (and species1)
  /// \todo THis function seems ugly and hard to maintain
  std::pair<SID, UDATA*> addGenome (const GENOME &g,
                                    Node_ptr species0, Node_ptr species1,
                                    SID sid0, SID sid1) {

    if (debug()) {
      std::cerr << "Attempting to add genome " << g.crossoverData().id
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
    SpeciesContribution contrib;
    std::map<SID, float> scores;

    // Register first species
    species.push_back(species0);
    contrib.emplace_back(sid0, 1 + (sid0 == sid1));

    // Register (if needed) second species
    assert((species1 == nullptr) == (sid0 == sid1 || sid1 == SID::INVALID));
    if (species1) {
      species.push_back(species1);
      contrib.emplace_back(sid1, 1);
    }

    // Find best top-level species
    for (uint i=0; i<species.size(); i++) {
      Node_ptr s = species[i];
      float score = speciesMatchingScore(g, s, dccache);
      if (bestScore < score) {
        bestSpecies = s;
        bestScore = score;
        bestSpeciesDCCache = dccache;
      }
      scores[s->id()] = score;
    }

    // Order the contributions to put the best 'parent' first
    std::sort(contrib.begin(), contrib.end(),
              [&scores] (const Contribution &lhs, const Contribution &rhs) {
                return scores.at(lhs.species) >= scores.at(rhs.species);
    });

    if (debug() >= 2) {
      std::cerr << "\ttop-level scores:";
      for (auto &it: scores)
        std::cerr << " {" << it.first << ", " << it.second << "}";
      std::cerr << std::endl;
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
    findBestDerived(g, species, bestSpecies, bestScore, bestSpeciesDCCache);

    // Belongs to subspecies ?
    if (bestScore > 0) {
      if (debug())
        std::cerr << "\tCompatible with " << bestSpecies->id()
                  << " (score=" << bestScore << ")" << std::endl;
      return updateSpeciesContents(g, bestSpecies, bestSpeciesDCCache, contrib);

    } else if (debug())
      std::cerr << "\tIncompatible with all subspecies (score=" << bestScore
                << ")" << std::endl;

    // Need to create new species
    if (Config::simpleNewSpecies()) {
      Node_ptr subspecies = makeNode(contrib);
      dccache.clear();
      if (debug())
        std::cerr << "Created new species " << subspecies->id() << std::endl;
      return updateSpeciesContents(g, subspecies, dccache,
                                   SpeciesContribution{});

    } else
      assert(false);

    assert(false);
    return {SID::INVALID, nullptr};
  }

  /// Insert \p g into node \p species, possibly changing the enveloppe.
  ///
  /// Callbacks:
  ///   - Callbacks_t::onGenomeEntersEnveloppe
  ///   - Callbacks_t::onGenomeLeavesEnveloppe
  UDATA* insertInto (uint step, const GENOME &g, Node_ptr species,
                     const DCCache &dccache, Callbacks *callbacks) {

    using op = _details::DistanceMap::key_type;
    const uint k = species->enveloppe.size();

    auto &dist = species->distances;

    UDATA *userData = nullptr;

    // Populate the enveloppe
    if (k < _enveloppeSize) {
      if (debug())  std::cerr << "\tAppend to the enveloppe" << std::endl;

      species->enveloppe.push_back(Node::EnvPoint::make(g));
      userData = species->enveloppe.back().userData.get();
      if (callbacks)  callbacks->onGenomeEntersEnveloppe(species->id(), g.crossoverData().id);
      for (uint i=0; i<k; i++)
        dist[{i, k}] = dccache.distances[i];

    // Better enveloppe point ?
    } else {
      assert(k == _enveloppeSize);
      std::vector<GID> ids (k);
      for (uint i=0; i<k; i++)  ids[i] = species->enveloppePointId(i);
      _details::EnveloppeContribution ec =
          computeContribution(dist, dccache.distances, g.crossoverData().id, ids);

      // Genome inside the enveloppe. Nothing to do
      if (!ec.better) {
        if (debug())
          std::cerr << "\t" << g.crossoverData().id << "'s contribution is too low ("
                    << ec.value << ")" << std::endl;

      // Replace closest enveloppe point with new one
      } else {
        typename Node::EnvPoint &ep = species->enveloppe[ec.than];
        auto ep_id = ep.genome.crossoverData().id;

        if (debug())
          std::cerr << "\t" << g.crossoverData().id << "'s contribution is better "
                    << "than enveloppe point " << ec.than << " (id: "
                    << ep_id << ", c = " << ec.value << ")" << std::endl;

        if (callbacks) {
          callbacks->onGenomeLeavesEnveloppe(species->id(), ep_id);
          callbacks->onGenomeEntersEnveloppe(species->id(), g.crossoverData().id);
        }

        ep.userData->removedFromEnveloppe();
        userData = ep.userData.get();
        *ep.userData = UDATA(ep_id);

        ep.genome = g;
        for (uint i=0; i<k; i++)
          if (i != ec.than)
            dist[op{i,ec.than}] = dccache.distances[i];
      }
    }

    species->data.count++;
    species->data.currentlyAlive++;
    species->data.lastAppearance = step;

    return userData;
  }

  /// Update species \p s by inserting genome \p g, updating the contributions
  /// and registering the GID>SID association
  std::pair<SID, UDATA*>
  updateSpeciesContents(const GENOME &g, Node_ptr s,
                        const DCCache &cache,
                        const SpeciesContribution &ctb) {

    UDATA *userData = insertInto(_step, g, s, cache, _callbacks);
    if (!ctb.empty()) updateContributions(s, ctb);
    _idToSpecies.insert(g.crossoverData().id, s->id());
    return std::make_pair(s->id(), userData);
  }

  /// Update species \p s contributions with the provided values
  void updateContributions (Node_ptr s, const SpeciesContribution &contrib,
                            bool fromFile = false) {
    Node *oldMC = s->parent(),
         *newMC = s->update(contrib, _nodes);

    // No node (except the primordial species which cannot be re-assigned)
    // should be parentless. Except when creating a node
    assert(s->id() == SID(0) || oldMC || contrib.empty());

    if (oldMC != newMC) {
      assert(newMC);

      // Parent changed. Update and notify
      if (oldMC)  oldMC->delChild(s);
      newMC->addChild(s);

      if (!fromFile) {
        updateElligibilities();

#ifndef NDEBUG
      /// Check that no other nodes have changed their parent
        checkMC();
#endif

        _callbacks->onMajorContributorChanged(s->id(),
                                              oldMC->id(), newMC->id());
      }
    }
  }

  /// Triggers a tree-wide update of all contributors elligibility
  /// \todo remove test
  void updateElligibilities (void) {
    for (auto &p: _nodes) {
      Node_ptr &n = p.second;
      Node *oldMC = n->parent(),
           *newMC = n->updateElligibilities(_nodes);

      (void)oldMC;
      (void)newMC;
      assert(!oldMC || oldMC == newMC);
    }
  }

#ifndef NDEBUG
  /// Debug function
  void checkMC (void) {
    for (auto &p: _nodes) {
      Node_ptr &n = p.second;
      if (!n->parent()) {
        if (n->id() != SID(0))
          throw std::logic_error("Parent-less node is not valid (i.e. not 0)");
        continue;
      }

      auto pc = n->parent()->children();
      if (std::find(pc.begin(), pc.end(), n) == pc.end())
        throw std::logic_error("Node is not attached to the correct parent");
    }
  }
#endif

  /// Delete species with an underfilled enveloppe to limit clutter
  void performStillbornTrimming (void) {
    static const auto &T = Config::stillbornTrimmingThreshold();
    static const auto &D = Config::stillbornTrimmingDelay();
    static const float MD = Config::stillbornTrimmingMinDelay();

    bool remove = true;
    for (auto it = _nodes.begin(); it != _nodes.end();
         remove ? it = _nodes.erase(it) : ++it ) {

      Node &s = *it->second;
      remove = false;

      // Ignore non-leaf nodes
      if (!s.children().empty())  continue;

      // Only process dead species
      if (!s.extinct()) continue;

      bool underfilled = (s.enveloppe.size() < T * _enveloppeSize);
      uint liveTime = s.data.lastAppearance - s.data.firstAppearance;
      uint deadTime = _step - s.data.lastAppearance;
      if (underfilled && std::max(MD, liveTime * D) < deadTime) {
        if (Config::DEBUG_STILLBORNS()) {
          std::cerr << "Removing species " << s.id() << " with enveloppe size of "
                    << s.enveloppe.size() << " / " << _enveloppeSize << " ("
                    << 100. * s.enveloppe.size() / _enveloppeSize << "%) and "
                    << "survival time of " << " max(" << MD << ", " << D << " * ("
                    << s.data.lastAppearance << " - " << s.data.firstAppearance
                    << ")) = " << std::max(MD, D * liveTime) << " < " << deadTime
                    << " = " << _step << " - " << s.data.lastAppearance
                    << std::endl;
        }

        if (s.parent()) s.parent()->delChild(it->second);  // Erase from parent
        _stillborns++;
        remove = true;
      }
    }
  }

// =============================================================================
// == Json conversion

private:
  /// Serialize Node \p n into a json
  static json toJson (const Node &n) {
    json j, jd, jc;

    for (const auto &d: n.distances)
      jd.push_back({d.first.first, d.first.second, d.second});

    for (const auto &c: n.children())
      jc.push_back(toJson(*c));

    j["id"] = n.id();
    j["data"] = n.data;
    j["envlp"] = n.enveloppe;
    j["contribs"] = n.contributors.data();
    j["dists"] = jd;
    j["children"] = jc;

    return j;
  }

  /// Rebuilds PTree hierarchy and internal structure based on the contents of
  /// json \p j
  Node_ptr rebuildHierarchy(const json &j) {
    Contributors c (j["id"], j["contribs"]);
    Node_ptr n = Node::make_shared(c);

    _nodes[n->id()] = n;

    n->data = j["data"];
    n->enveloppe = j["envlp"].get<decltype(Node::enveloppe)>();
    const json jd = j["dists"];
    const json jc = j["children"];

    using op = _details::DistanceMap::key_type;
    for (const auto &d: jd)
      n->distances[op{d[0], d[1]}] = d[2];

    for (const auto &c: jc)
      rebuildHierarchy(c);

    return n;
  }

public:

  /// Serialise PTree \p pt into a json
  /// \arg complete Whether to include all data required to resume a simulation
  /// or only those used in displaying/analysing
  static void toJson (json &j, const PhylogenicTree &pt, bool complete) {
    j["_step"] = pt._step;
    j["_envSize"] = pt._enveloppeSize;
    j["_stillborns"] = pt._stillborns;
    j["tree"] = toJson(*pt._root);
    if (complete) {
      j["nextSID"] = pt._nextNodeID;
      j["i2s"] = pt._idToSpecies;
    }
  }

  /// Deserialise PTree \p pt from json \p j
  /// \arg complete Whether to load all data required to resume a simulation
  /// or only those used in displaying/analysing.
  static void fromJson (const json &j, PhylogenicTree &pt, bool complete) {
    pt._step = j["_step"];
    pt._stillborns = j["_stillborns"];
    pt._enveloppeSize = j["_envSize"];
    if (Config::enveloppeSize() != pt._enveloppeSize)
      utils::doThrow<std::invalid_argument>(
        "Current configuration file specifies an enveloppe size of ",
        Config::enveloppeSize(), " whereas the provided PTree was built with ",
        pt._enveloppeSize);

    pt._root = pt.rebuildHierarchy(j["tree"]);

    if (complete) {
      pt._nextNodeID = j["nextSID"];
      pt._idToSpecies = j["i2s"].get<IdToSpeciesMap>();
    }

    // Ensure correct parenting
    for (auto &n: pt._nodes)
      pt.updateContributions(n.second, {}, true);

#ifndef NDEBUG
    pt.checkMC();
#endif
  }

  /// Stores itself at the given location
  bool saveTo (const std::string &filename, bool complete = false) const {
    std::ofstream ofs (filename);
    if (!ofs) {
      std::cerr << "Unable to open '" << filename << "' for writing"
                << std::endl;
      return false;
    }

    json j;
    toJson(j, *this, complete);
    ofs << j.dump(2);
    return true;
  }

  /// \returns a phylogenic tree rebuilt from data at the given location
  static PhylogenicTree readFrom (const std::string &filename, bool complete = false) {
    std::ifstream ifs (filename);
    if (!ifs)
      throw std::invalid_argument ("Unable to open '" + filename
                                   + "' for reading");

    else {
      PhylogenicTree pt;
      json j = json::parse(utils::readAll(filename));
      fromJson(j, pt, complete);
      return pt;
    }
  }
};

} // end of namespace phylogeny

#endif // _PHYLOGENIC_TREE_H_
