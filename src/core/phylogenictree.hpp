#ifndef _PHYLOGENIC_TREE_H_
#define _PHYLOGENIC_TREE_H_

#include <vector>
#include <map>
#include <memory>
#include <fstream>

#include <cassert>
#include <iostream>

#include "kgd/external/json.hpp"
#include "treetypes.hpp"
#include "ptreeconfig.h"
#include "crossover.h"

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

  /// \brief Called when a species' major contributor has changed
  ///
  /// Provides the id of the species whose major contributor just changed
  /// as well as the id of the previous and new MC
  void onMajorContributorChanged (SID sid, SID oldMC, SID newMC);
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

namespace _details {

/// Distance & compatibilities cache
struct DCCache {
  /// Cache collection of distances
  std::vector<float> distances;

  /// Cache collection of compatibilities
  std::vector<float> compatibilities;

  /// Remove all contents
  void clear (void) {
    distances.clear(),  compatibilities.clear();
  }

  /// Prepare exactly \p n units of storage space
  void reserve (uint n) {
    distances.reserve(n), compatibilities.reserve(n);
  }

  /// Append values
  void push_back (float d, float c) {
    distances.push_back(d), compatibilities.push_back(c);
  }

  /// \returns the size of the cache
  size_t size (void) const {
    assert(distances.size() == compatibilities.size());
    return distances.size();
  }
};


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

/// Helper alias to a map whose keys are ordered so that
/// \f$\forall i,j \in M i<j\f$
using DistanceMap = std::map<ordered_pair<uint>, float>;

/// Description of the contribution of a genome to a species enveloppe
struct EnveloppeContribution {
  bool better;  ///< Should an enveloppe point be replaced
  uint than;    ///< Index of the enveloppe point to replace
  float value;  ///< Confidence of the replacement pertinence
};

/// Helper alias to a genome identificator
using GID = genotype::BOCData::GID;

/// Computes whether or not the considered species would be better described by
/// replacing a point from the current enveloppe (with distance map \p edist)
/// by an incoming genome (with distances \p gdist)
EnveloppeContribution computeContribution(const DistanceMap &edist,
                                          const std::vector<float> &gdist,
                                          GID gid, const std::vector<GID> &ids);

} // end of namespace _details

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

  /// \copydoc TreeTypes::SID
  using SID = TreeTypes::SID;

  /// Helper alias for the Parent enumeration
  using Parent = genotype::BOCData::Parent;

  /// Specialization used by this tree. Uses CRTP
  using Callbacks = Callbacks_t<PhylogenicTree<GENOME>>;

  /// Helper alias for the configuration data
  using Config = config::PTree;

  /// Helper alias to the type used to cache distance/compatibilities values
  using DCCache = _details::DCCache;

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
    TreeTypes::LivingSet aliveSpecies;
    for (IT it = begin; it != end; ++it)
      aliveSpecies.insert(_idToSpecies.at(geneticID(*it)));

    // Update internal data
    for (SID sid: aliveSpecies)
      nodeAt(sid)->data.lastAppearance = step;
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

    // Retrieve parent's species
    /// \todo this should not fail on initialization data (SID == -1)
    SID mSID = _idToSpecies.parentSID(g, Parent::MOTHER),
        fSID = _idToSpecies.parentSID(g, Parent::FATHER);

    Node *s0 = nullptr, *s1 = nullptr;
    if (mSID == SID::INVALID && fSID == SID::INVALID)
      s0 = _root;

    else if (fSID == SID::INVALID || mSID == fSID)
      s0 = _nodes[mSID];

    else {
      s0 = _nodes[mSID];
      s1 = _nodes[fSID];
    }

    Node *species = addGenome(g, s0, s1, {mSID, fSID});
    return species->id;
  }

  /// Remove \p g from this PTree (and update relevant internal data)
  SID delGenome (const GENOME &g) {
    auto sid = _idToSpecies.remove(g);
    if (Config::DEBUG())
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

  struct Node;

  /// Smart pointer (shared) to a species node
  using Node_ptr = std::shared_ptr<Node>;

  /// Helper type to the cached linear node collection
  using Nodes = TreeTypes::sidvector<Node_ptr>;

  /// Contributor field for a species node
  class NodeContributor {
    SID _speciesID;  ///< Reference to the contributor
    uint _count; ///< Number of contributions

  public:
    /// Constructor
    NodeContributor(Node_ptr species, uint initialCount)
      : _speciesID(species->id), _count(initialCount) {}

    /// Accessor to the contributor reference
    SID speciesID (void) const {
      return _speciesID;
    }

    /// Increment the number of contributions
    NodeContributor& operator+=(uint k) {
      _count += k;
      return *this;
    }

    /// Compare according to the respective number of contributions
    friend bool operator< (const NodeContributor &lhs, const NodeContributor &rhs) {
      // bigger contributions go first in the array
      return !(lhs._count < rhs._count);
    }
  };

  /// Sorted collection of contributors for a species node
  class Contributors {
    /// The buffer containing the individual contributions
    std::vector<NodeContributor> vec;

    /// The associated node
    Node *node;

    /// The main contributor
    Node *main;

    /// Update the main contributor cached variable
    void updateMain (const Nodes &nodes) {
      assert(!vec.empty());

      uint i=0;
      SID sid = SID::INVALID;
      while(sid != node->id && i<vec.size())
        sid = vec[i++].speciesID();
      main = (sid == SID::INVALID ? nullptr : nodes[sid].get());
    }

  public:
    /// Alias for the data structure containing the contributing SIDs
    using Contribution = std::multiset<SID>;

    /// Constructor. Registers the node whose contributor collection it manages
    Contributors (Node *n) : node(n), main(nullptr) {}

    /// \return The main contributor, i.e. with highest number of contributions
    Node* getMain (void) const {
      return main;
    }

    /// Register new contributions, updates internal data and returns the new
    /// main contributor
    Node* update (Contribution sids, const Nodes &nodes) {
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
        auto it = sids.begin();
        SID sid = *it;
        uint k = 1;
        while (*(++it) == sid)  k++;
        sids.erase(sid);

        vec.emplace_back(nodes[sid], k);
      }

      // sort by decreasing contribution
      std::stable_sort(vec.begin(), vec.end());

      updateMain(nodes);
      return main;
    }
  };

  /// \copydoc Contributors::Contribution
  using SpeciesContribution = typename Contributors::Contribution;

  /// Species node
  struct Node {
    SID id; ///< Species identificator
    SpeciesData data; ///< Species additionnal data

    /// Collection of contributors the this species' gene pool
    Contributors contributors;

    /// Subspecies of this node
    std::vector<Node_ptr> children;

    /// Collection of borderoids (in opposition to centroids)
    std::vector<GENOME> enveloppe;

    /// Cache map for the intra-enveloppe distances
    _details::DistanceMap distances;

    /// Create a node with \p id and \p parent
    Node (SID id) : id(id), contributors(this) {}

    /// \returns the main contributor for this species (excluding itself)
    Node* parent (void) {
      return contributors.getMain();
    }

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
  Nodes _nodes;

  /// Genome to species lookup table
  IdToSpeciesMap _idToSpecies;

  /// Pointer to the callbacks object. Null by default
  mutable Callbacks *_callbacks;

  uint _hybrids;  ///< Number of hybrids encountered thus far \todo WIP
  uint _step; ///< Current timestep for this tree

  /// Create a smart pointer to a node created on-the-fly under \p parent
  /// Callbacks:
  ///   - Callbacks_t::onNewSpecies
  Node_ptr makeNode (const SpeciesContribution &initialContrib) {
    Node_ptr p = std::make_shared<Node>(nextNodeID());
    p->data.firstAppearance = _step;
    p->data.lastAppearance = _step;
    p->data.count = 1;
    p->contributors.update(initialContrib, _nodes);

    _nodes.push_back(p);

    Node *parent = p->parent();
    if (parent) parent->children.push_back(p);
    if (_callbacks)  _callbacks->onNewSpecies(parent ? parent->id : SID::INVALID, p->id);

    return p;
  }

  /// \copydoc nodeAt
  auto& nodeAt (SID i) {
    return const_cast<Node_ptr&>(std::as_const(*this).nodeAt(i));
  }

  /// Find the appropriate place for \p g in the subtree(s) rooted at
  ///  \p species0 (and species1)
  Node* addGenome (const GENOME &g, Node *species0, Node *species1,
                   const SpeciesContribution &contrib) {

    if (Config::DEBUG()) {
      std::cerr << "Attempting to add genome " << g.id()
                << " to species ";
      if (!species1)
        std::cerr << species0->id;
      else
        std::cerr << "either " << species0->id << " or " << species1->id;
      std::cerr << std::endl;
    }

    DCCache dccache, bestSpeciesDCCache;
    Node *bestSpecies = nullptr;
    float bestScore = -std::numeric_limits<float>::max();

    std::vector<Node*> species;
    species.push_back(species0);
    if (species1) species.push_back(species0);

    // Find best top-level species
    for (Node *s: species) {
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

    if (Config::DEBUG()) {
      std::cerr << "\tIncompatible with ";
      if (!species1)  std::cerr << species0->id;
      else  std::cerr << "both " << species0->id << " and " << species1->id;
      std::cerr << std::endl;
    }

    // Find best derived species
    uint i=0, k = 0, max = species0->children.size();
    if (species1) max = std::max(max, species1->children.size());
    while (bestScore <= 0 && i < max) {
      Node_ptr &subspecies = species[k]->children[i];
      float score = speciesMatchingScore(g, subspecies, dccache);
      if (bestScore < score) {
        bestSpecies = subspecies;
        bestScore = score;
        bestSpeciesDCCache = dccache;
      }

      if (!species1 || species[1-k]->size() <= i)
        i++;
      else {
        if (k==1) i++;
        k = 1 - k;
      }
    }

    // Belongs to subspecies ?
    if (bestScore > 0)
      return updateSpeciesContents(g, bestSpecies, bestSpeciesDCCache, contrib);

    // Need to create new species
    if (Config::simpleNewSpecies()) {
      Node_ptr subspecies = makeNode(contrib);
      dccache.clear();

      insertInto(_step, g, subspecies, dccache, _callbacks);
      _idToSpecies.insert(g.id(), subspecies->id);
      return subspecies;

    } else
      assert(false);

    assert(false);
    return nullptr;
  }

  /// \return Whether \p g is similar enough to \p species
  friend float speciesMatchingScore (const GENOME &g, Node_ptr species,
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
    if (k < Config::enveloppeSize()) {
      if (Config::DEBUG())  std::cerr << "\tAppend to the enveloppe" << std::endl;

      species->enveloppe.push_back(g);
      if (callbacks)  callbacks->onGenomeEntersEnveloppe(species->id, g.id());
      for (uint i=0; i<k; i++)
        dist[{i, k}] = dccache.distances[i];

    // Better enveloppe point ?
    } else {
      assert(k == Config::enveloppeSize());
      std::vector<GID> ids (k);
      for (uint i=0; i<k; i++)  ids[i] = species->enveloppe[i].id();
      _details::EnveloppeContribution ec =
          computeContribution(dist, dccache.distances, g.id(), ids);

      // Genome inside the enveloppe. Nothing to do
      if (!ec.better) {
        if (Config::DEBUG())
          std::cerr << "\t" << g.id() << "'s contribution is too low ("
                    << ec.value << ")" << std::endl;

      // Replace closest enveloppe point with new one
      } else {
        if (Config::DEBUG())
          std::cerr << "\t" << g.id() << "'s contribution is better "
                    << "than enveloppe point " << ec.than << " (id: "
                    << species->enveloppe[ec.than].id()
                    << ", c = " << ec.value << ")" << std::endl;

        if (callbacks) {
          callbacks->onGenomeLeavesEnveloppe(species->id, species->enveloppe[ec.than].id());
          callbacks->onGenomeEntersEnveloppe(species->id, g.id());
        }
        species->enveloppe[ec.than] = g;
        for (uint i=0; i<k; i++)
          if (i != ec.than)
            dist[{i,ec.than}] = dccache.distances[i];
      }
    }

    species->data.count++;
    species->data.lastAppearance = step;
  }

  /// Update species \p s by inserting genome \p g, updating the contributions
  /// and registering the GID>SID association
  Node_ptr updateSpeciesContents(const GENOME &g, Node_ptr s,
                                 DCCache &cache,
                                 const SpeciesContribution &ctb) {

    insertInto(_step, g, s, cache, _callbacks);
    updateContributions(s, ctb);
    _idToSpecies.insert(g.id(), s->id);
    return s;
  }

  /// Update species \p s contributions with the provided values
  void updateContributions (Node_ptr s, const SpeciesContribution &contrib) {
    Node *oldMC = s->parent(),
         *newMC = s->contributors.update(contrib, _nodes);

    if (oldMC != newMC) {
      // Parent changed. Update and notify
      if (oldMC) {
        auto &v = oldMC->children;
        v.erase(std::remove(v.begin(), v.end(), s), v.end());
      }

      newMC->children.push_back(s);

      _callbacks->onMajorContributorChanged(s->id, oldMC->id, newMC->id);
    }
  }

// =============================================================================
// == Json conversion

private:
  /// Rebuilds PTree hierarchy and internal structure based on the contents of
  /// json \p j
  Node_ptr rebuildHierarchy(Node_ptr parent, const nlohmann::json &j) {
    Node_ptr n = makeNode(parent);
    /// \todo careful not to forget to maintain contributors

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

    /// \todo careful not to forget to maintain contributors

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
