#ifndef KGD_PHYLOGENIC_TREE_H
#define KGD_PHYLOGENIC_TREE_H

#include <vector>
#include <map>
#include <memory>
#include <fstream>
#include <bitset>

#include <cassert>
#include <iostream>

#include "../ptreeconfig.h"

#include "treetypes.h"
#include "node.hpp"
#include "callbacks.hpp"

/*!
 * \file phylogenetictree.hpp
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
template <typename GENOME, typename UDATA>
class PhylogeneticTree {
  /// Helper lambda for debug printing
  static constexpr auto debug = [] {
    return config::PTree::DEBUG_LEVEL() * config::PTree::DEBUG_PTREE();
  };

// =============================================================================
// == Helper types
public:
  /// Helper alias to the genome type template parameter
  using Genome = GENOME;

  /// Helper alias to the genome type template parameter
  using UserData = UDATA;

  /// Helper alias to a species node
  using Node = phylogeny::Node<Genome, UserData>;

  /// \copydoc Node::Ptr
  using Node_ptr = typename Node::Ptr;

  /// \copydoc Node::Collection
  using Nodes = typename Node::Collection;

  /// Specialization used by this tree. Uses CRTP
  using Callbacks = Callbacks_t<PhylogeneticTree<Genome, UserData>>;

  /// Helper alias for the configuration data
  using Config = config::PTree;

  /// Helper alias to the type used to cache distance/compatibilities values
  using DCCache = _details::DCCache;

  /// \copydoc Contributors::Contributions
  using SpeciesContribution = typename Contributors::Contributions;

  /// \copydoc phylogeny::InsertionResult
  using InsertionResult = phylogeny::InsertionResult<UserData>;

// =============================================================================
// == Resource management (creation, destruction, copy)

  /// Create an empty PTree
  PhylogeneticTree(void) {
    _nextNodeID = SID(0);
    _rsetSize = Config::rsetSize();
    _stillborns = 0;
    _step = 0;
    _root = nullptr;
    _callbacks = nullptr;
  }

  /// Constructs a deep copy of that PTree
  PhylogeneticTree (const PhylogeneticTree &that) {
    _nextNodeID = that._nextNodeID;

    _root = deepcopy(that._root);
    updateElligibilities();

    _callbacks = nullptr;

    _rsetSize = that._rsetSize;
    _stillborns = that._stillborns;
    _step = that._step;
  }

  /// Assigns that PTree to this one
  PhylogeneticTree& operator= (PhylogeneticTree that) {
    swap(*this, that);
    return *this;
  }

  /// Nothing to do. All is based on smart-pointers.
  ~PhylogeneticTree (void) {}

private:
  /// Performs a deepcopy of that_n node and all descendants into this PTree
  Node_ptr deepcopy (const Node_ptr &that_n) {
    Node_ptr this_n = Node::make_shared(that_n->contributors);

    this_n->data = that_n->data;
    this_n->rset = that_n->rset;
    this_n->distances = that_n->distances;

    _nodes[this_n->id()] = this_n;

    for (const Node_ptr &that_c: that_n->children())
      this_n->addChild(deepcopy(that_c));

    return this_n;
  }

  /// Swap contents of the provided phylogenetic trees
  friend void swap (PhylogeneticTree &lhs, PhylogeneticTree &rhs) {
    using std::swap;
    swap(lhs._nextNodeID, rhs._nextNodeID);
    swap(lhs._root, rhs._root);
    swap(lhs._nodes, rhs._nodes);
    swap(lhs._callbacks, rhs._callbacks);
    swap(lhs._rsetSize, rhs._rsetSize);
    swap(lhs._stillborns, rhs._stillborns);
    swap(lhs._step, rhs._step);
  }

public:

// =============================================================================
// == Accessors

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
    if (i == SID::INVALID)
      utils::Thrower("SID::INVALID is (by definition) invalid");

    auto it = _nodes.find(i);
    if (it == _nodes.end())
      utils::Thrower("No node found for species ", i);

    return it->second;
  }

  /// \return the user data for enveloppe point \p gid or nullptr if it is a
  /// regular individual
  UserData* getUserData (const PID &pid) const {
    const auto &species = nodeAt(pid.sid);
    for (const auto &ep: species->rset)
      if (ep.genome.genealogy().self.gid == pid.gid)
        return ep.userData.get();
    return nullptr;
  }

  /// \return the current timestep for this PTree
  uint step (void) const {
    return _step;
  }

  /// Access current value without modifying it
  SID nextNodeID (void) const {
    return _nextNodeID;
  }

  /// Access current set of alive species ids
  const LivingSet& aliveSpecies (void) const {
    return _aliveSpecies;
  }

protected:
  /// \copydoc nodeAt
  auto& nodeAt (SID i) {
    return const_cast<Node_ptr&>(std::as_const(*this).nodeAt(i));
  }

public:

// =============================================================================
// == Modifiers

  /// Sets the callbacks used by this ptree
  void setCallbacks (Callbacks *c) const { _callbacks = c; }

  /// Sets the current timestep for this PTree
  void setStep (uint step) {
    _step = step;
  }

protected:
  /// Wraps incrementation of the species identificator counter
  SID nextNodeID (void) {
    using SID_t = std::underlying_type<SID>::type;
    SID curr = _nextNodeID;
    _nextNodeID = SID(SID_t(_nextNodeID)+1);
    return curr;
  }

public:

// =============================================================================
// == Core API

  /// Update the set of still-alive species based on the list of still-alive
  /// genomes [\p begin,\p end[ extracted through \p geneticID
  ///
  /// Callbacks:
  ///   - Callbacks_t::onStepped
  ///
  /// \tparam IT Iterator to the begin/end of the population list
  /// \tparam F Functor for extracting the genome id from an iterator
  template <typename IT, typename F>
  void step (uint step, IT begin, IT end, F sidExtractor) {
    // Determine which species are still alive
    _aliveSpecies.clear();
    for (IT it = begin; it != end; ++it)
      _aliveSpecies.insert(sidExtractor(*it));

    // Update internal data
    for (SID sid: _aliveSpecies)
      nodeAt(sid)->data.lastAppearance = step;
    _step = step;

    static const auto &T = Config::stillbornTrimmingPeriod();
    if ((T > 0) && (_step % T) == 0)  performStillbornTrimming();

    // Potentially notify outside world
    if (_callbacks) _callbacks->onStepped(step, _aliveSpecies);
  }

  /// Insert \p g into this PTree
  /// \return The species \p g was added to and, if it was also added to the
  /// enveloppe, a pointer to the associated user data structure
  /// (nullptr otherwise).
  InsertionResult addGenome (const Genome &g) {
    // Ensure that the root exists
    if (!_root) {
      _root = makeNode(SpeciesContribution{});
      return updateSpeciesContents(g, _root, DCCache{}, SpeciesContribution{});
    }

    // Retrieve parent's species
    const Genealogy &genealogy = g.genealogy();
    SID mSID = genealogy.mother.sid,
        fSID = genealogy.father.sid;

    Node_ptr s0 = nullptr, s1 = nullptr;
    if (mSID == SID::INVALID && fSID == SID::INVALID)
      s0 = _root;

    else if (fSID == SID::INVALID || mSID == fSID)
      s0 = nodeAt(mSID);

    else {
      s0 = nodeAt(mSID);
      s1 = nodeAt(fSID);
    }

    // Remove (now obsolete) candidacies
    if (s0->data.pendingCandidates > 0)  s0->data.pendingCandidates--;
    if (s1 && s1->data.pendingCandidates > 0)  s1->data.pendingCandidates--;

    auto ret = addGenome(g, s0, s1, mSID, fSID);

    _stats.insertions++;

    if (Config::DEBUG_LEVEL())  std::cerr << std::endl;
    return ret;
  }

  /// Remove \p g from this PTree (and update relevant internal data)
  void delGenome (const Genome &g) {
    SID sid = g.genealogy().self.sid;

    if (debug())
      std::cerr << "New last appearance of species " << sid << " is " << _step
                << std::endl;

    SpeciesData &data = nodeAt(sid)->data;
    data.lastAppearance = _step;
    data.currentlyAlive--;

    _stats.deletions++;
  }

  /// Register candidate for future insertion attempt in either (sub)species
  void registerCandidate (const Genealogy &g) {
    performCandidacyRegistration(g, +1);
  }

  /// Unregister candidate that will not, after all, attempt insertion in either
  /// (sub)species.
  /// \warning Implies a previous call to registerCandidate() with the same \p g
  void unregisterCandidate (const Genealogy &g) {
    performCandidacyRegistration(g, -1);
  }

// =============================================================================
// == Stats management

  /// Ad-oc structure for printing the stats header
  struct StatsHeader {
    /// Prints the stats header
    friend std::ostream& operator<< (std::ostream &os, const StatsHeader&) {
      return os << " PTInsertions PTDeletions PTComparisons PTBranching";
    }
  };

  /// Structure for storing statistics about the phylogenetic dynamics
  struct Stats {
    uint insertions = 0;  ///< Number of genomes inserted
    uint deletions = 0;   ///< Number of genomes removed
    uint comparisons = 0; ///< Number of representatives tested
    uint branching = 0;   ///< Number of subspecies at root points

    /// Inserts provided stats in a default fashion
    friend std::ostream& operator<< (std::ostream &os, const Stats &s) {
      return os << " " << s.insertions << " " << s.deletions << " "
                << s.comparisons << " " << s.branching;
    }

  } _stats; ///< Field storing the phylogenetic dynamics

  /// Resets the statistics
  void resetStats (void) {
    _stats = Stats{};
  }

  /// \return the phylogenetic dynamics statistics
  const Stats& stats (void) const {
    return _stats;
  }

// =============================================================================
// == Member variables

private:
  /// Identificator for the next species
  SID _nextNodeID;

protected:
  /// The PTree root. Null until the first genome is inserted
  Node_ptr _root;

  /// Nodes collection for logarithmic access
  Nodes _nodes;

  /// Set of currently alive species
  LivingSet _aliveSpecies;

  uint _rsetSize;  ///< Number of enveloppe points
  uint _stillborns; ///< Number of stillborn species removed
  uint _step; ///< Current timestep for this tree

  /// Pointer to the callbacks object. Null by default
  mutable Callbacks *_callbacks;

// =============================================================================
// == Helper functions

  /// Create a smart pointer to a node created on-the-fly with contributors
  /// as described in \p initialContrib
  /// Callbacks:
  ///   - Callbacks_t::onNewSpecies
  Node_ptr makeNode (const SpeciesContribution &contrib) {

    SID id = nextNodeID();
    Contributors c (id);

    Node_ptr p = Node::make_shared(c);
    assert(p);

    p->data.firstAppearance = _step;
    p->data.lastAppearance = _step;
    p->data.count = 0;
    p->data.currentlyAlive = 0;
    p->data.pendingCandidates = 0;

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

  /// \todo remove one
  /// \return Whether \p g is similar enough to \p species
  static float speciesMatchingScoreSimicontinuous (const Genome &g,
                                                   Node_ptr species,
                                                   DCCache &dccache,
                                                   Stats &stats) {
    uint k = species->rset.size();

    dccache.clear();
    dccache.reserve(k);

    uint matable = 0;
    for (const auto &ep: species->rset) {
      double d = distance(g, ep.genome);
      double c = std::min(g.compatibility(d), ep.genome.compatibility(d));

      stats.comparisons++;

      if (c >= Config::compatibilityThreshold()) matable++;
      dccache.push_back(d, c);
    }

    assert(dccache.size() == k);
    return matable - Config::similarityThreshold() * k;
  }

  /// \todo remove one
  /// \return Whether \p g is similar enough to \p species
  static float speciesMatchingScoreContinuous (const Genome &g,
                                               Node_ptr species,
                                               DCCache &dccache,
                                               Stats &stats) {
    uint k = species->rset.size();

    dccache.clear();
    dccache.reserve(k);

    float avgCompat = 0;
    for (const auto &ep: species->rset) {
      double d = distance(g, ep.genome);
      double c = std::min(g.compatibility(d), ep.genome.compatibility(d));

      stats.comparisons++;

      avgCompat += c;
      dccache.push_back(d, c);
    }

    assert(dccache.size() == k);
    return avgCompat / float(k) - Config::avgCompatibilityThreshold();
  }

  /// \todo remove
  /// Proxy for delegating score computation to the appropriate function
  /// \see Config::FULL_CONTINUOUS
  static float speciesMatchingScore (const Genome &g, Node_ptr species,
                                     DCCache &dccache, Stats &stats) {
    auto f =
      Config::DEBUG_FULL_CONTINUOUS() ?
          speciesMatchingScoreContinuous
        : speciesMatchingScoreSimicontinuous;
    return f(g, species, dccache, stats);
  }

  /// Finds the best derived species amongst the list of parents
  void findBestDerived (const Genome &g, const std::vector<Node_ptr> &species,
                        Node_ptr &bestSpecies, float &bestScore,
                        DCCache &bestSpeciesDCCache) {

    DCCache dccache;

    const auto S = species.size();

    using it_t = decltype(std::declval<Node>().children().rbegin());
    std::vector<it_t> its, ends;
    its.reserve(S);
    ends.reserve(S);
    for (const Node_ptr &sp: species) {
      its.push_back(sp->children().crbegin());
      ends.push_back(sp->children().crend());
    }

    const uint allDone = [&species] {
      uint tmp = 0;
      for (uint i=0; i<species.size(); i++) tmp |= (1<<i);
      return tmp;
    }();

    uint done = 0;
    int k = 0;

    if (debug() >= 2) std::cerr << "\tComputing scores:\n";
    while (done != allDone) {
      auto &it = its[k];
      if (it == ends[k]) {
        done |= (1<<k);

      } else {
        _stats.branching++;

        const Node_ptr &subspecies = *it;
        float score = speciesMatchingScore(g, subspecies, dccache, _stats);

        if (debug() >= 2)
          std::cerr << "\t\t" << subspecies->id() << ": " << score << std::endl;

        if (bestScore < score) {
          bestSpecies = subspecies;
          bestScore = score;
          bestSpeciesDCCache = dccache;
        }

        if (bestScore > 0)
          return;

        ++it;
      }

      k = (k + 1) % S;
    }
  }

  /// Find the appropriate place for \p g in the subtree(s) rooted at
  ///  \p species0 (and species1)
  /// \todo THis function seems ugly and hard to maintain
  InsertionResult addGenome (const Genome &g,
                             Node_ptr species0, Node_ptr species1,
                             SID sid0, SID sid1) {

    if (debug()) {
      std::cerr << "Attempting to add genome " << g.genealogy().self.gid
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
      float score = speciesMatchingScore(g, s, dccache, _stats);
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
  UserData* insertInto (uint step, const Genome &g, Node_ptr species,
                     const DCCache &dccache, Callbacks *callbacks) {

    using op = _details::DistanceMap::key_type;
    const uint k = species->rset.size();

    auto &dist = species->distances;

    UserData *userData = nullptr;

    // Populate the enveloppe
    if (k < _rsetSize) {
      if (debug())  std::cerr << "\tAppend to the enveloppe" << std::endl;

      species->rset.push_back(Node::Representative::make(g));
      userData = species->rset.back().userData.get();
      species->rset.back().timestamp = _step;
      if (callbacks)  callbacks->onGenomeEntersEnveloppe(species->id(),
                                                         g.genealogy().self.gid);
      for (uint i=0; i<k; i++)
        dist[{i, k}] = dccache.distances[i];

    // Better enveloppe point ?
    } else {
      assert(k == _rsetSize);
      std::vector<GID> ids (k);
      for (uint i=0; i<k; i++)  ids[i] = species->representativeId(i);
      _details::EnveloppeContribution ec =
          computeContribution(dist, dccache.distances, g.genealogy().self.gid, ids);

      // Genome inside the enveloppe. Nothing to do
      if (!ec.better) {
        if (debug())
          std::cerr << "\t" << g.genealogy().self.gid << "'s contribution is too low ("
                    << ec.value << ")" << std::endl;

      // Replace closest enveloppe point with new one
      } else {
        typename Node::Representative &ep = species->rset[ec.than];
        auto ep_id = ep.genome.genealogy().self.gid;

        if (debug())
          std::cerr << "\t" << g.genealogy().self.gid << "'s contribution is better "
                    << "than enveloppe point " << ec.than << " (id: "
                    << ep_id << ", c = " << ec.value << ")" << std::endl;

        if (callbacks) {
          callbacks->onGenomeLeavesEnveloppe(species->id(), ep_id);
          callbacks->onGenomeEntersEnveloppe(species->id(), g.genealogy().self.gid);
        }

        ep.userData->removedFromEnveloppe();
        userData = ep.userData.get();
        *ep.userData = UserData(ep_id);

        ep.genome = g;
        for (uint i=0; i<k; i++)
          if (i != ec.than)
            dist[op{i,ec.than}] = dccache.distances[i];

        ep.timestamp = _step;
      }
    }

    species->data.count++;
    species->data.currentlyAlive++;
    species->data.lastAppearance = step;

    return userData;
  }

  /// Update species \p s by inserting genome \p g, updating the contributions
  /// and registering the GID>SID association in the genome's dedicated field
  InsertionResult
  updateSpeciesContents(const Genome &g, Node_ptr s,
                        const DCCache &cache,
                        const SpeciesContribution &ctb) {

    UserData *userData = insertInto(_step, g, s, cache, _callbacks);
    if (!ctb.empty()) updateContributions(s, ctb);
    return InsertionResult{s->id(), userData};
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

  /// Actually updates the candidacy values
  void performCandidacyRegistration (const Genealogy &g, int dir) {
    SID mSID = g.mother.sid, fSID = g.father.sid;
    if (mSID != SID::INVALID)
      nodeAt(mSID)->data.pendingCandidates += dir;
    if (mSID != fSID && fSID != SID::INVALID)
      nodeAt(fSID)->data.pendingCandidates += dir;
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

//#pragma GCC push_options
//#pragma GCC optimize ("O0")
//#warning performStillbornTrimming() optimisation disabled
  /// Delete species with an underfilled enveloppe to limit clutter
  void performStillbornTrimming (void) {
    static const auto &T = Config::stillbornTrimmingThreshold();
    static const auto &D = Config::stillbornTrimmingDelay();
    static const float MD = Config::stillbornTrimmingMinDelay();

    if (Config::DEBUG_STILLBORNS())
      std::cerr << "Performing stillborn trimming for step "
                << _step << std::endl;

    bool remove = false;
    for (auto it = _nodes.begin(); it != _nodes.end();
         remove ? it = _nodes.erase(it) : ++it, remove = false ) {

      auto p = it->second;
      assert(p);

      Node &s = *p;

      // Ignore non-leaf nodes
      if (!s.children().empty())  continue;

      // Only process dead species
      if (!s.extinct()) continue;

      bool underfilled = (s.rset.size() < T * _rsetSize);
      uint liveTime = s.data.lastAppearance - s.data.firstAppearance;
      uint deadTime = _step - s.data.lastAppearance;
      if (underfilled && std::max(MD, liveTime * D) < deadTime) {
        if (Config::DEBUG_STILLBORNS()) {
          std::cerr << "Removing species " << s.id() << " with enveloppe size of "
                    << s.rset.size() << " / " << _rsetSize << " ("
                    << 100. * s.rset.size() / _rsetSize << "%) and "
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
//#pragma GCC pop_options

// =============================================================================
// == Generic printing

  /// Stream \p pt to \p os. Mostly for debugging purpose: output is quickly
  /// unintelligible
  friend std::ostream& operator<< (std::ostream &os, const PhylogeneticTree &pt) {
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
    j["envlp"] = n.rset;
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
    n->rset = j["envlp"].get<decltype(Node::rset)>();
    const json &jd = j["dists"];
    const json &jc = j["children"];

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
  static void toJson (json &j, const PhylogeneticTree &pt) {
    j["_step"] = pt._step;
    j["_envSize"] = pt._rsetSize;
    j["_stillborns"] = pt._stillborns;
    j["alive"] = pt._aliveSpecies;
    j["tree"] = toJson(*pt._root);
    j["nextSID"] = pt._nextNodeID;
  }

  /// Deserialise PTree \p pt from json \p j
  /// \arg complete Whether to load all data required to resume a simulation
  /// or only those used in displaying/analysing.
  static void fromJson (const json &j, PhylogeneticTree &pt) {
    pt._step = j["_step"];
    pt._stillborns = j["_stillborns"];
    pt._rsetSize = j["_envSize"];
    if (Config::rsetSize() != pt._rsetSize)
      utils::Thrower(
        "Current configuration file specifies an enveloppe size of ",
        Config::rsetSize(), " whereas the provided PTree was built with ",
        pt._rsetSize);

    pt._root = pt.rebuildHierarchy(j["tree"]);
    pt._aliveSpecies = j["alive"].get<LivingSet>();
    pt._nextNodeID = j["nextSID"];

    // Ensure correct parenting
    for (auto &n: pt._nodes)
      pt.updateContributions(n.second, {}, true);

#ifndef NDEBUG
    pt.checkMC();
#endif
  }

  /// Asserts that two phylogenetic trees are equal
  friend void assertEqual (const PhylogeneticTree &lhs,
                           const PhylogeneticTree &rhs, bool deepcopy) {
    using utils::assertEqual;
    assertEqual(lhs._root, rhs._root, deepcopy);
    assertEqual(lhs._nodes, rhs._nodes, deepcopy);
    assertEqual(lhs._aliveSpecies, rhs._aliveSpecies, deepcopy);

    assertEqual(lhs._nextNodeID, rhs._nextNodeID, deepcopy);
    assertEqual(lhs._rsetSize, rhs._rsetSize, deepcopy);
    assertEqual(lhs._stillborns, rhs._stillborns, deepcopy);
    assertEqual(lhs._step, rhs._step, deepcopy);
  }

  /// Stores itself at the given location
  bool saveTo (const stdfs::path &filename) const {
    std::ofstream ofs (filename);
    if (!ofs) {
      std::cerr << "Unable to open '" << filename << "' for writing"
                << std::endl;
      return false;
    }

    saveTo(ofs, -1);
    return true;
  }

  /// Stores itself in the provided stream
  void saveTo (std::ostream &os, uint ident = 2) const {
    json j;
    toJson(j, *this);
    os << j.dump(ident);
  }

  /// \returns a phylogenic tree rebuilt from data at the given location
  static PhylogeneticTree readFrom (const std::string &filename) {
    std::ifstream ifs (filename);
    if (!ifs)
      throw std::invalid_argument ("Unable to open '" + filename
                                   + "' for reading");

    else {
      PhylogeneticTree pt;
      json j = json::parse(utils::readAll(filename));
      fromJson(j, pt);
      return pt;
    }
  }
};

} // end of namespace phylogeny

#endif // KGD_PHYLOGENIC_TREE_H
