#ifndef NODE_HPP
#define NODE_HPP

/*!
 * \file node.hpp
 *
 * Contains the definition for a single species node in the phylogenic tree
 */

#include "speciesdata.hpp"
#include "speciescontributors.h"

namespace phylogeny {

/// Species node
template <typename GENOME, typename UDATA>
struct Node {
  /// Helper alias to the type used for a pointer to node
  using Ptr = std::shared_ptr<Node>;

  /// Helper alias to a collection of nodes
  using Collection = std::map<SID, Ptr>;

  /// Stores the data relative to an enveloppe point
  struct Representative {
    uint timestamp; ///< Insertion date
    GENOME genome;  ///< The genome for this representant
    std::unique_ptr<UDATA> userData; ///< Associated user managed statistics

    /// Default constructor
    Representative (void) = default;

    /// Copy constructs a representative. User data is deep copied.
    Representative (const Representative &that)
      : genome(that.genome),
        userData(std::make_unique<UDATA>(*that.userData)) {}

    /// Defaulted move constructor
    Representative (Representative &&) = default;

    /// Assign another representative to this
    Representative& operator= (Representative that) {
      swap(*this, that);
      return *this;
    }

    /// Creates the enveloppe point for genome \p g and default-initialize
    /// the associated user data
    static Representative make (const GENOME &g) {
      return Representative(g);
    }

    /// Serialize enveloppe point \p p into a json
    friend void to_json (json &j, const Representative &p) {
      j = { p.genome, *p.userData };
    }

    /// Deserialize enveloppe point \p p from a json
    friend void from_json (const json &j, Representative &p) {
      p.genome = j[0].get<GENOME>();
      p.userData = std::make_unique<UDATA>(phylogeny::GID::INVALID);
      *p.userData = j[1].get<UDATA>();
    }

    /// Asserts that two enveloppe points are equal
    friend void assertEqual (const Representative &lhs,
                             const Representative &rhs, bool deepcopy) {
      using utils::assertEqual;
      assertEqual(lhs.genome, rhs.genome, deepcopy);
      assertEqual(lhs.userData, rhs.userData, deepcopy);
    }

  private:
    /// Creates a representative of the provided genome
    Representative (const GENOME &genome)
      : genome(genome),
        userData(std::make_unique<UDATA>(genome.genealogy().self.gid)) {}

    /// Swaps contents of the two representatives
    void swap (Representative &lhs, Representative &rhs) {
      using std::swap;
      swap(lhs.genome, rhs.genome);
      swap(lhs.userData, rhs.userData);
    }
  };

private:
  /// Prevents constructor access from outside the class
  struct cookie {};

  /// Reference to the species' parent (main contributor)
  Node *_parent;

  /// Subspecies of this node
  std::vector<Ptr> _children;

public:
  SpeciesData data; ///< Species additionnal data

  /// Collection of contributors the this species' gene pool
  Contributors contributors;

  /// Collection of borderoids (in opposition to centroids)
  std::vector<Representative> rset;

  /// Cache map for the intra-enveloppe distances
  _details::DistanceMap distances;

  /// Creates a node from a contributors collection (hidden from user. use the
  /// make_shared version)
  explicit Node (Contributors &&contribs, const cookie&)
    : _parent(nullptr), contributors(contribs) {}

  /// \returns a pointer to a newly allocated node created from the provided
  /// arguments
  template <typename ...ARGS>
  static Ptr make_shared (ARGS... args) {
    return std::make_shared<Node>(std::forward<ARGS>(args)..., cookie{});
  }

  /// \returns the species identificator for this node
  SID id (void) const {
    return contributors.getNodeID();
  }

  /// \returns the main contributor for this species (excluding itself)
  Node* parent (void) {
    return _parent;
  }

  /// \returns the main contributor for this species (excluding itself)
  const Node* parent (void) const {
    return _parent;
  }

  /// \returns the collection of subspecies root at this node
  const auto& children (void) const {
    return _children;
  }

  /// \returns the subspecies at index \p i
  auto& child (size_t i) {
    return _children[i];
  }

  /// \returns the genome of representative \p i
  auto representativeGenome (uint i) const {
    return rset[i].genome;
  }

  /// \returns the genetic identificator for representative \p i
  auto representativeId (uint i) const {
    return representativeGenome(i).genealogy().self.gid;
  }

  /// \returns whether this species still has some members in the simulation
  bool extinct (void) const {
    return data.currentlyAlive == 0 && data.pendingCandidates == 0;
  }

  /// Adds subspecies \p child to this node
  void addChild (Ptr child) {
    _children.push_back(child);
  }

  /// Removes subspecies \p child from this node
  void delChild (Ptr child) {
    _children.erase(std::remove(_children.begin(), _children.end(), child),
                    _children.end());
  }

  /// Helper function generating a lambda binded to the provided collection
  /// \p nodes
  static auto elligibilityTester (const Collection &nodes) {
    using namespace std::placeholders;
    return std::bind(&Contributors::elligibile<Collection>,
                     _1, _2, nodes);
  }

  /// Updates the species contributions manager and the species' main parent
  /// \returns the new species' main parent
  Node* update (Contributors::Contributions sids, const Collection &nodes) {
    SID mainSID = contributors.update(sids, elligibilityTester(nodes));
    return updateParent(mainSID, nodes);
  }

  /// Triggers a tree-wide recomputation of the elligibilities of all node
  /// contributors. Possible chain-reaction (multiple parent modifications)
  ///
  /// \returns the current (possibly changed?) parent
  Node* updateElligibilities (const Collection &nodes) {
    SID mainSID = contributors.updateElligibilities(elligibilityTester(nodes));
    return updateParent(mainSID, nodes);
  }


  /// Stream this node. Mostly for debugging purpose.
  friend std::ostream& operator<< (std::ostream &os, const Node &n) {
    std::string spacing = "> ";
    const Node *p = &n;
    while ((p = p->_parent))  spacing += "  ";

    os << spacing << "[" << n.id << "] ( ";
    for (const Representative &p: n.rset)    os << p.genome.id() << " ";
    os << ")\n";

    for (const Ptr &ss: n._children)  os << *ss.get();

    return os;
  }

  /// Dump this node, in dot format.
  void logTo (std::ostream &os) const {
    os << "\t" << id() << ";\n";
    for (const Ptr &n: _children) {
      os << "\t" << id() << " -> " << n->id() << ";\n";
      n->logTo(os);
    }
  }

  /// Asserts that two phylogenetic nodes are equal
  friend void assertEqual (const Node &lhs, const Node &rhs, bool deepcopy) {
    using utils::assertEqual;

    assertEqual(bool(lhs._parent), bool(rhs._parent), deepcopy);
    if (lhs._parent && rhs._parent)
      assertEqual(lhs._parent->id(), rhs._parent->id(), deepcopy);

    assertEqual(lhs.data, rhs.data, deepcopy);
    assertEqual(lhs.contributors, rhs.contributors, deepcopy);
    assertEqual(lhs.rset, rhs.rset, deepcopy);
    assertEqual(lhs.distances, rhs.distances, deepcopy);

    assertEqual(lhs._children, rhs._children, deepcopy);
  }

private:
  /// Updates the parent with the, possibily null, species identified by \p sid
  Node* updateParent(SID sid, const Collection &nodes) {
    return _parent = (sid == SID::INVALID) ? nullptr : nodes.at(sid).get();
  }
};

} // end of namespace phylogeny

#endif // NODE_HPP
