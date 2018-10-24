#ifndef NODE_HPP
#define NODE_HPP

/*!
 * \file node.hpp
 *
 * Contains the definition for a single species node in the phylogenic tree
 */

#include "enumvector.hpp"
#include "speciesdata.hpp"
#include "speciescontributors.h"

namespace phylogeny {

/// Species node
template <typename GENOME>
struct Node {
  /// Helper alias to the type used for a pointer to node
  using Ptr = std::shared_ptr<Node>;

  /// Helper alias to a collection of nodes
  using Collection = enumvector<SID, Ptr>;

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
  std::vector<GENOME> enveloppe;

  /// Cache map for the intra-enveloppe distances
  _details::DistanceMap distances;

  /// Creates a node from a contributors collection (hidden from user. use the
  /// make_shared version)
  explicit Node (Contributors &&contribs, const cookie&)
    : contributors(contribs) {}

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
  Node* parent (void) { return _parent;  }

  /// \returns the collection of subspecies root at this node
  const auto& children (void) const {
    return _children;
  }

  /// \returns the subspecies at index \p i
  auto& child (size_t i) {
    return _children[i];
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
  Node* update (Contributors::Contribution sids, const Collection &nodes) {
    SID mainSID = contributors.update(sids, elligibilityTester(nodes));
    return _parent = (mainSID == SID::INVALID) ? nullptr : nodes[mainSID].get();
  }

  /// Triggers a tree-wide recomputation of the elligibilities of all node
  /// contributors. Possible chain-reaction (multiple parent modifications)
  void updateElligibilities (const Collection &nodes) {
    assert(false);
  }

  /// Stream this node. Mostly for debugging purpose.
  friend std::ostream& operator<< (std::ostream &os, const Node &n) {
    std::string spacing = "> ";
    const Node *p = &n;
    while ((p = p->_parent))  spacing += "  ";

    os << spacing << "[" << n.id << "] ( ";
    for (const GENOME &g: n.enveloppe)    os << g.id() << " ";
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
};

} // end of namespace phylogeny

#endif // NODE_HPP
