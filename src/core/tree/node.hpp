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

  /// Prevents constructor access from outside the class
  struct cookie {};

  SpeciesData data; ///< Species additionnal data

  /// Collection of contributors the this species' gene pool
  Contributors contributors;

  /// Reference to the species' parent (main contributor)
  Node *parent;

  /// Subspecies of this node
  std::vector<Ptr> children;

  /// Collection of borderoids (in opposition to centroids)
  std::vector<GENOME> enveloppe;

  /// Cache map for the intra-enveloppe distances
  _details::DistanceMap distances;

  /// Creates a node from a contributors collection (hidden from user. use the
  /// make_shared version)
  explicit Node (Contributors &&contribs, const Collection &nodes, const cookie&)
    : contributors(contribs) {
    update({}, nodes);
  }

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
  Node* getParent (void) { return parent;  }

  /// Updates the species contributions manager and the species' main parent
  /// \returns the new species' main parent
  Node* update (Contributors::Contribution sids, const Collection &nodes) {
    SID mainSID = contributors.update(sids);
    if (mainSID == SID::INVALID)
          parent = nullptr;
    else  parent = nodes[mainSID].get();
    return parent;
  }

  /// Stream this node. Mostly for debugging purpose.
  friend std::ostream& operator<< (std::ostream &os, const Node &n) {
    std::string spacing = "> ";
    const Node *p = &n;
    while ((p = p->parent))  spacing += "  ";

    os << spacing << "[" << n.id << "] ( ";
    for (const GENOME &g: n.enveloppe)    os << g.id() << " ";
    os << ")\n";

    for (const Ptr &ss: n.children)  os << *ss.get();

    return os;
  }

  /// Dump this node, in dot format.
  void logTo (std::ostream &os) const {
    os << "\t" << id << ";\n";
    for (const Ptr &n: children) {
      os << "\t" << id << " -> " << n->id << ";\n";
      n->logTo(os);
    }
  }
};

} // end of namespace phylogeny

#endif // NODE_HPP
