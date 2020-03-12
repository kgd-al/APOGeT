#ifndef _PT_GRAPH_BUILDER_HPP_
#define _PT_GRAPH_BUILDER_HPP_

/*!
 * \file ptgraphbuilder.h
 *
 * Contains the declarations for the externalized phylogenic tree graphical
 * interface builder
 */

#include <QGraphicsScene>
#include <QGraphicsItem>

#include "../core/tree/phylogenetictree.hpp"

namespace gui {

/// User-controlled variations on the base phylogeny viewer
struct ViewerConfig {
  Q_GADGET
public:

  /// Minimal survival a species must have to be shown
  uint minSurvival = 0;

  /// Minimal enveloppe fullness a species must have to be shown
  float minEnveloppe = 0;

  /// Maximal range for survival monitoring
  uint clippingRange = -1;

  /// Whether only show paths leading to still alive species
  bool survivorsOnly = false;

  /// Whether to display nodes
  bool showNames = true;

  /// Whether to display hybridism graph on overlay
  bool showHybrids = false;

  /// Whether to keep the scene fully in view
  bool autofit = true;

  /// Whether to keep a screenshot per step
  bool screenshots = false;

  /// Tree rasterized radius when rendering to file
  float rasterRadius = -1;

  /// Values for the node/timelines color
  enum Colors {
    NONE = 0, ///< No color
    SURVIVORS = 1,  ///< Color survivor paths in red
    CUSTOM = 2  ///< Color path from specific species
  };
  Q_ENUM(Colors)
  Colors color = Colors::SURVIVORS; ///< Current color mode

  /// Definition of a custom coloring
  struct ColorSpec {
    phylogeny::SID sid; ///< SID of the colored species
    QColor color; ///< Color defined by the user
    bool enabled; ///< Is it currently active?
  };

  /// \cond internal
  /// Comparison structure for custom color specifications
  struct ColorSpecCMP {
    using is_transparent = void;  ///< Allows comparing with member types

    /// Compare two species id
    bool operator() (phylogeny::SID lhs, phylogeny::SID rhs) const {
      return lhs < rhs;
    }

    /// Compare two species id
    bool operator() (phylogeny::SID lhs, const ColorSpec &rhs) const {
      return operator() (lhs, rhs.sid);
    }

    /// Compare two species id
    bool operator() (const ColorSpec &lhs, phylogeny::SID rhs) const {
      return operator() (lhs.sid, rhs);
    }

    /// Compare two species id
    bool operator() (const ColorSpec &lhs, const ColorSpec &rhs) {
      return operator() (lhs.sid, rhs.sid);
    }
  };
  /// \endcond

  /// Mapping for custom coloring
  std::set<ColorSpec, ColorSpecCMP> colorSpecs;
};

/// \cond internal

struct PhylogenyViewer_base;

/// Helper alias to the species identificator
using SID = phylogeny::SID;

/// Constant pointer to a phylogenic tree viewer instance
using VTree = PhylogenyViewer_base *const;

struct GUIItems;

/// Cache structure used while building the PTree graph
struct PTreeBuildingCache {
  VTree tree;  ///< The tree being built

  const ViewerConfig &config; ///< The configuration data
  const uint time;  ///< The current timestamp

  GUIItems &items;  ///< The graphics items cache data
};

struct PolarCoordinates;
struct Path;
struct Timeline;

/// A species node
class Node : public QGraphicsItem {

  /// Whether or not this species is still alive at the tree's current timestep
  bool _alive;

  /// Whether or not this species (or one of its descendant) is still alive at
  /// the tree's current timestep
  bool _onSurvivorPath;

public:
  /// Enumeration encoding a node's visibility
  enum Visibility {
       SHOW_NAME = 1 << 0,
          PARENT = 1 << 1,
       SURVIVORS = 1 << 2,
    MIN_SURVIVAL = 1 << 3,
    MIN_FULLNESS = 1 << 4,
      CLIP_RANGE = 1 << 5
  };
  Q_DECLARE_FLAGS(Visibilities, Visibility)

  Visibilities visibilities;  ///< This node's current visibility values

  /// The tree owning this node
  PhylogenyViewer_base *const treeBase;

  const SID id;  ///< The identificator of the associated species node
  Node *parent;   ///< The parent node (if any)

  /// Helper alias for the ptree's species data
  using Data = phylogeny::SpeciesData;
  const Data &data; ///< The data of the associated species

  uint rset; ///< Size of the associated species' R-Set
  uint children;  ///< Number of subspecies

  const QString sid; ///< String representation of the node's identificator
  Path *path; ///< The graphic item connecting this node to its parent (if any)
  Timeline *timeline; ///< The graphic item depicting this node's lifetime

  /// The graphic items corresponding to the associated species 'children'
  QVector<Node*> subnodes;

  /// Current border color
  QPen coloredPen;

  /// Build a graphic node out of a potential parent and PTree data
  template <typename PN>
  Node (VTree tree, Node *parent, const PN &n)
    : treeBase(tree), id(n.id()), parent(parent), data(n.data),
      sid(QString::number(std::underlying_type<SID>::type(id))),
      path(nullptr), timeline(nullptr) {

    rset = n.rset.size();
    children = n.children().size();

    _alive = false;
    setOnSurvivorPath(false);

    autoscale();
    setAcceptHoverEvents(true);
  }

  /// Recompute all cached data
  void invalidate (const QPointF &newPos);

  /// Format species data for use in the tooltip
  QString computeTooltip (void) const;

  /// \see computeTooltip
  void updateTooltip(void) {
    setToolTip(computeTooltip());
  }

  /// Recompute personnal scale
  void autoscale (void);

  /// \returns Whether the associated species has outlived timestamp \p time
  bool isStillAlive(uint time) const {
    return data.lastAppearance >= time;
  }

  /// \returns Whether the associated species is still, currently, alive
  bool alive (void) const {
    return _alive;
  }

  /// \returns Whether the associated species or any of its descendant is still,
  /// currently, alive
  bool onSurvivorPath (void) const {
    return _onSurvivorPath;
  }

  /// Store whether the current node is still alive and notify its hierarchy
  void updateNode (bool alive);

  /// Store whether the current has any alive descendant and sets the ZValue
  /// accordingly
  void setOnSurvivorPath (bool osp);

  /// Update color based on the corresponding config values
  void updateColor (void);

  /// \returns Whether this node has sufficient visiblity values
  bool subtreeVisible (void) const {
    static constexpr auto mask =
        (SURVIVORS | MIN_FULLNESS | MIN_SURVIVAL | CLIP_RANGE | PARENT);
    return (visibilities & mask) == mask;
  }

  /// \returns Whether this node should be paint with respect to its visibility
  /// values
  bool shouldPaint (void) const {
    return visibilities | SHOW_NAME;
  }

  /// Sets visibility value \p v to \p visible
  void setVisible (Visibility v, bool visible);

  /// \returns the timestep at which this species first appeared
  uint appearance (void) const {
    return data.firstAppearance;
  }

  /// \returns the timestep at which this species last appeared
  uint disappearance (void) const {
    return data.lastAppearance;
  }

  /// \returns the number of timestep the associated species has lived for
  uint survival (void) const {
    return data.lastAppearance - data.firstAppearance;
  }

  /// \returns the ratio of the associated species' enveloppe points with
  /// respect to the target amount defined in config::PTree::enveloppeSize
  float fullness (void) const {
    return float(rset) / config::PTree::rsetSize();
  }

  /// Triggers a callback when this species node is hovered
  void hoverEnterEvent(QGraphicsSceneHoverEvent*) override;

  /// Triggers a callback when this species node is no longer hovered
  void hoverLeaveEvent(QGraphicsSceneHoverEvent*) override;

  /// Requests display of the species details
  void mouseDoubleClickEvent(QGraphicsSceneMouseEvent *e) override;

  /// Request contextual menu display
  void contextMenuEvent(QGraphicsSceneContextMenuEvent *e) override;

  /// \returns this graphics item bounding box
  QRectF boundingRect(void) const override;

  /// \returns this graphics item real shape
  QPainterPath shape(void) const override;

  /// Paints this node through the provided painter
  void paint (QPainter *painter, const QStyleOptionGraphicsItem*, QWidget*) override;
};
Q_DECLARE_OPERATORS_FOR_FLAGS(Node::Visibilities)

/// Graphics item connecting a child Node to its parent's timeline
struct Path : public QGraphicsItem {
  Node *start;  ///< The source Node (parent)
  Node *end;    ///< The target Node (child)

  /// The shape used to paint this path
  QPainterPath _shape;

  /// Create a path object between parent Node \p start and \p end
  Path(Node *start, Node *end);

  /// \copydoc Node::invalidate
  void invalidatePath (void);

  /// \copydoc Node::boundingRect
  QRectF boundingRect() const override;

  /// \returns this graphics item shape
  QPainterPath shape (void) const override {
    return _shape;
  }

  /// \copydoc Node::paint
  void paint(QPainter *painter, const QStyleOptionGraphicsItem*, QWidget*) override;
};

/// Graphics item representing a Node lifespan
struct Timeline : public QGraphicsItem {
  Node *node; ///< The parent node

  /// The point collection describing the path and survivor state
  ///  - 0: parent node position
  ///  - 1: position at which last seen alive (parent node or a descendant)
  ///  - 2: end-of-life position
  QPointF points[3];

  /// Colors for the segments points[0]-points[1] and points[1]-points[2]
  QColor colors [2];

  /// Create a timeline associated with \p node
  Timeline(Node *node);

  /// \copydoc Path::invalidatePath
  void invalidatePath (void);

  /// \copydoc Node::boundingRect
  QRectF boundingRect() const override {
    return shape().boundingRect();
  }

  /// \copydoc Path::shape
  QPainterPath shape (void) const override;

  /// \copydoc Node::paint
  void paint(QPainter *painter, const QStyleOptionGraphicsItem*, QWidget*) override;
};

/// Displays species tracking data
struct Tracker : public QGraphicsItem {  
  VTree tree; ///< The tree whose species are tracked

  /// A single tracked (or ancestor of tracked) species
  struct TrackedSpecies {
    const Node *species;  ///< The species in question
    QList<TrackedSpecies*> descendants; ///< Its (in)direct descendants

    QPainterPath path; ///< The region of direct influence
    QPainterPath hollowedPath;  ///< The region of exclusive influence

    QColor color; ///< The color of influence

    ~TrackedSpecies(void) {
      for (auto *ts: descendants) delete ts;
    }
  };

  /// The root of the tracked hierarchy tree
  TrackedSpecies *commonAncestor;

public:
  /// Builds a species tracking drawer
  Tracker (VTree tree);

  /// \returns the same bounding rect as the graph's bounds
  QRectF boundingRect(void) const override;

  /// Recomputes tracking data
  void updateTracking (void);

  /// Paints the paths for the various tracked species
  void paint (QPainter *painter, const QStyleOptionGraphicsItem*,
              QWidget*) override;

private:
  /// Recursive painter for a given tracked species
  void paint (QPainter *painter, const TrackedSpecies *ts);
};

/// Graphics item displaying a node's contributor
struct Contributors : public QGraphicsItem {
  VTree tree; ///< The tree whose contributions this draws

  /// Identificator of a path between two nodes
  struct PathID {
    QPointF from; ///< Source
    QPointF to;   ///< Destination

    /// Compare two points in lexicographic order
    static bool lower (const QPointF &lhs, const QPointF &rhs) {
      return lhs.x() != rhs.x() ? lhs.x() < rhs.x() : lhs.y() < rhs.y();
    }

    /// Compare two path identificators in lexicographic order
    friend bool operator< (const PathID &lhs, const PathID &rhs) {
      return lhs.from != rhs.from ?
        lower(lhs.from, rhs.from) : lower(lhs.to, rhs.to);
    }
  };

  /// Describes a path portion
  struct Path {
    QPainterPath path; ///< The Qt path
    float width;  ///< The path width
  };

  QMap<PathID, Path> paths;  ///< The paths connecting to the contributors
  QVector<QPair<QPointF, QString>> labels;  ///< The labels for each contributor

  Node *species;  ///< The species whose contributors are shown (or null)

  /// Builds a contributors drawer
  Contributors (VTree tree);

  /// Show the drawer for the provided node
  void show (SID sid, const GUIItems &items,
             const phylogeny::Contributors &contribs);

  /// Hide the drawer
  void hide (void);

  /// \returns the same bounding rect as the graph's bounds
  QRectF boundingRect(void) const override;

  /// Paints the paths to the various contributors
  void paint (QPainter *painter, const QStyleOptionGraphicsItem*,
              QWidget*) override;

private:
  /// \returns the identificator for the given path
  PathID pathID (const QPainterPath &p);

  /// Creates a new entry for p's id or update existing weight
  void addOrUpdate (const QPainterPath &p, float w);

  /// Register a, potentially subdivided, vertical path between the two nodes
  void verticalPath (Node *n0, Node *n1, float w);

  /// Register the horizontal path from \p n to its parent's timeline
  /// if \p vertical, also register the vertical section
  void makePath (Node *n, float w, bool vertical);
};

/// Graphics item dimming out parts of the tree
struct Dimmer : public QGraphicsItem {

  VTree tree; ///< The tree whose contributions this draws

  /// Path used to dim out parts of the tree
  QPainterPath dimPath;

  /// Builds a contributors drawer
  Dimmer (VTree tree);

  /// Sets the path used for dimming out the tree
  void setDimmingPath (const QPainterPath &path) {
    dimPath = path;
    update();
  }

  /// \returns the same bounding rect as the graph's bounds
  QRectF boundingRect(void) const override;

  /// Paints the paths to the various contributors
  void paint (QPainter *painter, const QStyleOptionGraphicsItem*, QWidget*) override;
};

/// Graphics item managing the graph's boundaries and legend
struct Border : public QGraphicsItem {
  VTree tree; ///< The tree whose border it is drawing
  bool empty; ///< Whether or not any data was inputed in the associated graph
  double radius;  ///< How many timesteps are registered in the associated ptree
  QPainterPath shape; ///< The legend axis
  QList<QPair<int, QPointF>> legend;  ///< The legend values

  QPen pen; ///< The pen used to stroke the legend axis
  QFont font; ///< The font used to type out the legend text
  QFontMetrics metrics; ///< Used to compute the size of the legend text

  /// Create border graphics item with given initial height
  Border (VTree tree, double radius);

  /// Set whether or not this item has data to display
  void setEmpty (bool empty) {
    bool wasEmpty = this->empty;
    this->empty = empty;
    if (wasEmpty != empty) {
      updateShape();
      update();
    }
  }

  /// Set current height (i.e. number of timesteps in the associated ptree)
  void setRadius (double r) {
    radius = r;
    updateShape();
    update();
  }

  /// \copydoc Node::boundingRect
  QRectF boundingRect(void) const {
    return shape.boundingRect().adjusted(0, -metrics.ascent(), 0, 0);
  }

  /// \copydoc Node::invalidate
  void updateShape (void);

  /// \copydoc Node::paint
  void paint(QPainter *painter, const QStyleOptionGraphicsItem*, QWidget*);
};

namespace details {
/// The set of pen available for painting
enum PenType {
  PATH_BASE, ///< Default pen
  PATH_SURVIVOR, ///< Paths leading to survivor species stand out
  PATH_CONTRIBUTOR, ///< Paths leading to a species' contributor
  BORDER_AXIS,  ///< Major grid axis
};
} // end of namespace details

/// Cache structure for easy management of the graph's various components
struct GUIItems {
  bool initialized; ///< Whether or not the items have been allocated

  QGraphicsScene *scene;  ///< Root of all the graphics items
  Border *border; ///< Border & legend manager
  Node *root; ///< Root of the graph's tree

  /// Species tracking drawer
  Tracker *tracker;

  /// Species contributions drawer
  Contributors *contributors;

  /// Clipping dimmer
  Dimmer *dimmer;

  QMap<SID, Node*> nodes;  ///< Lookup table for the graphics nodes

  QMap<details::PenType, QPen> pens;  ///< Collection of pens
};

/// Helper structure managing the construction of a PTree's associated graph
struct PTGraphBuilder {
  using Config = ViewerConfig;  ///< \copydoc gui::ViewerConfig

  /// Helper alias for the cache regrouping all building material
  using Cache = PTreeBuildingCache;

  /// Helper alias for the coordinate type
  using Coordinates = PolarCoordinates;

  /// Helper alias for the pen collection
  using PenSet = decltype(GUIItems::pens);

  /// Builds the pen set used for drawing stuff
  /// (precise, concise documentation is key)
  static PenSet buildPenSet (void);

  /// Updates pens to match the current state of the tree
  static void updatePenSet (float radius, PenSet &pens);

  /// \returns the appropriate width for a node drawn in a tree of radius
  static float nodeWidth (float radius);

  /// \returns the appropriate width for a path drawn in a tree of radius
  static float pathWidth (float baseWidth, float radius);

  /// \returns the appropriate width for a end-of-path decoration in a tree of
  /// radius
  static float plopRadius (float baseWidth, float radius);

  /// \returns the appropriate font size for a tree of radius
  static float fontSize (float radius);

  /// Parse the \p pt PTree and build the associated graph complete with nodes,
  /// paths and legend
  template <typename GENOME, typename UDATA>
  static void fillScene (const phylogeny::PhylogeneticTree<GENOME, UDATA> &pt,
                         Cache &cache) {

    cache.items.border = new Border(cache.tree, cache.time);
    cache.items.scene->addItem(cache.items.border);

    if (auto root = pt.root())
      addSpecies(nullptr, *root, cache);

    cache.items.border->setEmpty(!bool(pt.root()));

    cache.items.tracker = new Tracker (cache.tree);
    cache.items.scene->addItem(cache.items.tracker);

    cache.items.contributors = new Contributors(cache.tree);
    cache.items.scene->addItem(cache.items.contributors);

    cache.items.dimmer = new Dimmer(cache.tree);
    cache.items.scene->addItem(cache.items.dimmer);

    cache.items.scene->setSceneRect(cache.items.border->boundingRect());

    cache.items.initialized = true;
  }

  /// Append a new Node to the graph based on the data contained in \p n
  template <typename PN>
  static void addSpecies(Node *parent, const PN &n, Cache &cache) {
    // Create node
    Node *gn = new Node (cache.tree, parent, n);

    // Update related cache values
    if (parent)
          parent->subnodes.push_front(gn);
    else  cache.items.root = gn;
    cache.items.nodes[gn->id] = gn;
    cache.items.scene->addItem(gn);

    // Generate path to parent if needed
    if (parent) {
      auto gp = new gui::Path(parent, gn);
      gn->path = gp;
      cache.items.scene->addItem(gp);
      gp->setVisible(gn->subtreeVisible());
    }

    // Create timeline object
    auto gt = new gui::Timeline(gn);
    gn->timeline = gt;
    cache.items.scene->addItem(gt);
    gt->setVisible(gn->subtreeVisible());

    // Process subspecies
    for (const auto &n_: n.children())
      addSpecies(gn, *n_, cache);

    // Manage visibility
    gn->updateNode(gn->isStillAlive(cache.time));
    gn->setVisible(Node::SHOW_NAME, cache.config.showNames);
    gn->setVisible(Node::SURVIVORS, !cache.config.survivorsOnly || gn->onSurvivorPath());
    gn->setVisible(Node::MIN_SURVIVAL, gn->survival() >= cache.config.minSurvival);
    gn->setVisible(Node::MIN_FULLNESS, gn->fullness() >= cache.config.minEnveloppe);
    gn->setVisible(Node::CLIP_RANGE, gn->appearance() <= cache.config.clippingRange);
    gn->setVisible(Node::PARENT, parent ? parent->subtreeVisible() : true);

    // Update rendering
    gn->updateColor();
  }

  /// Recompute all graphics items positions (nodes, paths, timelines)
  static void updateLayout (GUIItems &items);

private:
  /// Recompute all graphics items position for the subtree rooted at \p localRoot
  static void updateLayout (Node *localRoot, PolarCoordinates &pc);
};

/// \endcond

} // end of namespace gui

#endif // _PT_GRAPH_BUILDER_HPP_
