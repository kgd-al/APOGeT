#ifndef _PT_GRAPH_BUILDER_HPP_
#define _PT_GRAPH_BUILDER_HPP_

#include <QGraphicsScene>
#include <QGraphicsItem>

#include <QTextStream>

#include "../core/phylogenictree.hpp"


#include <QDebug>


namespace gui {

struct ViewerConfig {
  uint minSurvival;
  float minEnveloppe;
  bool showNames;
  bool circular;
  bool autofit;
};

struct GUIItems;
struct PTreeBuildingCache {
  const ViewerConfig &config;
  const uint time;

  GUIItems &items;
};

struct PolarCoordinates;

struct Path;
class Node : public QGraphicsItem {
  bool _alive;
  bool _onSurvivorPath;

public:
  enum Visibility {
    SHOW_NAME = 1,
    PARENT = 2,
    MIN_SURVIVAL = 4,
    MIN_FULLNESS = 8
  };
  Q_DECLARE_FLAGS(Visibilities, Visibility)
  Visibilities visibilities;

  const uint id;
  Node *parent;

  using Data = phylogeny::SpeciesData;
  const Data &data;

  uint enveloppe;
  uint children;

  Path *path;
  QVector<Node*> subnodes;

  template <typename PN>
  Node (Node *parent, const PN &n) : id(n.id), parent(parent), data(n.data) {
    enveloppe = n.enveloppe.size();
    children = n.children.size();

    _alive = false;
    _onSurvivorPath = false;

    autoscale();
  }

  void updateTooltip(void);
  void autoscale (void);

  bool isStillAlive(uint time) const {
    return data.lastAppearance >= time;
  }

  bool alive (void) const {
    return _alive;
  }

  bool onSurvivorPath (void) const {
    return _onSurvivorPath;
  }

  void updateNode (bool alive, uint time);

  void setOnSurvivorPath (bool osp);

  bool subtreeVisible (void) const {
    static constexpr auto mask = (MIN_FULLNESS | MIN_SURVIVAL | PARENT);
    return (visibilities & mask) == mask;
  }

  bool shouldPaint (void) const {
    return visibilities | SHOW_NAME;
  }

  void setVisible (Visibility v, bool visible);

  uint survival (void) const {
    return data.lastAppearance - data.firstAppearance;
  }

  float fullness (void) const {
    return float(enveloppe) / PTreeConfig::enveloppeSize();
  }

  QRectF boundingRect(void) const;

  void paint (QPainter *painter, const QStyleOptionGraphicsItem*, QWidget*);
};
Q_DECLARE_OPERATORS_FOR_FLAGS(Node::Visibilities)

struct Path : public QGraphicsItem {
  using Cache = PTreeBuildingCache;

  Node *_start, *_end;
  QPainterPath _shape, _survivor;
  QPen _pen;

  Path(Node *start, Node *end);

  void updatePath (void);
  void updatePath (uint time);

  QRectF boundingRect() const;

  QPainterPath shape (void) const override {
    return _shape.united(_survivor);
  }

  void paint(QPainter *painter, const QStyleOptionGraphicsItem *, QWidget *);
};

struct Border : public QGraphicsItem {
  bool empty;
  double height;
  QPainterPath shape;
  QList<QPair<int, QPointF>> legend;

  QPen pen;
  QFont font;
  QFontMetrics metrics;

  Border (double height);

  void setEmpty (bool empty) {
    bool wasEmpty = this->empty;
    this->empty = empty;
    if (wasEmpty != empty) {
      updateShape();
      update();
    }
  }

  void setHeight (double h) {
    height = h;
    updateShape();
    update();
  }

  QRectF boundingRect(void) const {
    return shape.boundingRect().adjusted(0, -metrics.ascent(), 0, 0);
  }

  void updateShape (void);

  void paint(QPainter *painter, const QStyleOptionGraphicsItem*, QWidget*);
};

struct GUIItems {
  QGraphicsScene *scene;
  Border *border;
  Node *root;

  QMap<uint, Node*> nodes;
};

struct PTGraphBuilder {
  using Config = ViewerConfig;
  using Cache = PTreeBuildingCache;
  using Coordinates = PolarCoordinates;

  template <typename GENOME>
  static void fillScene (const phylogeny::PhylogenicTree<GENOME> &pt, Cache &cache) {
    if (auto root = pt.root())
      addSpecies(nullptr, *root, cache);

    cache.items.border = new Border(cache.time);
    cache.items.scene->addItem(cache.items.border);
    cache.items.border->setEmpty(bool(pt.root()));

    cache.items.scene->setSceneRect(cache.items.border->boundingRect());
  }

  template <typename PN>
  static bool addSpecies(Node *parent, const PN &n, Cache &cache) {
    Node *gn = new Node (parent, n);
    uint survival = gn->survival();
    float fullness = gn->fullness();
    bool survivor = gn->isStillAlive(cache.time);

    if (parent)
          parent->subnodes.push_front(gn);
    else  cache.items.root = gn;
    cache.items.nodes[gn->id] = gn;
    cache.items.scene->addItem(gn);

    for (auto it=n.children.rbegin(); it!=n.children.rend(); ++it)
      survivor |= addSpecies(gn, *(*it), cache);

    auto gp = new gui::Path(parent, gn);
    gn->path = gp;
    cache.items.scene->addItem(gp);
    gp->setVisible(gn->subtreeVisible());

    gn->setVisible(Node::SHOW_NAME, cache.config.showNames);
    gn->setVisible(Node::MIN_SURVIVAL, survival >= cache.config.minSurvival);
    gn->setVisible(Node::MIN_FULLNESS, fullness >= cache.config.minEnveloppe);
    gn->setVisible(Node::PARENT, parent ? parent->subtreeVisible() : true);
    gn->updateNode(survivor, cache.time);

    return survivor;
  }

  static void updateLayout (GUIItems &items);

private:
  static void updateLayout (Node *localRoot, PolarCoordinates &pc);
};

} // end of namespace gui

#endif // _PT_GRAPH_BUILDER_HPP_
