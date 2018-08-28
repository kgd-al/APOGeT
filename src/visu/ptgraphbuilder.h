#ifndef _PTREE_INTROSPECTER_HPP_
#define _PTREE_INTROSPECTER_HPP_

#include <QGraphicsScene>
#include <QGraphicsItem>

#include <QTextStream>

#include "../core/phylogenictree.hpp"

namespace gui {

struct ViewerConfig {
  uint minSurvival;
  float minEnveloppe;
  bool showNames;
  bool circular;
  bool autofit;
};

struct PTreeBuildingCache {
  const ViewerConfig &config;
  const uint time;

  QGraphicsScene *scene;

  struct PolarCoordinates {
    const double originalWidth, widthWithLegend;

    ///< Also inverted to cope with Qt coordinate system
    const double phase;

    uint nextX;

    PolarCoordinates (double width);
    static float xCoord (uint i);
    double toPrimaryAngle (const QPointF &p) const;
    QPointF operator() (uint time);
  } polarCoordinate;
};

class Node : public QGraphicsItem {
  bool _stillAlive;

public:
  const uint id;

  using Data = phylogeny::SpeciesData;
  const Data &data;

  float fullness;
  uint children;

  template <typename PN>
  Node (const QPointF pos, const PN &n) : id(n.id), data(n.data) {
    fullness = n.enveloppe.size() / PTreeConfig::enveloppeSize();
    children = n.children.size();

    setPos(pos);
    updateTooltip();
  }

  void updateTooltip(void);

  bool isStillAlive(uint time) const {
    return data.lastAppearance >= time;
  }

  bool isStillAlive (void) const {
    return _stillAlive;
  }

  void setStillAlive (bool alive) {
    _stillAlive = alive;
    setZValue(_stillAlive ? 2 : 1);
  }

  uint survival (void) const {
    return data.lastAppearance - data.firstAppearance;
  }

  QRectF boundingRect(void) const;

  void paint (QPainter *painter, const QStyleOptionGraphicsItem*, QWidget*);
};

struct Path : public QGraphicsItem {
  using Cache = PTreeBuildingCache;

  Node *_start, *_end;
  QPainterPath _shape, _survivor;
  QPen _pen;

  Path(Node *start, Node *end, const Cache &cache);

  QRectF boundingRect() const;
  double length (const QPointF &p) {
    return sqrt(p.x()*p.x() + p.y()*p.y());
  }

  QPainterPath shape (void) const override {
    return _shape.united(_survivor);
  }

  void paint(QPainter *painter, const QStyleOptionGraphicsItem *, QWidget *);
};

struct Border : public QGraphicsItem {
  double height;
  QPainterPath shape;
  QList<QPair<int, QPointF>> legend;

  QPen pen;
  QFont font;
  QFontMetrics metrics;

  Border (double height);

  QRectF boundingRect(void) const {
    return shape.boundingRect().adjusted(0, -metrics.ascent(), 0, 0);
  }

  void paint(QPainter *painter, const QStyleOptionGraphicsItem*, QWidget*);
};

} // end of namespace gui

namespace phylogeny {

template <typename GENOME>
struct PTreeIntrospecter {
  using PT = PhylogenicTree<GENOME>;
  using PN = typename PT::Node;

  using Config = gui::ViewerConfig;
  using Cache = gui::PTreeBuildingCache;


  static void fillScene (const PT &pt, QGraphicsScene *scene, const Config &config) {
    uint step = pt.step();
    double width = Cache::PolarCoordinates::xCoord(pt.width()-2);

    Cache c {
      config, step,
      scene,
      {width}
    };

    if (pt._root)
      addSpecies(nullptr, *pt._root, c);

    auto border = new gui::Border(step);
    scene->addItem(border);

    scene->setSceneRect(border->boundingRect());
  }

  static bool addSpecies(gui::Node *parent, const PN &n, Cache &cache) {
    QPointF pos;
    if (parent)
      pos += cache.polarCoordinate(n.data.firstAppearance);

    gui::Node *gn = new gui::Node (pos, n);
    uint survival = gn->survival();
    float fullness = gn->fullness();
    bool survivor = gn->isStillAlive(cache.time);

    gn->setScale(fullness);

    bool visible = cache.config.showNames;
    visible &= (survival >= cache.config.minSurvival);
    visible &= (fullness >= cache.config.minEnveloppe);
    gn->setVisible(visible);

    cache.scene->addItem(gn);

    for (auto it=n.children.rbegin(); it!=n.children.rend(); ++it)
      survivor |= addSpecies(gn, *(*it), cache);

    gn->setStillAlive(survivor);
    cache.scene->addItem(new gui::Path(parent, gn, cache));

    return survivor;
  }

  template <typename F>
  static void update (QGraphicsScene *scene, F f) {
    for (auto *i: scene->items()) f(i);
  }
};

} // end of namespace phylogeny

#endif // _PTREE_INTROSPECTER_HPP_
