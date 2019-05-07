#include <QTextStream>

#include <QPainter>
#include <QToolTip>

#include <QVector3D>
#include <qmath.h>

/// \todo remove
#include <QDebug>

#include "phylogenyviewer.h"
#include "ptgraphbuilder.h"
#include "graphicutils.h"

namespace gui {

// ============================================================================
// == Constants
// ============================================================================

static constexpr bool debugDrawAABB = false;

// == Legend ========================================================

///< Note that the coordinate system is inverted (y points downward)
static constexpr float LEGEND_PHASE = -M_PI / 2;
static constexpr float LEGEND_SPACE = M_PI / 12;
static constexpr uint LEGEND_TICKS = 4;

// == Nodes style ===================================================

static constexpr float NODE_RADIUS = 10;
static constexpr float NODE_MARGIN = 2;
static constexpr float NODE_SIZE = 2 * (NODE_RADIUS + NODE_MARGIN);

static constexpr float END_POINT_SIZE = NODE_RADIUS / 4;

// == Z-values ======================================================

static constexpr int NODE_SURVIVOR_LEVEL = 11;
static constexpr int NODE_EXTINCT_LEVEL = 10;

static constexpr int DIMMER_LEVEL = 0;
static constexpr int CONTRIBUTORS_LEVEL = 0;

static constexpr int PATH_SURVIVOR_LEVEL = -5;
static constexpr int TIMELINE_SURVIVOR_LEVEL = -6;

static constexpr int PATH_EXTINCT_LEVEL = -10;
static constexpr int TIMELINE_EXTINCT_LEVEL = -11;

static constexpr int STRACKING_LEVEL = -20;

static constexpr int BOUNDS_LEVEL = -30;

// == Paint style ===================================================

static constexpr float AXIS_WIDTH = 1;
static constexpr float PATH_WIDTH = 1.5;

static constexpr Qt::GlobalColor PATH_DEFAULT_COLOR = Qt::darkGray;
static constexpr Qt::GlobalColor PATH_SURVIVOR_COLOR = Qt::red;
static constexpr Qt::GlobalColor PATH_CONTRIBUTOR_COLOR = Qt::green;


// ============================================================================
// == Coordinate computation
// ============================================================================

struct PolarF {
  double a, r;
};

double angle (const QPointF &p) {
  return std::atan2(p.y(), p.x());
}

double radius (const QPointF &p) {
  return std::sqrt(p.x()*p.x() + p.y()*p.y());
}

PolarF toPolar (const QPointF &p) {
  return PolarF { angle(p), radius(p) };
}

QPointF toCartesian (double a, double r) {
  return r * QPointF(cos(a), sin(a));
}

QPointF toCartesian (const PolarF &p) {
  return toCartesian(p.a, p.r);
}

/// Generates and manages polar coordinates for the nodes/paths
struct PolarCoordinates {

  /// The angular phase used in coordinate computation.
  /// Inverted to cope with Qt coordinate system.
  static constexpr float phase = LEGEND_PHASE + LEGEND_SPACE/2.;

  const double width; ///< Width of the graph (legend included)

  uint nextX; ///< X coordinate of the next node

  /// \returns The angle for \p in the range [phase,2&pi;+phase]
  static double primaryAngle (double a) {
    while (a < phase) a += 2 * M_PI;
    while (a > 2 * M_PI + phase) a -= 2 * M_PI;
    return a;
  }

  /// \returns The angle for \p in the range [phase,2&pi;+phase]
  static double primaryAngle (const QPointF &p) {
    if (p.isNull()) return phase;
    return primaryAngle(atan2(p.y(), p.x()));
  }

  /// \returns The euclidian distance of \p p with the origin
  static double length (const QPointF &p) {
    return sqrt(p.x()*p.x() + p.y()*p.y());
  }

  /// \returns the coordinate of node number \p i
  static float xCoord (uint i) {  return i * NODE_SIZE; }

  /// Creates a polar coordinates object with the specified unscaled \p width
  PolarCoordinates (double width)
    : width(2 * M_PI * width / (2 * M_PI - LEGEND_SPACE)),
      nextX(0) {}

  /// \returns the position of the next point
  QPointF operator() (uint time) {
    double a = phase;
    if (width > 0)  a += 2. * M_PI * xCoord(nextX++) / width;
  //  qDebug() << "theta(" << nextX-1 << "): " << a * 180. / M_PI;
    return toCartesian(a, time);
  }
};


// ============================================================================
// == Utilities
// ============================================================================

QPainterPath& addArc (QPainterPath &p, const QPointF &p1, int sign = 1) {
  double a0 = PolarCoordinates::primaryAngle(p.currentPosition());
  double a1 = PolarCoordinates::primaryAngle(p1);
  double r1 = PolarCoordinates::length(p1);

  p.arcTo(QRect(-r1, -r1, 2*r1, 2*r1), -qRadiansToDegrees(a0),
          qRadiansToDegrees(sign*(a0 - a1)));

  return p;
}

QPainterPath makeArc (const QPointF &p0, const QPointF &p1) {
  QPainterPath path;
  path.moveTo(p0);
  return addArc(path, p1);
}

QPointF timelineAnchor (const Node *n) {
  return toCartesian(
    PolarCoordinates::primaryAngle(n->parent->scenePos()),
    PolarCoordinates::length(n->scenePos()));
}


// ============================================================================
// == Graph node
// ============================================================================

QRectF Node::boundingRect(void) const {
  return QRectF(-.5*NODE_SIZE, -.5*NODE_SIZE, NODE_SIZE, NODE_SIZE);
}

void Node::invalidate(const QPointF &newPos) {
  setPos(newPos);
  if (path) path->invalidatePath();
  timeline->invalidatePath();
  update();
}

QString Node::computeTooltip (void) const {
  QString tooltip;
  QTextStream qss (&tooltip);
  qss << "Node " << std::underlying_type<SID>::type(id) << "\n"
      << "Enveloppe: " << 100 * fullness() << "%\n"
      << "Appeared at " << data.firstAppearance << "\n"
      << "Disappeared at "
        << (_alive ? "-" : QString::number(data.lastAppearance))
        << "\n"
      << data.count << " individuals\n"
      << children << " subspecies";
  return tooltip;
}

void Node::autoscale(void) {
  setScale(fullness() * PTGraphBuilder::nodeWidth(treeBase->radius()));
  updateTooltip();
  update();
}

void Node::setVisible (Visibility v, bool visible) {
  visibilities.setFlag(v, visible);

  if (v != SHOW_NAME) {
    // Update own visibility as well as related paths'
    visible = subtreeVisible();
    if (path) path->setVisible(visible);
    timeline->setVisible(visible);
    QGraphicsItem::setVisible(visible);

    // Propagate to children
    for (Node *n: subnodes)
      n->setVisible(PARENT, visible);
  }
}

void Node::updateNode (bool alive) {
  _alive = alive;

  // Notify hierarchy
  if (_alive) {
    Node *n = this;
    do {
      n->setOnSurvivorPath(true);
      n = n->parent;
    } while (n && !n->_onSurvivorPath);
  }

  if (path) path->invalidatePath();
  timeline->invalidatePath();
  updateTooltip();
  autoscale();
}

void Node::setOnSurvivorPath (bool osp) {
  static constexpr int levels [3][2] = {
    { NODE_EXTINCT_LEVEL, NODE_SURVIVOR_LEVEL },
    { TIMELINE_EXTINCT_LEVEL, TIMELINE_SURVIVOR_LEVEL },
    { PATH_EXTINCT_LEVEL, PATH_SURVIVOR_LEVEL },
  };
  bool s = _onSurvivorPath = osp;
  setZValue(levels[0][s]);
  if (timeline) timeline->setZValue(levels[1][s]);
  if (path) path->setZValue(levels[2][s]);
}

void Node::hoverEnterEvent(QGraphicsSceneHoverEvent*) {
  treeBase->hoverEvent(id, true);
}

/// Triggers a callback when this species node is no longer hovered
void Node::hoverLeaveEvent(QGraphicsSceneHoverEvent*) {
  treeBase->hoverEvent(id, false);
}

void Node::mouseDoubleClickEvent(QGraphicsSceneMouseEvent *e) {
  treeBase->doubleClickEvent(*this, e);
}

void Node::updateColor(void) {
  const auto &config = treeBase->config();
  coloredPen.setColor(PATH_DEFAULT_COLOR);
  if (config.color == ViewerConfig::SURVIVORS && _onSurvivorPath)
    coloredPen.setColor(PATH_SURVIVOR_COLOR);
  else if (config.color == ViewerConfig::CUSTOM) {
    /// TODO could be improved. See GUIItems::nodes
    auto it = config.colorSpecs.find(id);
    if (it != config.colorSpecs.end())
      coloredPen.setColor(it->color);
  }
  update();
  timeline->update();
  if (path) path->update();
}

void Node::paint (QPainter *painter, const QStyleOptionGraphicsItem*, QWidget*) {
  if (debugDrawAABB) {
    painter->save();
    QPen pen = painter->pen();
    pen.setColor(Qt::blue);
    pen.setWidthF(0);
    painter->setPen(pen);
    painter->drawRect(boundingRect());
    painter->restore();
  }

  if (visibilities.testFlag(SHOW_NAME)) {    
    float pw = painter->pen().widthF();
    QRectF r (boundingRect().center() - QPoint(NODE_RADIUS, NODE_RADIUS),
              2*QSizeF(NODE_RADIUS, NODE_RADIUS));
    r.adjust(pw, pw, -pw, -pw);

    painter->save();
      QFont f = painter->font();
      f.setPixelSize(12);
      painter->setFont(f);

      painter->setClipRect(boundingRect());

      // Test to see correct size
      QRectF bounds;
      painter->setPen(Qt::transparent);
      painter->drawText(r, Qt::AlignCenter, sid, &bounds);

      painter->setBrush(Qt::white);
      painter->setPen(coloredPen);
      painter->drawEllipse(boundingRect().center(), NODE_RADIUS, NODE_RADIUS);

      // Scale and do paint
      float s = r.width() / bounds.width();
      painter->setPen(Qt::black);
      if (s < 1) {
        painter->scale(s, s);
        r = QTransform::fromScale(1 / s, 1 / s).mapRect(r);
      }
      painter->drawText(r, Qt::AlignCenter, sid);

    painter->restore();
  }
}


// ============================================================================
// == Path between a parent and child node
// ============================================================================

Path::Path(Node *start, Node *end) : start(start), end(end) {
  setZValue(PATH_EXTINCT_LEVEL);
  assert(start && end);
  assert(start->treeBase == end->treeBase);
}

void Path::invalidatePath(void) {
  prepareGeometryChange();

  _shape = QPainterPath();
  _shape.setFillRule(Qt::WindingFill);

  _shape.addPath(makeArc(start->scenePos(), end->scenePos()));
  _shape.addEllipse(_shape.pointAtPercent(1), END_POINT_SIZE, END_POINT_SIZE);

  update();
}

QRectF Path::boundingRect() const {
  qreal extra = (PATH_WIDTH + 20) / 2.0;

  return shape().boundingRect().normalized()
      .adjusted(-extra, -extra, extra, extra);
}

void Path::paint(QPainter *painter, const QStyleOptionGraphicsItem*, QWidget*) {
  if (debugDrawAABB) {
    painter->save();
    QPen pen = painter->pen();
    pen.setColor(Qt::red);
    pen.setWidthF(0);
    painter->setPen(pen);
    painter->drawRect(boundingRect());
    painter->restore();
  }

  QPen pen = end->treeBase->pathPen(details::PATH_BASE);
  pen.setColor(end->coloredPen.color());
  painter->setPen(pen);
  painter->drawPath(_shape);
}


// ============================================================================
// == Timeline for a node
// ============================================================================

Timeline::Timeline(Node *node) : node(node) {
  setZValue(TIMELINE_EXTINCT_LEVEL);
}

void Timeline::invalidatePath(void) {
  prepareGeometryChange();

  points[0] = node->scenePos();
  double a = PolarCoordinates::primaryAngle(points[0]);

  points[2] = node->data.lastAppearance * QPointF(cos(a), sin(a));

  if (node->alive()) // Survivor timeline goes all the way
    points[1] = points[2];

  else if (!node->onSurvivorPath())  // No survivor part for this path
    points[1] = points[0];

  else {
    // Compute location of last child survivor
    double l = node->data.firstAppearance;
    for (const Node *gn: node->subnodes) {
      if (!gn->onSurvivorPath())  continue;

      double l_ = gn->data.firstAppearance;
      if (l < l_) l = l_;
    }

    points[1] = l * QPointF(cos(a), sin(a));
  }

  update();
}

QPainterPath Timeline::shape (void) const {
  QPainterPath path;
  path.moveTo(points[0]);
  path.lineTo(points[2]);
  path.addEllipse(points[2], END_POINT_SIZE, END_POINT_SIZE);
  return QPainterPathStroker (node->treeBase->pathPen(details::PATH_BASE)).createStroke(path);
//  return QPainterPath();
}

void Timeline::paint(QPainter *painter, const QStyleOptionGraphicsItem *, QWidget *) {
  if (debugDrawAABB) {
    painter->save();
    QPen pen = painter->pen();
    pen.setColor(Qt::green);
    pen.setWidthF(0);
    painter->setPen(pen);
    painter->drawRect(boundingRect());
    painter->restore();
  }

  if (points[0] != points[1]) {
    QPen pen = node->treeBase->pathPen(details::PATH_BASE);
    pen.setColor(node->coloredPen.color());
    painter->setPen(pen);
    painter->drawLine(points[0], points[1]);
  }

  if (points[1] != points[2]) {
    painter->setPen(node->treeBase->pathPen(details::PATH_BASE));
    painter->drawLine(points[1], points[2]);
  }

  painter->setBrush(painter->pen().color());
  painter->drawEllipse(points[2], END_POINT_SIZE, END_POINT_SIZE);
}


// ============================================================================
// == Species tracking drawer
// ============================================================================

Tracker::Tracker (VTree tree) : tree(tree), commonAncestor(nullptr) {
  setZValue(STRACKING_LEVEL);
}

QRectF Tracker::boundingRect(void) const {
  return tree->boundingRect();
}

struct Span {
  double a;
  uint r;

  Span (double angle, uint radius)
    : a(PolarCoordinates::primaryAngle(angle)), r(radius) {}

  Span (const Node *n) : Span(toPolar(n->pos()).a, n->disappearance()) {}

  static Span extract (const Node *n) {
    Span s (n);
    for (const Node *c: n->subnodes)
      if (c->isVisible())
        s = max(s, extract(c));
    return s;
  }

  friend Span max (const Span &lhs, const Span &rhs) {
    return Span(std::max(lhs.a, rhs.a), std::max(lhs.r, rhs.r) );
  }
};

QVector3D toV3D (const QColor &c) {
  return QVector3D(c.redF(), c.greenF(), c.blueF());
}

QColor toColor (const QVector3D &v) {
  return QColor::fromRgbF(v[0], v[1], v[2]);
}

QPainterPath buildPath (const Node *n) {
  QPainterPath path;

  QPointF startC = n->pos();
  path.moveTo(startC);

  Span span = Span::extract(n);
  path.lineTo(toCartesian(angle(n->timeline->points[2]), span.r));
  addArc(path, toCartesian(span.a, span.r), 1);
  path.lineTo(toCartesian(angle(path.currentPosition()), radius(startC)));
  addArc(path, startC);

  return path;
}

struct AncestryNode {
  const Node *node;

  AncestryNode *parent = nullptr;
  std::set<AncestryNode*> children;

  bool monitored = false;
};

template<typename SN, typename AN>
void buildAncestries (SN &specsNodes, AN &anodes,
                      AncestryNode **root, const Node *n) {

  if (n->parent)  buildAncestries(specsNodes, anodes, root, n->parent);
  if (anodes.find(n) == anodes.end()) {
    AncestryNode *an = new AncestryNode;
    an->node = n;

    if (n->parent) {
      an->parent = anodes.value(n->parent, nullptr);
      an->parent->children.insert(an);

    } else
      *root = an;

    auto it = specsNodes.find(n);
    if (it != specsNodes.end()) {
      an->monitored = true;
      specsNodes.erase(it);
    }

    anodes[n] = an;
  }
}

/// TODO Remove
void debugPrintAncestry (const AncestryNode *n, uint depth) {
  {
    auto q = qDebug().noquote().nospace();
    q << QString(depth*2, ' ');
    if (n->monitored)  q << "[";
    q << n->node->sid;
    if (n->monitored)  q << "]";
  }
  for (const AncestryNode *c: n->children)  debugPrintAncestry(c, depth+1);
}

AncestryNode* simplify (AncestryNode *n) {
  if (!n->monitored && n->children.size() == 1) {
    AncestryNode *c = *n->children.begin();
    c->parent = n->parent;
    if (n->parent)  n->parent->children.insert(c);
    return simplify(c);
  } else {
    decltype(n->children) children;
    for (AncestryNode *c: n->children)
      children.insert(simplify(c));
    n->children = children;
    return n;
  }
}

Tracker::TrackedSpecies* buildRenderingTree (const AncestryNode *n,
                                             const ViewerConfig &config) {
  const auto &specs = config.colorSpecs;

  auto *ts = new Tracker::TrackedSpecies;
  ts->species = n->node;

  ts->path = buildPath(ts->species);

  QVector3D color;
  for (const AncestryNode *c: n->children) {
    auto d = buildRenderingTree(c, config);
    ts->descendants.append(d);
    color += toV3D(d->color);
  }

  auto it = specs.find(ts->species->id);
  if (it != specs.end())  ts->color = it->color;
  else if (!ts->descendants.empty()) {
    color /= ts->descendants.size();
    qDebug() << "Color for " << ts->species->sid << ": " << color;
    ts->color = toColor(color);
  }

  return ts;
}

void Tracker::updateTracking(void) {
  const auto &config = tree->config();
  if (config.color == ViewerConfig::CUSTOM) {
    const auto &specs = config.colorSpecs;
    const auto &nodes = tree->items().nodes;

    // Erase tracking data
    if (commonAncestor) delete commonAncestor;

    // Find common ancestries
    std::set<const Node*> specsNodes;
    for (const auto &spec: specs)
      if (spec.enabled)
        specsNodes.insert(nodes.value(spec.sid));

    AncestryNode *root = nullptr;
    QMap<const Node*, AncestryNode*> anodes;
    while (!specsNodes.empty())
      buildAncestries(specsNodes, anodes, &root, *specsNodes.begin());

    qDebug() << "Obtained ancestry:";
    debugPrintAncestry(root, 0);

    root = simplify(root);

    qDebug() << "Simplified ancestry:";
    debugPrintAncestry(root, 0);

    commonAncestor = buildRenderingTree(root, config);
  }
}

void Tracker::paint (QPainter *painter, const TrackedSpecies *ts) {
  painter->save();
    QPen pen = tree->pathPen(details::PATH_BASE);
    pen.setWidthF(.5 * pen.widthF());
    pen.setColor(ts->color);
    painter->setPen(pen);
    QColor fillColor = ts->color;
    fillColor.setAlphaF(.25);
    painter->setBrush(fillColor);

    painter->drawPath(ts->path);
  painter->restore();

  for (const auto *ts_: ts->descendants)  paint(painter, ts_);
}

void Tracker::paint (QPainter *painter,
                     const QStyleOptionGraphicsItem*, QWidget*) {
  if (!commonAncestor)  return;
  painter->save();
//    painter->setCompositionMode(QPainter::CompositionMode_Source);
    paint(painter, commonAncestor);
  painter->restore();
}

// ============================================================================
// == Species contributions drawer
// ============================================================================

Contributors::Contributors (VTree tree) : tree(tree), species(nullptr) {
  setZValue(CONTRIBUTORS_LEVEL);
}

QRectF Contributors::boundingRect(void) const {
  return tree->boundingRect();
}

Contributors::PathID Contributors::pathID (const QPainterPath &p) {
  return PathID{p.pointAtPercent(0), p.pointAtPercent(1)};
}

void Contributors::addOrUpdate (const QPainterPath &p, float w) {
  if (p.pointAtPercent(0) == p.pointAtPercent(1))
    return;

  auto id = pathID(p);
  auto pit = paths.find(id);
  if (pit == paths.end())
    paths.insert(id, {p, w});
  else
    pit.value().width += w;
}

void Contributors::verticalPath(Node *n0, Node *n1, float w) {
  assert(!n1 || n0->parent == n1 || n0 == n1->parent || n0->parent == n1->parent);

  Node *parent = nullptr;
  if (!n1) parent = n0->parent;
  else if (n0->parent == n1) parent = n1;
  else if (n0 == n1->parent) parent = n0;
  else  parent = n0->parent;

  using T = QPair<uint,QPointF>;
  const auto &nodes = parent->subnodes;

  int ni0 = nodes.size()-1;
  if (parent == n0->parent) ni0 = nodes.indexOf(n0);
  assert(ni0 >= 0);

  int ni1 = nodes.size()-1;
  if (n1 && parent == n1->parent) ni1 = nodes.indexOf(n1);
  assert(ni1 >= 0);

  QVector<T> points;
  int m = std::min(ni0, ni1), M = std::max(ni0, ni1);
  for (int ni = m; ni <= M; ni++)
    points.append(T{uint(nodes[ni]->id), timelineAnchor(nodes[ni])});

  if (!n1 || n0->parent != n1->parent)
    points.append(T{uint(parent->id), parent->scenePos()});

  for (int j=0; j<points.size()-1; j++) {
    QPainterPath path;
    path.moveTo(points[j].second);
    path.lineTo(points[j+1].second);

    addOrUpdate(path, w);
  }
}

void Contributors::makePath (Node *n, float w, bool vertical) {
  QPainterPath path;
  path.addPath(makeArc(n->parent->scenePos(), n->scenePos()));
  addOrUpdate(path, w);

  if (vertical)
    verticalPath(n, nullptr, w);
}

void Contributors::show (SID sid, const GUIItems &items,
                         const phylogeny::Contributors &contribs) {

  const float R = tree->radius();

  // Retrieve graphic item of given species
  Node *n = species = items.nodes.value(sid);
  assert(n->id == sid);

  // Retrieve all of its ancestors
  std::vector<Node*> np;
  {
    Node *p = n;
    while (p) np.push_back(p), p = p->parent;
  }

  // Clean previous drawing
  paths.clear();
  labels.clear();

  // Total contribution count (excluding itself)
  float totalWidth = 0;
  for (auto &c: contribs)
    if (c.speciesID() != sid)  totalWidth += c.count();

  uint unaccounted = 0;

  // Start parsing individual contributions
  for (auto &c: contribs) {
    if (c.speciesID() == sid) continue;

    float w = c.count() / totalWidth;
    Node *nc = items.nodes.value(c.speciesID());
    if (nc) {
      // Store label
      QPoint labelPos = nc->scenePos().toPoint();
      labelPos.setX(labelPos.x() + .025*R);
      QString label = QString::number(100 * w, 'f', 2) + "%";
      labels.append(QPair<QPointF, QString>(labelPos, label));

      Node *n_ = nullptr; // iterator

      // Find path from contributor to common ancestor
      n_ = nc;
      Node *ca = nc;
      auto it = std::find(np.begin(), np.end(), n_);
      while (it == np.end()) {  // Common ancestor not found
        it = std::find(np.begin(), np.end(), n_->parent);
        makePath(n_, w, it == np.end());

        ca = n_;
        n_ = n_->parent;
        assert(n_);
      }
      assert(ca);

      Node *commonAncestor = n_;

      // Find path from node to common ancestor
      n_ = n;
      Node *na = n;
      while (n_ != commonAncestor) {
        makePath(n_, w, n_->parent != commonAncestor);

        na = n_;
        n_ = n_->parent;
        assert(n_);
      }
      assert(na);

      // Connect paths
      verticalPath(ca, na, w);
    } else
      unaccounted += c.count();
  }

  if (unaccounted > 0) {
    QPoint labelPos = n->scenePos().toPoint();
    labelPos.setX(labelPos.x() + .025*R);
    QString label = QString::number(100. * unaccounted / totalWidth, 'f', 2)
                    + "% unaccounted";
    labels.append(QPair<QPointF, QString>(labelPos, label));
  }

  QGraphicsItem::show();
  update();
}

void Contributors::hide (void) {
  species = nullptr;
  paths.clear();
  QGraphicsItem::hide();
}

void Contributors::paint (QPainter *painter,
                          const QStyleOptionGraphicsItem*, QWidget*) {

  if (paths.empty())  return;

  QPen pen = tree->pathPen(details::PATH_CONTRIBUTOR);
  QColor c = pen.color();
  float R = tree->radius();

  painter->save();

    // Dim-out the rest of the graph
    painter->setBrush(QColor::fromRgbF(0,0,0,.5));
    painter->drawEllipse({0,0}, R, R);
    painter->setBrush(Qt::transparent);

    // Draw path with fading color
    for (const Path &p: paths) {
      pen.setColor(QColor::fromHsvF(c.hsvHueF(), p.width, .5 + .5 * p.width));
      painter->setPen(pen);
      painter->drawPath(p.path);
    }

    // Draw labels
    pen.setColor(Qt::black);
    QFont f = painter->font();
    f.setPixelSize(PTGraphBuilder::fontSize(R));
    painter->setFont(f);
    painter->setPen(pen);
    painter->setBackground(Qt::white);
    painter->setBackgroundMode(Qt::OpaqueMode);
    for (const auto &l: labels) {
      painter->drawText(l.first, l.second);
    }

  painter->restore();
}

// ============================================================================
// == Tree dimmer
// ============================================================================

Dimmer::Dimmer (VTree tree) : tree(tree) {
  setZValue(DIMMER_LEVEL);
}

QRectF Dimmer::boundingRect(void) const {
  return tree->boundingRect();
}

void Dimmer::paint (QPainter *painter,
                          const QStyleOptionGraphicsItem*, QWidget*) {
  if (dimPath.elementCount() > 0) {
    painter->setBrush(QColor::fromRgbF(0,0,0,.5));
    painter->drawPath(dimPath);
  }
}


// ============================================================================
// == Graph borders and legends
// ============================================================================

Border::Border (VTree tree, double height)
  : tree(tree), radius(height),
    pen(Qt::gray, 1, Qt::DashLine),
    font("Courrier", 20),
    metrics(font) {

  empty = (height == 0);
  updateShape();
  setZValue(BOUNDS_LEVEL);
}

void Border::updateShape(void) {
  prepareGeometryChange();
  shape = QPainterPath();
  legend.clear();

  if (empty) {
    QString msg = "Waiting for input";
    QRect bounds = metrics.boundingRect(msg);
    shape.addText(-bounds.center(), font, msg);

  } else {
    QPointF p = radius * QPointF(cos(LEGEND_PHASE), sin(LEGEND_PHASE));

    shape.moveTo(0,0);
    shape.lineTo(p);

    for (uint i=1; i<=LEGEND_TICKS; i++) {
      double v = double(i) / LEGEND_TICKS;
      double h = radius * v;
      shape.addEllipse({0,0}, h, h);

      legend << QPair(i, v * p);
    }
  }
}

void format (float &n, float &d, float e) {
  float n_ = std::floor(n / e);
  d = 10 * (n - e * n_) / e;
  n = n_;
}

QString prettyNumber (float n) {
  if (n < 1e3)        return QString::number(n);
  else if (n > 1e12)  return QString::number(n, 'e', 3);

  float d = 0;
  QString u = "";
  if (n < 1e6)        format(n, d, 1e3), u = "K";
  else if (n < 1e9)   format(n, d, 1e6), u = "M";
  else if (n < 1e12)  format(n, d, 1e9), u = "G";

  QString res = QString::number(n) + u;
  if (d > 0)  res += " " + QString::number(d);
  return res;
}

/// \todo Stops working above 2'000'000
void Border::paint(QPainter *painter, const QStyleOptionGraphicsItem*, QWidget*) {
  painter->setPen(tree->pathPen(details::BORDER_AXIS));
  painter->setBrush(Qt::white);

  font.setPixelSize(PTGraphBuilder::fontSize(radius));

  metrics = QFontMetrics(font);
  painter->setFont(font);

  if (radius > 0) {
    QList<QPair<QRectF, QString>> texts;
    QRegion clip (-radius, -radius, 2*radius, 2*radius);

    // Compute legend values and clip-out text areas
    for (auto &p: legend) {
      double v = p.first / double(LEGEND_TICKS);
      double h = radius * v;
      QString text = prettyNumber(h);
      QRectF textRect = QRectF(metrics.boundingRect(text)).translated(p.second);
      textRect.translate(-.5*textRect.width(), .5*textRect.height() - metrics.descent());
      clip = clip.subtracted(QRegion(textRect.toAlignedRect()));
      texts << QPair(textRect, text);
    }

    // Draw disks in reverse order (so that they overlap correctly)
    painter->save();
    painter->setPen(Qt::transparent);
    for (auto it=legend.rbegin(); it!=legend.rend(); ++it) {
      const auto &p  = *it;
      double v = p.first / double(LEGEND_TICKS);
      double h = radius * v;

      QRectF rect (-h, -h, 2*h, 2*h);

      painter->setBrush(gui::mix(QColor(0, 255 * (1 - v), 255 * v), QColor(Qt::white), 16./255.));
      painter->drawEllipse(rect);
    }
    painter->restore();

    // Draw 'axis'
    painter->save();
    painter->setBrush(Qt::transparent);
    painter->setClipRegion(clip);
    painter->drawPath(shape);
    painter->restore();

    // Draw tic values
    for (auto &p: texts)
      painter->drawText(p.first, Qt::AlignCenter, p.second);

  } else {
    // Just draw the 'pending' message
    painter->setBrush(Qt::gray);
    painter->drawPath(shape);
  }
}


// ============================================================================
// == Graph builder
// ============================================================================

PTGraphBuilder::PenSet PTGraphBuilder::buildPenSet (void) {
  PenSet map;

  QPen base (PATH_DEFAULT_COLOR, PATH_WIDTH,
             Qt::SolidLine, Qt::RoundCap, Qt::RoundJoin);

  QPen survivor = base;
  survivor.setColor(PATH_SURVIVOR_COLOR);

  QPen contributor = base;
  contributor.setColor(PATH_CONTRIBUTOR_COLOR);

  QPen border (Qt::gray, 1, Qt::DashLine);

  map[details::PATH_BASE] = base;
  map[details::PATH_SURVIVOR] = survivor;
  map[details::PATH_CONTRIBUTOR] = contributor;
  map[details::BORDER_AXIS] = border;
  return map;
}

void PTGraphBuilder::updatePenSet(float radius, PenSet &pens) {
  pens[details::PATH_BASE].setWidthF(pathWidth(PATH_WIDTH, radius));
  pens[details::PATH_SURVIVOR].setWidthF(pens[details::PATH_BASE].widthF());
  pens[details::PATH_CONTRIBUTOR].setWidthF(2.*pens[details::PATH_BASE].widthF());
  pens[details::BORDER_AXIS].setWidthF(pathWidth(AXIS_WIDTH, radius));
}

float PTGraphBuilder::nodeWidth(float radius) {
  return radius / (NODE_SIZE * 20);
}

float PTGraphBuilder::pathWidth(float baseWidth, float radius) {
  return baseWidth * radius / 400.;
}

float PTGraphBuilder::fontSize(float radius) {
  return std::max(1.f, radius / 50);
}

void PTGraphBuilder::updateLayout (Node *localRoot, PolarCoordinates &pc) {
  if (localRoot->subtreeVisible()) {
    localRoot->invalidate(pc(localRoot->data.firstAppearance));

    for (gui::Node *n: localRoot->subnodes)
      updateLayout(n, pc);
  }
}

void PTGraphBuilder::updateLayout (GUIItems &items) {
  uint visible = 0;
  for (const gui::Node *n: items.nodes)  visible += n->subtreeVisible();

  if (visible > 0) {
    double width = PolarCoordinates::xCoord(visible);
    PolarCoordinates pc (width);
    updateLayout(items.root, pc);
  }
}

} // end of namespace gui
