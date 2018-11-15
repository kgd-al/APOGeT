#include <QTextStream>

#include <QPainter>
#include <QToolTip>
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
static constexpr int NODE_EXCTINCT_LEVEL = 10;

static constexpr int CONTRIBUTORS_LEVEL = 0;

static constexpr int PATHS_LEVEL = -10;
static constexpr int TIMELINES_LEVEL = -11;
static constexpr int BOUNDS_LEVEL = -20;

// == Paint style ===================================================

static constexpr float AXIS_WIDTH = 1;
static constexpr float PATH_WIDTH = 1.5;

static constexpr Qt::GlobalColor PATH_DEFAULT_COLOR = Qt::darkGray;
static constexpr Qt::GlobalColor PATH_SURVIVOR_COLOR = Qt::red;
static constexpr Qt::GlobalColor PATH_CONTRIBUTOR_COLOR = Qt::green;


// ============================================================================
// == Coordinate computation
// ============================================================================

/// Generates and manages polar coordinates for the nodes/paths
struct PolarCoordinates {

  /// The angular phase used in coordinate computation.
  /// Inverted to cope with Qt coordinate system.
  static constexpr float phase = LEGEND_PHASE + LEGEND_SPACE/2.;

  const double width; ///< Width of the graph (legend included)

  uint nextX; ///< X coordinate of the next node

  /// \returns The angle for \p in the range [phase,2&pi;+phase]
  static double primaryAngle (const QPointF &p) {
    if (p.isNull()) return phase;
    double a = atan2(p.y(), p.x());
    while (a < phase) a += 2 * M_PI;
    while (a > 2 * M_PI + phase) a -= 2 * M_PI;
    return a;
  }

  /// \returns The euclidian distance of \p p with the origin
  static double length (const QPointF &p) {
    return sqrt(p.x()*p.x() + p.y()*p.y());
  }

  static QPointF toCartesian (double a, double r) {
    return r * QPointF(cos(a), sin(a));
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

QPainterPath makeArc (const QPointF &p0, const QPointF &p1) {
  double a0 = PolarCoordinates::primaryAngle(p0);
  double a1 = PolarCoordinates::primaryAngle(p1);
  double r1 = PolarCoordinates::length(p1);

  QPainterPath path;
  path.moveTo(p1);
  path.arcTo(QRect(-r1, -r1, 2*r1, 2*r1), -qRadiansToDegrees(a1),
             qRadiansToDegrees(a1 - a0));

  return path;
}

QPointF timelineAnchor (const Node *n) {
  return PolarCoordinates::toCartesian(
    PolarCoordinates::primaryAngle(n->parent->scenePos()),
    PolarCoordinates::length(n->scenePos()));
}


// ============================================================================
// == Graph node
// ============================================================================

QRectF Node::boundingRect(void) const {
  return QRectF(-NODE_SIZE, -NODE_SIZE, 2*NODE_SIZE, 2*NODE_SIZE);
}

void Node::invalidate(const QPointF &newPos) {
  setPos(newPos);
  if (path) path->invalidatePath();
  timeline->invalidatePath();
}

void Node::updateTooltip (void) {
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
  setToolTip(tooltip);
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

  timeline->invalidatePath();
  updateTooltip();
  autoscale();
}

void Node::setOnSurvivorPath (bool osp) {
  bool s = _onSurvivorPath = osp;
  setZValue(s ? NODE_SURVIVOR_LEVEL : NODE_EXCTINCT_LEVEL);
  if (path) path->setZValue(PATHS_LEVEL + s);
}

void Node::hoverEnterEvent(QGraphicsSceneHoverEvent*) {
  treeBase->hoverEvent(id, true);
}

/// Triggers a callback when this species node is no longer hovered
void Node::hoverLeaveEvent(QGraphicsSceneHoverEvent*) {
  treeBase->hoverEvent(id, false);
}

void Node::paint (QPainter *painter, const QStyleOptionGraphicsItem*, QWidget*) {
  if (visibilities.testFlag(SHOW_NAME)) {
    QRectF r (boundingRect().center() - QPoint(NODE_RADIUS, NODE_RADIUS),
              2*QSizeF(NODE_RADIUS, NODE_RADIUS));

    painter->save();
      painter->setBrush(Qt::white);
      painter->setPen(_onSurvivorPath ? Qt::red : Qt::black);
      painter->drawEllipse(boundingRect().center(), NODE_RADIUS, NODE_RADIUS);
      painter->setPen(Qt::black);
      QFont f = painter->font();
      f.setPixelSize(12);
      painter->setFont(f);
      painter->drawText(r, Qt::AlignCenter, sid);
    painter->restore();
  }
}


// ============================================================================
// == Path between a parent and child node
// ============================================================================

Path::Path(Node *start, Node *end) : start(start), end(end) {
  setZValue(PATHS_LEVEL);
  assert(start && end);
  assert(start->treeBase == end->treeBase);
}

void Path::invalidatePath(void) {
  prepareGeometryChange();

  _shape = QPainterPath();
  _shape.setFillRule(Qt::WindingFill);

  _shape.addPath(makeArc(start->scenePos(), end->scenePos()));
  _shape.addEllipse(_shape.pointAtPercent(1), END_POINT_SIZE, END_POINT_SIZE);
}

QRectF Path::boundingRect() const {
  qreal extra = (PATH_WIDTH + 20) / 2.0;

  return shape().boundingRect().normalized()
      .adjusted(-extra, -extra, extra, extra);
}

void Path::paint(QPainter *painter, const QStyleOptionGraphicsItem*, QWidget*) {
  painter->setPen(start->treeBase->pathPen(
                    end->onSurvivorPath() ? details::PATH_SURVIVOR
                                          : details::PATH_BASE));
  painter->drawPath(_shape);
}


// ============================================================================
// == Timeline for a node
// ============================================================================

Timeline::Timeline(Node *node) : node(node) {
  setZValue(TIMELINES_LEVEL);
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
}

QPainterPath Timeline::shape (void) const {
//  QPainterPath path;
//  path.moveTo(points[0]);
//  path.lineTo(points[1]);
//  path.lineTo(points[2]);
//  path.addEllipse(points[2], END_POINT_SIZE, END_POINT_SIZE);
//  return QPainterPathStroker (BASE_PEN).createStroke(path);
  return QPainterPath();
}

void Timeline::paint(QPainter *painter, const QStyleOptionGraphicsItem *, QWidget *) {
  if (points[0] != points[1]) {
    painter->setPen(node->treeBase->pathPen(details::PATH_SURVIVOR));
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

  // Total contribution count (excluding itself)
  float totalWidth = 0;
  for (auto &c: contribs)
    if (c.speciesID() != sid)  totalWidth += c.count();

  paths.clear();
  labels.clear();
  for (auto &c: contribs) {
    if (c.speciesID() == sid) continue;

    float w = c.count() / totalWidth;
    Node *nc = items.nodes.value(c.speciesID());

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

  if (!config::PTree::winningPathOnly()) {
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
