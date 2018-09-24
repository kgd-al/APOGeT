#include <QTextStream>

#include <QPainter>
#include <qmath.h>

#include "ptgraphbuilder.h"
#include "graphicutils.h"

#include <QDebug>

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
static constexpr float NODE_SPACING = NODE_SIZE / 2;

static constexpr float END_POINT_SIZE = NODE_RADIUS / 4;

// == Z-values ======================================================
static constexpr int NODE_SURVIVOR_LEVEL = 2;
static constexpr int NODE_EXCTINCT_LEVEL = 1;

static constexpr int PATHS_LEVEL = -1;
static constexpr int TIMELINES_LEVEL = -2;
static constexpr int BOUNDS_LEVEL = -10;

// == Paint style ===================================================

static constexpr float PATH_WIDTH = 2.5;

static constexpr Qt::GlobalColor PATH_DEFAULT_COLOR = Qt::darkGray;
static constexpr Qt::GlobalColor PATH_HIGHLIGHT_COLOR = Qt::red;

static const QPen BASE_PEN (PATH_DEFAULT_COLOR, PATH_WIDTH,
                            Qt::SolidLine, Qt::RoundCap, Qt::RoundJoin);

static const QPen HIGH_PEN = [] {
  QPen p = BASE_PEN;
  p.setColor(PATH_HIGHLIGHT_COLOR);
  return p;
}();


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

  /// \returns the coordinate of node number \p i
  static float xCoord (uint i) {
    return i * NODE_SIZE + (i > 0 ? i-1 : 0) * NODE_SPACING;
  }

  /// Creates a polar coordinates object with the specified unscaled \p width
  PolarCoordinates (double width)
    : width(2 * M_PI * width / (2 * M_PI - LEGEND_SPACE)),
      nextX(0) {}

  /// \returns the position of the next point
  QPointF operator() (uint time) {
    double a = phase;
    if (width > 0)  a += 2. * M_PI * xCoord(nextX++) / width;
  //  qDebug() << "theta(" << nextX-1 << "): " << a * 180. / M_PI;
    return time * QPointF(cos(a), sin(a));
  }
};


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
  setScale(fullness());
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
}

void Node::setOnSurvivorPath (bool osp) {
  _onSurvivorPath = osp;
  setZValue(_onSurvivorPath ? NODE_SURVIVOR_LEVEL : NODE_EXCTINCT_LEVEL);
}

void Node::paint (QPainter *painter, const QStyleOptionGraphicsItem*, QWidget*) {
  if (visibilities.testFlag(SHOW_NAME)) {
    QRectF r (boundingRect().center() - QPoint(NODE_RADIUS, NODE_RADIUS),
              2*QSizeF(NODE_RADIUS, NODE_RADIUS));

    painter->setBrush(Qt::white);
    painter->setPen(_onSurvivorPath ? Qt::red : Qt::black);
    painter->drawEllipse(boundingRect().center(), NODE_RADIUS, NODE_RADIUS);
    painter->setPen(Qt::black);
    painter->drawText(r, Qt::AlignCenter, sid);
  }
}


// ============================================================================
// == Path between a parent and child node
// ============================================================================

Path::Path(Node *start, Node *end) : _start(start), _end(end) {}

void Path::invalidatePath(void) {
  prepareGeometryChange();

  _shape = QPainterPath();
  _shape.setFillRule(Qt::WindingFill);

  double a0 = PolarCoordinates::primaryAngle(_start->pos());

  const QPointF p1 = _end->scenePos();
  double a1 = PolarCoordinates::primaryAngle(p1);
  double r1 = PolarCoordinates::length(p1);

  _shape.moveTo(p1);
  _shape.arcTo(QRect(-r1, -r1, 2*r1, 2*r1), -qRadiansToDegrees(a1), qRadiansToDegrees(a1 - a0));
  _shape.addEllipse(_shape.pointAtPercent(1), END_POINT_SIZE, END_POINT_SIZE);
}

QRectF Path::boundingRect() const {
  qreal extra = (PATH_WIDTH + 20) / 2.0;

  return shape().boundingRect().normalized()
      .adjusted(-extra, -extra, extra, extra);
}

void Path::paint(QPainter *painter, const QStyleOptionGraphicsItem *, QWidget *) {
  painter->setPen(_end->onSurvivorPath() ? HIGH_PEN : BASE_PEN);
  painter->drawPath(_shape);
}


// ============================================================================
// == Timeline for a node
// ============================================================================

Timeline::Timeline(Node *node) : _node(node) {
  setZValue(TIMELINES_LEVEL);
}

void Timeline::invalidatePath(void) {
  prepareGeometryChange();

  _points[0] = _node->scenePos();
  double a = PolarCoordinates::primaryAngle(_points[0]);

  _points[2] = _node->data.lastAppearance * QPointF(cos(a), sin(a));

  if (_node->alive()) // Survivor timeline goes all the way
    _points[1] = _points[2];

  else if (!_node->onSurvivorPath())  // No survivor part for this path
    _points[1] = _points[0];

  else {
    // Compute location of last child survivor
    double l = _node->data.firstAppearance;
    for (const Node *gn: _node->subnodes) {
      if (!gn->onSurvivorPath())  continue;

      double l_ = gn->data.firstAppearance;
      if (l < l_) l = l_;
    }

    _points[1] = l * QPointF(cos(a), sin(a));
  }
}

QPainterPath Timeline::shape (void) const {
  QPainterPath path;
  path.moveTo(_points[0]);
  path.lineTo(_points[1]);
  path.lineTo(_points[2]);
  path.addEllipse(_points[2], END_POINT_SIZE, END_POINT_SIZE);
  return QPainterPathStroker (BASE_PEN).createStroke(path);
}

void Timeline::paint(QPainter *painter, const QStyleOptionGraphicsItem *, QWidget *) {
  if (_points[0] != _points[1]) {
    painter->setPen(HIGH_PEN);
    painter->drawLine(_points[0], _points[1]);
  }

  if (_points[1] != _points[2]) {
    painter->setPen(BASE_PEN);
    painter->drawLine(_points[1], _points[2]);
  }

  painter->setBrush(painter->pen().color());
  painter->drawEllipse(_points[2], END_POINT_SIZE, END_POINT_SIZE);
}


// ============================================================================
// == Graph borders and legends
// ============================================================================

Border::Border (double height)
  : height(height),
    pen(Qt::gray, 1, Qt::DashLine),
    font("Courrier", 20),
    metrics(font) {

  empty = (height == 0);
  updateShape();
  setZValue(-10);
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
    QPointF p = height * QPointF(cos(LEGEND_PHASE), sin(LEGEND_PHASE));

    shape.moveTo(0,0);
    shape.lineTo(p);

    for (uint i=1; i<=LEGEND_TICKS; i++) {
      double v = double(i) / LEGEND_TICKS;
      double h = height * v;
      shape.addEllipse({0,0}, h, h);

      legend << QPair(i, v * p);
    }
  }
}

void Border::paint(QPainter *painter, const QStyleOptionGraphicsItem*, QWidget*) {
  painter->setPen(pen);
  painter->setBrush(Qt::white);

  float factor = height / 50;
  float f = std::max(1.f, factor);
  font.setPointSizeF(f);

  metrics = QFontMetrics(font);
  painter->setFont(font);

  if (!config::PTree::winningPathOnly()) {
    if (height > 0) {
      QList<QPair<QRectF, QString>> texts;
      QRegion clip (-height, -height, 2*height, 2*height);

      // Compute legend values and clip-out text areas
      for (auto &p: legend) {
        double v = p.first / double(LEGEND_TICKS);
        double h = height * v;
        QString text = QString::number(h);
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
        double h = height * v;

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

  qDebug() << "Updating layout for" << visible << "items";

  if (visible > 0) {
    double width = PolarCoordinates::xCoord(visible-1);
    PolarCoordinates pc (width);
    updateLayout(items.root, pc);
  }
}

} // end of namespace gui
