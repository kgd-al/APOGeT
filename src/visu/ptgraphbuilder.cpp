#include <QPainter>
#include <qmath.h>

#include "ptgraphbuilder.h"
#include "graphicutils.h"

#include <QDebug>

namespace gui {

// ============================================================================
// == Constants
// ============================================================================

///< Note that the coordinate system is inverted (y points downward)
static constexpr float LEGEND_PHASE = -M_PI/2;
static constexpr float LEGEND_SPACE = M_PI / 2;
static constexpr uint LEGEND_TICKS = 4;

static constexpr float NODE_RADIUS = 10;
static constexpr float NODE_MARGIN = 2;
static constexpr float NODE_SIZE = 2 * (NODE_RADIUS + NODE_MARGIN);
static constexpr float NODE_SPACING = NODE_SIZE / 2;

static constexpr float END_POINT_SIZE = NODE_RADIUS / 4;


// ============================================================================
// == Coordinate computation
// ============================================================================

PTreeBuildingCache::PolarCoordinates::PolarCoordinates (double width)
  : originalWidth(width),
    widthWithLegend(2 * M_PI * width / (2 * M_PI - LEGEND_SPACE)),
    phase(LEGEND_PHASE + LEGEND_SPACE/2.),
    nextX(0) {}

float PTreeBuildingCache::PolarCoordinates::xCoord (uint i) {
  return i * NODE_SIZE + (i > 0 ? i-1 : 0) * NODE_SPACING;
}

double PTreeBuildingCache::PolarCoordinates::toPrimaryAngle (const QPointF &p) const {
  if (p.isNull()) return phase;
  double a = atan2(p.y(), p.x());
  while (a < phase) a += 2 * M_PI;
  while (a > 2 * M_PI + phase) a -= 2 * M_PI;
  return a;
}

QPointF PTreeBuildingCache::PolarCoordinates::operator() (uint time) {
  double a = phase + 2. * M_PI * xCoord(nextX++) / widthWithLegend;
  return time * QPointF(cos(a), sin(a));
}


// ============================================================================
// == Graph node
// ============================================================================

QRectF Node::boundingRect(void) const {
  return QRectF(-NODE_SIZE, -NODE_SIZE, 2*NODE_SIZE, 2*NODE_SIZE);
}

void Node::updateTooltip (void) {
  QString tooltip;
  QTextStream qss (&tooltip);
  qss << "Node " << id << "\n"
      << "Enveloppe: " << fullness << "\n"
      << "Appeared at " << data.firstAppearance << "\n"
      << "Disappeared at " << data.lastAppearance << "\n"
      << data.count << " individuals\n"
      << children << " subspecies";
  setToolTip(tooltip);
}

void Node::paint (QPainter *painter, const QStyleOptionGraphicsItem*, QWidget*) {
  QRectF r (boundingRect().center() - QPoint(NODE_RADIUS, NODE_RADIUS),
            2*QSizeF(NODE_RADIUS, NODE_RADIUS));

  painter->setBrush(Qt::white);
  painter->setPen(_stillAlive ? Qt::red : Qt::black);
  painter->drawEllipse(boundingRect().center(), NODE_RADIUS, NODE_RADIUS);
  painter->setPen(Qt::black);
  painter->drawText(r, Qt::AlignCenter, QString::number(id));
}


// ============================================================================
// == Graph path
// ============================================================================

Path::Path(Node *start, Node *end, const Cache &cache) : _start(start), _end(end) {
  _pen = QPen(Qt::darkGray, 2.5, Qt::SolidLine, Qt::RoundCap, Qt::RoundJoin);

  _shape = QPainterPath();
  _shape.setFillRule(Qt::WindingFill);

  _survivor = QPainterPath();
  _survivor.setFillRule(Qt::WindingFill);

  double S = _end->isStillAlive();
  setZValue(S ? -1 : -2);

  const QPointF p = _end->scenePos();
  double a1 = cache.polarCoordinate.toPrimaryAngle(p);
  double r = length(p);

  const QPointF p2 = p + (_end->survival()) * QPointF(cos(a1), sin(a1));

  if (_start && (S || !PTreeConfig::winningPathOnly())) {
    double a0 = cache.polarCoordinate.toPrimaryAngle(_start->pos());
    QPainterPath &path = S ? _survivor : _shape;
    path.moveTo(p);
    path.arcTo(QRect(-r, -r, 2*r, 2*r), -qRadiansToDegrees(a1), qRadiansToDegrees(a1 - a0));
    path.addEllipse(path.pointAtPercent(1), END_POINT_SIZE, END_POINT_SIZE);
  }

  QPointF p1 = p;
  if (S) {
//    double l = length(p);
//    for (const auto &pn: _end->node.children) {
//      if (pn2gn.find(pn->id) == pn2gn.end())  continue;

//      Node *gn = pn2gn[pn->id];
//      if (!gn->survivor)  continue;

//      QPointF p_ = gn->pos();
//      double l_ = length(p_);
//      if (l < l_) {
//        p1 = p_;
//        l = l_;
//      }
//    }

//    if (p1 == p) p1 = p2;
//    else
//      p1 = length(p1) * QPointF(cos(a1), sin(a1));

//    assert(!_end->survivor || _end->node.children.empty() || p != p1);
//    _survivor.moveTo(p);
//    _survivor.lineTo(p1);
  }

  bool atEnd = length(p1) < length(p2);
  if (!atEnd || !PTreeConfig::winningPathOnly()){
    QPainterPath &path = atEnd ? _shape : _survivor;
    path.moveTo(p1);
    path.lineTo(p2);
    path.addEllipse(p2, END_POINT_SIZE, END_POINT_SIZE);
  }
}

QRectF Path::boundingRect() const {
  qreal extra = (_pen.width() + 20) / 2.0;

  return shape().boundingRect().normalized()
      .adjusted(-extra, -extra, extra, extra);
}

void Path::paint(QPainter *painter, const QStyleOptionGraphicsItem *, QWidget *) {
  painter->setPen(_pen);
  painter->drawPath(_shape);

  QPen pen (_pen);
  pen.setColor(Qt::red);
  painter->setPen(pen);
  painter->drawPath(_survivor);
}


// ============================================================================
// == Graph borders and legends
// ============================================================================

Border::Border (double height)
  : height(height),
    pen(Qt::gray, 1, Qt::DashLine),
    font("Courrier", 20),
    metrics(font) {

  double a = -.5 * M_PI;
  QPointF p = height * QPointF(cos(a), sin(a));

  shape = QPainterPath();
  shape.moveTo(0,0);
  shape.lineTo(p);

  for (uint i=1; i<=LEGEND_TICKS; i++) {
    double v = double(i) / LEGEND_TICKS;
    double h = height * v;
    shape.addEllipse({0,0}, h, h);

    legend << QPair(i, v * p);
  }

  setZValue(-10);
}

void Border::paint(QPainter *painter, const QStyleOptionGraphicsItem*, QWidget*) {
  painter->setPen(pen);
  painter->setBrush(Qt::white);

  float factor = height / 50;
  float f = std::max(1.f, factor);
  font.setPointSizeF(f);

  metrics = QFontMetrics(font);
  painter->setFont(font);

  if (!PTreeConfig::winningPathOnly()) {
    QRegion clip (-height, -height, 2*height, 2*height);
    QList<QPair<QRectF, QString>> texts;
    for (auto &p: legend) {
      double v = p.first / double(LEGEND_TICKS);
      double h = height * v;
      QString text = QString::number(h);
      QRectF textRect = QRectF(metrics.boundingRect(text)).translated(p.second);
      textRect.translate(-.5*textRect.width(), .5*textRect.height() - metrics.descent());
      clip = clip.subtracted(QRegion(textRect.toAlignedRect()));
      texts << QPair(textRect, text);
    }

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

    painter->save();
    painter->setBrush(Qt::transparent);
    painter->setClipRegion(clip);
    painter->drawPath(shape);
    painter->restore();

    for (auto &p: texts)
      painter->drawText(p.first, Qt::AlignCenter, p.second);
  }
}

} // end of namespace gui
