#ifndef _PTREE_INTROSPECTER_HPP_
#define _PTREE_INTROSPECTER_HPP_

#include <QTextStream>
#include <qmath.h>

#include <QGraphicsItem>

#include "../core/phylogenictree.hpp"

#include "graphicutils.h"

namespace phylogeny {

template <typename GENOME>
struct PTreeIntrospecter {
  using PT = PhylogenicTree<GENOME>;
  using PN = typename PT::Node;

  struct Node : public QGraphicsItem {
    static constexpr int S = 10;
    static constexpr int M = 2;

    QRectF subtree;
    const PN &node;
    bool survivor;

    Node (const QPointF pos, const PN &n) : node(n) {
      setPos(pos);
      subtree = boundingRect().translated(pos);

      QString tooltip;
      QTextStream qss (&tooltip);
      qss << "Node " << n.id << "\n"
          << "Enveloppe: " << n.enveloppe.size() << " / " << PTreeConfig::enveloppeSize() << "\n"
          << "Appeared at " << n.data.firstAppearance << "\n"
          << "Disappeared at " << n.data.lastAppearance << "\n"
//          << "x bounds: [" << n.data.xmin << "-" << n.data.xmax << "]\n"
          << n.data.count << " individuals\n"
          << n.children.size() << " subspecies";
      setToolTip(tooltip);
    }

    static double maxWidth (void) {
      return 2 * S + 2 * M;
    }

    QRectF boundingRect(void) const {
      return QRectF(-S, -S, 2*S, 2*S).adjusted(-M, -M, M, M);
    }

    void paint (QPainter *painter, const QStyleOptionGraphicsItem*, QWidget*) {
      QRectF r (boundingRect().center() - QPoint(S,S), 2*QSizeF(S,S));
      painter->setBrush(Qt::white);
      painter->setPen(survivor ? Qt::red : Qt::black);
      painter->drawEllipse(boundingRect().center(), S, S);
      painter->setPen(Qt::black);
      painter->drawText(r, Qt::AlignCenter, QString::number(node.id));

      if (false) {
        painter->setPen(Qt::blue);
        painter->drawRect(boundingRect());
      }

      if (false) {
        painter->setPen(Qt::red);
        painter->save();
        painter->translate(-pos());
        painter->drawRect(subtree);
        painter->restore();
  //      qDebug() << id << subtree;
      }
    }
  };

  struct Path : public QGraphicsItem {
    Node *_start, *_end;
    QPainterPath _shape, _survivor;
    QPen _pen;

    Path(Node *start, Node *end) : _start(start), _end(end) {
      _pen = QPen(Qt::darkGray, 2.5, Qt::SolidLine, Qt::RoundCap, Qt::RoundJoin);
      setZValue(-1);

      _shape = linearShape();
    }

    void toCircular (double da, const QMap<uint, Node*> &pn2gn) {
      _shape = QPainterPath();
      _shape.setFillRule(Qt::WindingFill);

      _survivor = QPainterPath();
      _survivor.setFillRule(Qt::WindingFill);

      double S = _end->survivor;

      const QPointF p = _end->scenePos();
      double a1 = toPrimaryAngle(p, da);
      double r = length(p);

      const QPointF p2 = p + (_end->node.data.lastAppearance - _end->node.data.firstAppearance) * QPointF(cos(a1), sin(a1));

      if (_start && (S || !PTreeConfig::winningPathOnly())) {
        double a0 = toPrimaryAngle(_start->pos(), da);
        QPainterPath &path = S ? _survivor : _shape;
        path.moveTo(p);
        path.arcTo(QRect(-r, -r, 2*r, 2*r), -qRadiansToDegrees(a1), qRadiansToDegrees(a1 - a0));
        path.addEllipse(path.pointAtPercent(1), Node::S/4, Node::S/4);
      }

      QPointF p1 = p;
      if (S) {
        double l = length(p);
        for (const auto &pn: _end->node.children) {
          if (pn2gn.find(pn->id) == pn2gn.end())  continue;

          Node *gn = pn2gn[pn->id];
          if (!gn->survivor)  continue;

          QPointF p_ = gn->pos();
          double l_ = length(p_);
          if (l < l_) {
            p1 = p_;
            l = l_;
          }
        }

        if (p1 == p) p1 = p2;
        else
          p1 = length(p1) * QPointF(cos(a1), sin(a1));

//        assert(!_end->survivor || _end->node.children.empty() || p != p1);
        _survivor.moveTo(p);
        _survivor.lineTo(p1);
      }

      bool atEnd = length(p1) < length(p2);
      if (!atEnd || !PTreeConfig::winningPathOnly()){
        QPainterPath &path = atEnd ? _shape : _survivor;
        path.moveTo(p1);
        path.lineTo(p2);
        path.addEllipse(p2, Node::S/4, Node::S/4);
      }
    }

    QRectF boundingRect() const {
      qreal extra = (_pen.width() + 20) / 2.0;

      return shape().boundingRect().normalized()
          .adjusted(-extra, -extra, extra, extra);
    }

    QPainterPath linearShape (void) {
      QPainterPath shape;
      const QPointF p = _end->scenePos();

      if (_start) {
        shape.addEllipse(QPointF(_start->scenePos().x(), p.y()), Node::S/4, Node::S/4);
        shape.moveTo(_start->scenePos().x(), p.y());
        shape.lineTo(p.x(), p.y());
      }

      shape.moveTo(p.x(), p.y());
      shape.lineTo(p.x(), _end->node.data.lastAppearance);
      shape.addEllipse(QPointF(p.x(), _end->node.data.lastAppearance), Node::S/4, Node::S/4);
      shape.setFillRule(Qt::WindingFill);

      return shape;
    }

    double length (const QPointF &p) {
      return sqrt(p.x()*p.x() + p.y()*p.y());
    }

    double toPrimaryAngle (const QPointF &p, double da) {
      if (p.isNull()) return da;
      double a = atan2(p.y(), p.x());
      while (a < da) a += 2 * M_PI;
      while (a > 2 * M_PI + da) a -= 2 * M_PI;
      return a;
    }

    QPainterPath shape (void) const override {
      return _shape.united(_survivor);
    }

    void paint(QPainter *painter, const QStyleOptionGraphicsItem *, QWidget *) {
      painter->setPen(_pen);
      painter->drawPath(_shape);

      QPen pen (_pen);
      pen.setColor(Qt::red);
      painter->setPen(pen);
      painter->drawPath(_survivor);
    }
  };

  struct Border : public QGraphicsItem {
    static constexpr int N = 4;
    double height;
    QPainterPath shape;
    QList<QPair<int, QPointF>> legend;

    QPen pen;
    QFont font;
    QFontMetrics metrics;

    Border (double height)
      : height(height),
        pen(Qt::gray, 1, Qt::DashLine),
        font("Courrier", 20),
        metrics(font) {

      double a = -.5 * M_PI;
      QPointF p = height * QPointF(cos(a), sin(a));

      shape = QPainterPath();
      shape.moveTo(0,0);
      shape.lineTo(p);

      for (int i=1; i<=N; i++) {
        double v = double(i) / N;
        double h = height * v;
        shape.addEllipse({0,0}, h, h);

        legend << QPair(i, v * p);
      }

      setZValue(-2);
    }

    QRectF boundingRect(void) const {
      return shape.boundingRect().adjusted(0, -metrics.ascent(), 0, 0);
    }

    void paint(QPainter *painter, const QStyleOptionGraphicsItem*, QWidget*) {
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
          double v = p.first / double(N);
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
          double v = p.first / double(N);
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
  };

  static float spacing (void) {
    return Node::maxWidth() / 2;
  }

  template <typename C>
  static void fillScene (const PT &pt, QGraphicsScene *scene, const C &config) {
    uint step = pt.step();
    double width = pt.width() * (Node::maxWidth() + spacing()) - spacing();

    if (pt._root)
      addSpecies(nullptr, *pt._root, scene, config, step);

    auto border = new Border(step);
    scene->addItem(border);

//    scene->setSceneRect(border->boundingRect());
  }

  template <typename C>
  static bool addSpecies(Node *parent, const PN &n, QGraphicsScene *scene, const C &config, uint T) {
    QString spacing;
    const PN *p = &n;
    while ((p = p->parent)) spacing += "  ";

    float survival = n.data.lastAppearance - n.data.firstAppearance;
    float size = n.enveloppe.size() / double(PTreeConfig::enveloppeSize());

    if (survival < config.minSurvival)  return false;
    if (size < config.minEnveloppe) return false;

    QPointF pos;
    if (parent) {
      pos += QPointF(
        parent->pos().x() + parent->subtree.width() + Node::M,
        float(n.data.firstAppearance)
      );
    }

    bool survivor = (n.data.lastAppearance >= T);
    Node *gn = new Node (pos, n);
    gn->setScale(size);
    gn->setVisible(config.showNames);

    scene->addItem(gn);

    for (auto it=n.children.rbegin(); it!=n.children.rend(); ++it)
      survivor |= addSpecies(gn, *(*it), scene, config, T);

    gn->survivor = survivor;
    scene->addItem(new Path(parent, gn));

    if (parent) parent->subtree = parent->subtree.united(gn->subtree);
    return survivor;
  }

  static void toCircular(QGraphicsScene *scene) {
    float originalWidth = scene->width();
    float width = originalWidth * 1.1;

    QMap<uint, Node*> pn2gn;
    for (auto *i: scene->items())
      if (Node *n = dynamic_cast<Node*>(i))
        pn2gn[n->node.id] = n;

    double da = -2 * M_PI * (.25 + (originalWidth + width) / (2 * width));

    for (auto *i: scene->items()) {
      if (Node *n = dynamic_cast<Node*>(i)) {
        QPointF pl = n->pos();
        double a = da + 2 * M_PI * pl.x() / width;
        double r = pl.y();
        QPointF pc (r * cos(a), r * sin(a));
        n->setPos(pc);

      } else if (Path *p = dynamic_cast<Path*>(i)) {
        p->toCircular(da, pn2gn);


      } else if (Border *b = dynamic_cast<Border*>(i)) {
        b->toCircular();
      }
    }
    scene->setSceneRect(scene->itemsBoundingRect());
  }
};

} // end of namespace phylogeny

#endif // _PTREE_INTROSPECTER_HPP_
