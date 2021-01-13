#ifndef KGD_APOGET_LAYER_HPP
#define KGD_APOGET_LAYER_HPP

/*!
 * \file layer.hpp
 *
 * Contains the definition a graphics scene layer item that can draw anything
 * through a QPainter object
 */

#include <cassert>
#include <functional>

#include <QGraphicsItem>

namespace gui {

/// Graphics scene layer item used to draw a bunch of stuff on the same plane.
/// Only one such drawing is allowed per graphic item
struct Layer : public QGraphicsItem {

  /// Helper alias to the drawer function type
  using Drawer = std::function<void(QPainter*)>;

  /// Describes a drawing
  struct Drawing {
    Drawer drawer;  ///< Drawer function object
    bool doDraw;    ///< Whether or not drawing is enabled
  };

  /// Storage space for the drawing objects
  QMap<QGraphicsItem*, Drawing> drawings;

  /// Whether or not to activate uniform alpha blending, i.e. setting, for each
  /// item, an alpha value of 1/n with n, the number of items to draw
  bool uniformAlphaBlending = false;

  /// Constructor
  Layer (QGraphicsItem *parent, int zvalue) : QGraphicsItem(parent) {
    setZValue(zvalue);
  }

  /// (De)Activates uniform alpha blending
  void setUniformAlphaBlending (bool uniform) {
    uniformAlphaBlending = uniform;
  }

  /// Register a new drawing object
  void addDrawing (QGraphicsItem *item, const Drawer &drawer) {
    drawings[item] = Drawing{drawer, true};
    update();
  }

  /// Change the visibility of an already registered drawing
  void showDrawing (QGraphicsItem *item, bool show) {
    assert(drawings.find(item) != drawings.end());
    drawings[item].doDraw = show;
    update();
  }

  /// Unregister a drawing object
  void removeDrawing (QGraphicsItem *item) {
    drawings.remove(item);
    update();
  }

  /// Overrides QGraphicsItem::boundingRect
  /// \returns a rectangle occupying a much space as its parent's
  QRectF boundingRect(void) const {
    return parentItem()->boundingRect();
  }

  /// Delegate drawing to all currently active drawers
  void paint (QPainter *painter, const QStyleOptionGraphicsItem*, QWidget*) {
    if (drawings.empty()) return;

    if (uniformAlphaBlending)
      setOpacity(1. / drawings.size());

    for (Drawing &d: drawings) {
      painter->save();
      if (d.doDraw) d.drawer(painter);
      painter->restore();
    }
  }
};

} // end of namespace gui

Q_DECLARE_METATYPE(gui::Layer*)

#endif // KGD_APOGET_LAYER_HPP
