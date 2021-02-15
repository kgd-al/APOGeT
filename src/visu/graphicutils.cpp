#include "graphicutils.h"

namespace gui {

QColor mix (const QColor &lhs, const QColor &rhs, double r) {
  return QColor(
    r * lhs.red() + (1-r) * rhs.red(),
    r * lhs.green() + (1-r) * rhs.green(),
    r * lhs.blue() + (1-r) * rhs.blue()
  );
}

QRectF centeredInto (const QRectF &outter, const QRectF &inner) {
  QSizeF requestedSize = outter.size();
  QSizeF actualSize = inner.size();
  float r = std::min(
    requestedSize.width() / actualSize.width(),
    requestedSize.height() / actualSize.height()
  );
  actualSize *= r;

  return QRect (
    (requestedSize.width() - actualSize.width()) / 2,
    (requestedSize.height() - actualSize.height()) / 2,
    actualSize.width(), actualSize.height()
  );
}

}
