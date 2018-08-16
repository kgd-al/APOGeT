#include "graphicutils.h"

namespace gui {

QColor mix (const QColor &lhs, const QColor &rhs, double r) {
  return QColor(
    r * lhs.red() + (1-r) * rhs.red(),
    r * lhs.green() + (1-r) * rhs.green(),
    r * lhs.blue() + (1-r) * rhs.blue()
  );
}

}
