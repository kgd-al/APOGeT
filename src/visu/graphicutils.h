#ifndef GRAPHICUTILS_H
#define GRAPHICUTILS_H

#include <QColor>
#include <QRectF>

/*!
 * \file graphicutils.h
 *
 * Contains utilities for dealing with Qt
 */

namespace gui {

/// \returns a mixture of \p lhs and \p rhs with ratio \p r
QColor mix (const QColor &lhs, const QColor &rhs, double r);

/// \returns a rectangle with the same ratio as \p inner expanded and centered
/// into \p outter
QRectF centeredInto (const QRectF &outter, const QRectF &inner);

}

#endif // GRAPHICUTILS_H
