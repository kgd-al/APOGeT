#ifndef GRAPHICUTILS_H
#define GRAPHICUTILS_H

#include <QColor>

/*!
 * \file graphicutils.h
 *
 * Contains utilities for dealing with Qt
 */

namespace gui {

/// \returns a mixture of \p lhs and \p rhs with ratio \p r
QColor mix (const QColor &lhs, const QColor &rhs, double r);

}

#endif // GRAPHICUTILS_H
