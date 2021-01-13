#ifndef KGD_SPECIESTRACKING_H
#define KGD_SPECIESTRACKING_H

#include <QStyledItemDelegate>
#include <QComboBox>
#include <QTableView>
#include <QStandardItemModel>
#include <QDialogButtonBox>
#include <QColorDialog>

#include "ptgraphbuilder.h"

namespace gui::species_tracking {
/// \cond internal
using Specs = decltype(ViewerConfig::colorSpecs);
using SID = decltype(ViewerConfig::ColorSpec::sid);
using Color = decltype(ViewerConfig::ColorSpec::color);

inline QString toString (SID sid) {
  return QString::number(uint(sid));
}

/// Manages visualisation and editing species identificators
class SIDDelegate : public QStyledItemDelegate {
  using Data = std::set<SID>; ///< Helper alias to a collection of species ID
  Data data;  ///< Species ID collection holder
  int longest;  ///< Length of the largest SID

public:
  /// Create a delegate from the provided data
  SIDDelegate(Data &&data, QWidget *parent = 0);

  /// \return the next species id to use
  auto nextSID (void) const {
    return *data.begin();
  }

  /// \return a large enough size to accomodate the largest provided SID
  QSize sizeHint(const QStyleOptionViewItem &/*option*/,
                 const QModelIndex &/*index*/) const override;

  /// Creates the editor (i.e. a combobox of the remaining SID)
  QWidget *createEditor(QWidget *parent, const QStyleOptionViewItem &/*option*/,
                        const QModelIndex &/*index*/) const override;

  /// Updates the editor with the current value
  void setEditorData(QWidget *editor, const QModelIndex &index) const override;

  /// Updates the underlying data model
  void setModelData(QWidget *editor, QAbstractItemModel *model,
                    const QModelIndex &index) const override;
};

/// Manages visualisation and editing of a species tracking color
class ColorDelegate : public QStyledItemDelegate {
  static const QList<QColor> defaultColors; ///< Gnuplot colors
public:
  /// Creates a color delegate
  ColorDelegate(QWidget *parent = 0) : QStyledItemDelegate(parent) {}

  /// \return the next color to use
  static auto nextColor (int i) {
    return defaultColors.at(i % defaultColors.size());
  }

  /// Sets up a color dialog to use the default colors
  static void setupColorDialog (QColorDialog &cd) {
    for (int i=0; i<defaultColors.size(); i++)
      cd.setCustomColor(i, defaultColors.at(i));
  }

  /// Paints as a colored button
  void paint(QPainter *painter, const QStyleOptionViewItem &option,
             const QModelIndex &index) const override;

  /// \returns a small size
  QSize sizeHint(const QStyleOptionViewItem &/*option*/,
                 const QModelIndex &/*index*/) const override {
    return QSize (20, 10);
  }

  /// Creates the editor (i.e. a QColorDialog)
  QWidget *createEditor(QWidget *parent, const QStyleOptionViewItem &/*option*/,
                        const QModelIndex &/*index*/) const override;

  /// Updates the eidtor with the current color
  void setEditorData(QWidget *editor, const QModelIndex &index) const override;

  /// Updates the underlying data model
  void setModelData(QWidget *editor, QAbstractItemModel *model,
                    const QModelIndex &index) const override;
};

/// A QTableView with limited vertical space use
class TightTableView : public QTableView {
  QSize sizeHint (void) const;  ///< The smallest height possible

public:
  /// \authors savolai & Ratah @ https://stackoverflow.com/a/50217711
  void verticalFit (void);
};

/// Manages species trackings definition (mapping SID with a color)
class Dialog : public QDialog {
  Q_OBJECT

  QStandardItemModel model; ///< The data model
  SIDDelegate sidDelegate;  ///< The species identificator manager
  ColorDelegate colorDelegate;  ///< The species color manager

  QDialogButtonBox *controls; ///< The control button

  /// Extract the collection of valid SID from a PTree
  auto validSIDs (const PhylogenyViewer_base *viewer) const;

  /// \return an initializer for the underlying data model
  auto columnFromSpec (const ViewerConfig::ColorSpec &cs);

  /// Reacts to data model changes
  void dataChanged (void);

public:
  /// Create a species tracking dialog from a pviewer and exiting color specs
  Dialog (PhylogenyViewer_base *viewer, const Specs &is);

  /// \return the update color specs
  Specs colorSelection (void) const;

signals:
  /// Warn connected objects the color specs have been modified
  void applied (void);
};

/// \endcond internal
} // end of namespace gui::species_tracking

#endif // KGD_SPECIESTRACKING_H
