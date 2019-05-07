#ifndef SPECIESTRACKING_H
#define SPECIESTRACKING_H

#include <QStyledItemDelegate>
#include <QComboBox>
#include <QTableView>
#include <QStandardItemModel>
#include <QDialogButtonBox>

#include "ptgraphbuilder.h"

namespace gui::species_tracking {
/// \cond internal
using Specs = decltype(ViewerConfig::colorSpecs);
using SID = decltype(ViewerConfig::ColorSpec::sid);
using Color = decltype(ViewerConfig::ColorSpec::color);

inline QString toString (SID sid) {
  return QString::number(uint(sid));
}

class SIDDelegate : public QStyledItemDelegate {
  using Data = std::set<SID>;
  Data data;
  int longest;

public:
  SIDDelegate(Data &&data, QWidget *parent = 0);

  auto nextSID (void) const {
    return *data.begin();
  }

  QSize sizeHint(const QStyleOptionViewItem &/*option*/,
                 const QModelIndex &/*index*/) const override;

  QWidget *createEditor(QWidget *parent, const QStyleOptionViewItem &/*option*/,
                        const QModelIndex &/*index*/) const override;

  void setEditorData(QWidget *editor, const QModelIndex &index) const override;

  void setModelData(QWidget *editor, QAbstractItemModel *model,
                    const QModelIndex &index) const override;
};

class ColorDelegate : public QStyledItemDelegate {
  static const QList<QColor> defaultColors;
public:
  ColorDelegate(QWidget *parent = 0) : QStyledItemDelegate(parent) {}

  static auto nextColor (int i) {
    return defaultColors.at(i % defaultColors.size());
  }

  void paint(QPainter *painter, const QStyleOptionViewItem &option,
             const QModelIndex &index) const override;

  QSize sizeHint(const QStyleOptionViewItem &/*option*/,
                 const QModelIndex &/*index*/) const override {
    return QSize (20, 10);
  }

  QWidget *createEditor(QWidget *parent, const QStyleOptionViewItem &/*option*/,
                        const QModelIndex &/*index*/) const override;

  void setEditorData(QWidget *editor, const QModelIndex &index) const override;

  void setModelData(QWidget *editor, QAbstractItemModel *model,
                    const QModelIndex &index) const override;
};

class TightTableView : public QTableView {
  QSize sizeHint (void) const;

public:
  /// \authors savolai & Ratah @https://stackoverflow.com/a/50217711
  void verticalFit (void);
};

class Dialog : public QDialog {
  Q_OBJECT

  QStandardItemModel model;
  SIDDelegate sidDelegate;
  ColorDelegate colorDelegate;

  QDialogButtonBox *controls;

  auto validSIDs (const PhylogenyViewer_base *viewer) const;

  auto columnFromSpec (const ViewerConfig::ColorSpec &cs);

  void dataChanged (void);

public:
  Dialog (PhylogenyViewer_base *viewer, const Specs &is);

  Specs colorSelection (void) const;

signals:
  void applied (void);
};

/// \endcond internal
} // end of namespace gui::species_tracking

#endif // SPECIESTRACKING_H
