#include <QDialog>
#include <QColorDialog>
#include <QHeaderView>
#include <QScrollBar>
#include <QPushButton>

#include <QDebug>

#include "speciestracking.h"

#include "phylogenyviewer.h"

namespace gui::species_tracking {

SIDDelegate::SIDDelegate(Data &&data, QWidget *parent)
  : QStyledItemDelegate(parent), data(data) {
  longest = 0;
  for (const auto &sid: data) {
    QString ssid = toString(sid);
    longest = std::max(longest, ssid.length());
  }
}

QSize SIDDelegate::sizeHint(const QStyleOptionViewItem &/*option*/,
               const QModelIndex &/*index*/) const {
  QComboBox cb;
  cb.addItem(QString(longest+1, ' '));
  return cb.sizeHint();
}

QWidget* SIDDelegate::createEditor(QWidget *parent,
                                   const QStyleOptionViewItem &/*option*/,
                                   const QModelIndex &/*index*/) const {
    QComboBox *cb = new QComboBox(parent);
    cb->setEditable(false);
    for (const auto &sid: data) cb->addItem(toString(sid));
    return cb;
  }

void SIDDelegate::setEditorData(QWidget *editor,
                                const QModelIndex &index) const {
  QComboBox *cb = qobject_cast<QComboBox *>(editor);
  Q_ASSERT(cb);
  // get the index of the text in the combobox that matches the current value of the item
  const QString currentText = index.data(Qt::EditRole).toString();
  const int cbIndex = cb->findText(currentText);
  // if it is valid, adjust the combobox
  if (cbIndex >= 0) cb->setCurrentIndex(cbIndex);
}

void SIDDelegate::setModelData(QWidget *editor, QAbstractItemModel *model,
                               const QModelIndex &index) const {
  QComboBox *cb = qobject_cast<QComboBox *>(editor);
  Q_ASSERT(cb);
  model->setData(index, cb->currentText(), Qt::EditRole);
}

void ColorDelegate::paint(QPainter *painter, const QStyleOptionViewItem &option,
           const QModelIndex &index) const {
  QStyleOptionButton button;
  button.rect = option.rect;
  button.state = QStyle::State_Enabled;
  button.state |= (option.state & QStyle::State_MouseOver);
  button.palette.setBrush(QPalette::Button,
                          index.data(Qt::BackgroundRole).value<QColor>());
  option.widget->style()->drawControl(QStyle::CE_PushButton, &button, painter);
}


QWidget* ColorDelegate::createEditor(QWidget *parent,
                                     const QStyleOptionViewItem &/*option*/,
                                     const QModelIndex &/*index*/) const {
  auto cd = new QColorDialog(parent);
  for (int i=0; i<defaultColors.size(); i++)
    cd->setCustomColor(i, defaultColors.at(i));
  return cd;
}

void ColorDelegate::setEditorData(QWidget *editor,
                                  const QModelIndex &index) const {
  QColorDialog *cd = qobject_cast<QColorDialog *>(editor);
  Q_ASSERT(cd);
  cd->setCurrentColor(index.data(Qt::BackgroundRole).value<Color>());
}

void ColorDelegate::setModelData(QWidget *editor, QAbstractItemModel *model,
                                 const QModelIndex &index) const {
  QColorDialog *cd = qobject_cast<QColorDialog *>(editor);
  Q_ASSERT(cd);
  model->setData(index, QVariant::fromValue(cd->currentColor()),
                 Qt::BackgroundRole);
}

const QList<QColor> ColorDelegate::defaultColors {
  QColor( 148,   0, 211 ),
  QColor(   0, 158, 115 ),
  QColor(  86, 180, 233 ),
  QColor( 230, 159,   0 ),
  QColor( 240, 228,  66 ),
  QColor(   0, 114, 178 ),
  QColor( 229,  30,  16 ),
  QColor(   0,   0,   0 ),
};

QSize TightTableView::sizeHint (void) const {
  int h = 0;
  for (int i=0; i<model()->rowCount(); i++) h += sizeHintForRow(i);
  return QSize(200, h);
}

void TightTableView::verticalFit (void) {
  static constexpr int margins = 2;
  int totalHeight = 0;

  int count = verticalHeader()->count();
  for (int i = 0; i < count; ++i)
    if (!verticalHeader()->isSectionHidden(i))
      totalHeight += verticalHeader()->sectionSize(i) + margins;

  if (!horizontalScrollBar()->isHidden())
    totalHeight += horizontalScrollBar()->height();

  if (!horizontalHeader()->isHidden())
    totalHeight += horizontalHeader()->height();

  setFixedHeight(totalHeight);
}

auto Dialog::validSIDs (const PhylogenyViewer_base *viewer) const {
  std::set<SID> visible;
  viewer->observeNodes([&visible] (const Node *n) {
    if (n->isVisible()) visible.insert(n->id);
  });
  return visible;
}

auto Dialog::columnFromSpec (const ViewerConfig::ColorSpec &cs) {
  QList<QStandardItem*> column;
  auto sid = new QStandardItem (toString(cs.sid));
  sid->setTextAlignment(Qt::AlignCenter);
  column << sid;
  auto color = new QStandardItem;
  color->setBackground(cs.color);
  column << color;
  auto enabled = new QStandardItem;
  enabled->setTextAlignment(Qt::AlignCenter); /// TODO This does nothing...
  enabled->setCheckable(true);
  enabled->setCheckState(cs.enabled ? Qt::Checked : Qt::Unchecked);
  column << enabled;
  return column;
}

Dialog::Dialog (PhylogenyViewer_base *viewer, const Specs &is)
  : QDialog(viewer), sidDelegate(validSIDs(viewer)) {

  setWindowTitle("Color Picker");

  auto *layout = new QVBoxLayout;
    auto hlayout = new QHBoxLayout;
      auto *tcontrols = new QDialogButtonBox(Qt::Vertical);
      auto *table = new TightTableView;
    controls = new QDialogButtonBox (
      QDialogButtonBox::Ok | QDialogButtonBox::Apply
    | QDialogButtonBox::Cancel);

  auto *add = tcontrols->addButton("+", QDialogButtonBox::ActionRole);
  connect(add, &QPushButton::clicked, [this, table] {
    model.appendColumn(columnFromSpec({
      sidDelegate.nextSID(),
      colorDelegate.nextColor(model.columnCount()),
      true
    }));
    table->resizeColumnToContents(model.columnCount() - 1);
    table->verticalFit();
    dataChanged();
  });

  auto *del = tcontrols->addButton("-", QDialogButtonBox::ActionRole);
  connect(del, &QPushButton::clicked, [this, table] {
    const auto selection = table->selectionModel()->selectedIndexes();
    for (const QModelIndex &i: selection) {
      if (i.row() == 0)
        model.removeColumn(i.column());
    }
    if (!selection.empty()) dataChanged();
  });

  setLayout(layout);
    layout->addLayout(hlayout);
      hlayout->addWidget(tcontrols);
      hlayout->addWidget(table);
    layout->addWidget(controls);

  tcontrols->setMaximumWidth(30);
  controls->button(QDialogButtonBox::Apply)->setEnabled(false);

  connect(controls, &QDialogButtonBox::accepted, this, &QDialog::accept);
  connect(controls, &QDialogButtonBox::rejected, this, &QDialog::reject);
  connect(controls->button(QDialogButtonBox::Apply), &QPushButton::clicked,
          [this] {
    controls->button(QDialogButtonBox::Apply)->setEnabled(false);
    emit applied();
  });

  for (const auto &cs: is)  model.appendColumn(columnFromSpec(cs));
  table->horizontalHeader()->hide();
  table->verticalHeader()->hide();
  table->setSizePolicy(QSizePolicy::Preferred, QSizePolicy::Fixed);
  table->setMinimumSize(100, 10);
  table->setModel(&model);
  table->setShowGrid(false);
  table->setItemDelegateForRow(0, &sidDelegate);
  table->setItemDelegateForRow(1, &colorDelegate);
  table->setHorizontalScrollBarPolicy(Qt::ScrollBarAlwaysOn);
  table->resizeColumnsToContents();
  table->resizeRowsToContents();
  table->verticalFit();
  connect(&model, &QAbstractItemModel::dataChanged, this, &Dialog::dataChanged);

  table->setSelectionMode(QAbstractItemView::SingleSelection);
  table->setSelectionBehavior(QAbstractItemView::SelectColumns);
  del->setEnabled(!table->selectionModel()->selectedColumns().empty());
  connect(table->selectionModel(), &QItemSelectionModel::selectionChanged,
          [del] (const QItemSelection &selected, const QItemSelection&) {
    del->setEnabled(!selected.empty());
  });
}

void Dialog::dataChanged(void) {
  controls->button(QDialogButtonBox::Apply)->setEnabled(true);
}

Specs Dialog::colorSelection(void) const {
  Specs specifications;
  for (int j=0; j<model.columnCount(); j++) {
    auto index = model.index(0, j);
    qDebug() << index.data() << ": "
             << index.sibling(1, j).data(Qt::BackgroundRole).value<Color>();
    specifications.insert({
      SID(index.data().toUInt()),
      model.index(1, j).data(Qt::BackgroundRole).value<Color>(),
      model.index(2, j).data(Qt::CheckStateRole).value<bool>()
    });
  }
  return specifications;
}
} // end of namespace gui::species_tracking
