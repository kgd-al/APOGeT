#include <QGraphicsItem>
#include <QGraphicsSceneMouseEvent>
#include <QStack>

#include <QToolBar>
#include <QLabel>
#include <QSlider>
#include <QCheckBox>
#include <QListWidget>

#include <QStyle>
#include <QStylePainter>
#include <QStyleOptionSlider>
#include <QToolTip>

#include <QPixmap>
#include <QFileDialog>

#ifndef NO_PRINTER
#include <QPrinter>
#endif

#ifndef NO_SVG
#include <QtSvg/QSvgGenerator>
#endif

#include <QDebug>


#include "phylogenyviewer.h"
#include "graphicutils.h"
#include "graphicsviewzoom.h"

/*!
 * \file phylogenyviewer.cpp
 *
 * Contains the implementation of the base phylogeny viewer
 */

namespace gui {

class OrientableLabel : public QLabel {
  Qt::Orientation _orientation;

public:
  template <typename... ARGS>
  OrientableLabel (Qt::Orientation o, ARGS&&... args)
    : QLabel(std::forward<ARGS>(args)...), _orientation(o) {}

  QSize minimumSizeHint(void) const override {
    QSize baseMinSize = QLabel::minimumSizeHint();
    if (_orientation == Qt::Vertical)
          return QSize(baseMinSize.height(), baseMinSize.width());
    else  return baseMinSize;
  }

  QSize sizeHint(void) const override {
    QSize baseSize = QLabel::sizeHint();
    if (_orientation == Qt::Vertical)
          return QSize(baseSize.height(), baseSize.width());
    else  return baseSize;
  }

  void setOrientation (Qt::Orientation o) {
    _orientation = o;
  }

  Qt::Orientation orientation (void) const {
    return _orientation;
  }

  void paintEvent(QPaintEvent *e) override {
    if (_orientation == Qt::Vertical) {
      QStylePainter painter(this);

      QStyleOptionFocusRect option;
      option.initFrom(this);
      option.backgroundColor = palette().color(QPalette::Background);

      painter.translate(0,sizeHint().height());
      painter.rotate(270);

      painter.drawText(QRect(QPoint(0, 0), QLabel::sizeHint()), Qt::AlignCenter, text());

    } else
      QLabel::paintEvent(e);
  }
};

/// Wrapper for a slider with a floating current value indicator
class FancySlider : public QSlider {

  /// The prefix used for displaying the current value
  const QString label;

public:
  /// Create a fancy slider with given \p orientation and prefix \p label
  FancySlider(Qt::Orientation orientation, const QString &label)
    : QSlider(orientation), label(label) {
    toolTip();
  }

  /// Update the current tooltip based on the #label and slider's current value
  QString toolTip(void) {
    QString tooltip = label + ": " + QString::number(value());
    setToolTip(tooltip);
    return tooltip;
  }

  /// Intercepted to provide a floating tooltip displaying the current value
  void sliderChange(QAbstractSlider::SliderChange change) override {
    QSlider::sliderChange(change);

    if (change == QAbstractSlider::SliderValueChange) {
      QStyleOptionSlider opt;
      initStyleOption(&opt);

      QRect sr = style()->subControlRect(QStyle::CC_Slider, &opt, QStyle::SC_SliderHandle, this);
      QPoint bottomRightCorner = sr.bottomLeft();

      QToolTip::showText(mapToGlobal( QPoint( bottomRightCorner.x(), bottomRightCorner.y() ) ),
                         toolTip(), this);
    }
  }
};

/// Helper template generating a parametered slider in one neat one-liner
template <typename O, typename F>
auto makeSlider (Qt::Orientation orientation, const QString label,
                 int min, int max, O *object, F callback) {
  QSlider *slider = new FancySlider(orientation, label);
  slider->setMinimum(min);
  slider->setMaximum(max);
  QObject::connect(slider, &QSlider::valueChanged, object, callback);

  return slider;
}

void PhylogenyViewer_base::constructorDelegate(uint steps, Direction direction) {
  // Create cache
  _items = {
    false,
    new QGraphicsScene(this),
    nullptr, nullptr, nullptr, nullptr,
    {},
    PTGraphBuilder::buildPenSet()
  };

  // Create view
  _view = new QGraphicsView(_items.scene, this);
  _items.scene->setBackgroundBrush(Qt::transparent);
  _view->setRenderHint(QPainter::Antialiasing, true);
  _view->setDragMode(QGraphicsView::ScrollHandDrag);
  _view->setBackgroundBrush(Qt::white);

  // Create layout
  auto *layout = new QBoxLayout(direction);
  Qt::Orientation orientation = Qt::Vertical;
  switch (direction) {
  case Direction::LeftToRight:
  case Direction::RightToLeft:
    orientation = Qt::Vertical;
    break;

  case Direction::TopToBottom:
  case Direction::BottomToTop:
    orientation = Qt::Horizontal;
    break;
  }

  // Create components
  QToolBar *toolbar = new QToolBar();
  toolbar->setOrientation(orientation);

  // Survival slider
  QString mSText = "Min. Survival";
  OrientableLabel *mSLabel = new OrientableLabel (orientation, mSText);
  QSlider *mSSlider = makeSlider(orientation, mSText,
                                 0, steps, this, &PhylogenyViewer_base::updateMinSurvival);
  mSSlider->setValue(_config.minSurvival);

  // Enveloppe slider
  QString mEText = "Min. enveloppe";
  OrientableLabel *mELabel = new OrientableLabel (orientation, mEText);
  QSlider *mESlider = makeSlider(orientation, mEText, 0, 100,
                                 this, &PhylogenyViewer_base::updateMinEnveloppe);
  mESlider->setValue(100 * _config.minEnveloppe);

  // Enveloppe slider
  QString cAText = "Clip at";
  OrientableLabel *cALabel = new OrientableLabel (orientation, cAText);
  QSlider *cASlider = makeSlider(orientation, cAText, 0, steps,
                                 this, &PhylogenyViewer_base::updateClippingRange);

  if (_config.clippingRange == uint(-1))
    cASlider->setValue(steps);
  else
    cASlider->setValue(_config.clippingRange);

  connect(this, &PhylogenyViewer_base::onTreeStepped,
          [mSSlider, cASlider] (uint step, const auto &) {
    mSSlider->setMaximum(step);

    bool wasAtMax = (cASlider->value() == cASlider->maximum());
    cASlider->setMaximum(step);
    if (wasAtMax) cASlider->setValue(step);
  });

  // Print action
  QAction *print = new QAction(style()->standardPixmap(QStyle::SP_DialogSaveButton), "Print", this);
  print->setShortcut(Qt::ControlModifier + Qt::Key_P);
  connect(print, &QAction::triggered, [this] { renderTo(""); });

  // Show names checkbox
  QCheckBox *showNames = new QCheckBox("Show names");
  showNames->setChecked(_config.showNames);
  connect(showNames, &QCheckBox::toggled, this, &PhylogenyViewer_base::toggleShowNames);

  // Autofit checkbox
  QCheckBox *autofit = new QCheckBox("AutoFit");
  autofit->setChecked(_config.autofit);
  connect(autofit, &QCheckBox::toggled, this, &PhylogenyViewer_base::makeFit);

  // Position widgets with respect to the orientation
  auto *sliderHolder = new QWidget;
  auto *sliderLayout = new QGridLayout;
  if (orientation == Qt::Horizontal) {
    sliderLayout->addWidget( mSLabel, 0, 0, Qt::AlignRight);
    sliderLayout->addWidget(mSSlider, 0, 1);
    sliderLayout->addWidget( mELabel, 1, 0, Qt::AlignRight);
    sliderLayout->addWidget(mESlider, 1, 1);
    sliderLayout->addWidget( cALabel, 2, 0, Qt::AlignRight);
    sliderLayout->addWidget(cASlider, 2, 1);

  } else {
    sliderLayout->addWidget( mSLabel, 0, 0, Qt::AlignBottom);
    sliderLayout->addWidget(mSSlider, 1, 0);
    sliderLayout->addWidget( mELabel, 0, 1, Qt::AlignBottom);
    sliderLayout->addWidget(mESlider, 1, 1);
    sliderLayout->addWidget( cALabel, 0, 2, Qt::AlignBottom);
    sliderLayout->addWidget(cASlider, 1, 2);
  }
  sliderHolder->setLayout(sliderLayout);
  toolbar->addWidget(sliderHolder);

  // Populate the rest of the toolbar
  toolbar->addWidget(showNames);
  toolbar->addWidget(autofit);
  toolbar->addAction(print);

  layout->addWidget(_view);
  layout->addWidget(toolbar);
  setLayout(layout);

  auto gvz = new Graphics_view_zoom(_view);
  connect(gvz, &Graphics_view_zoom::zoomed, [autofit] {
    autofit->setChecked(false);
  });

  setWindowTitle("Phenotypic tree");
}

void PhylogenyViewer_base::render(uint step) {
  QString filename;
  QTextStream qss(&filename);
  qss << "snapshots/ptree_step" << step << ".png";
  static int i=0;
  qDebug() << "[" << i++ << "] saved " << filename;
  renderTo(filename);
}

void PhylogenyViewer_base::resizeEvent(QResizeEvent*) {
  makeFit(_config.autofit);
}


// ============================================================================
// == Config update
// ============================================================================

void PhylogenyViewer_base::updateMinSurvival(uint v) {
  _config.minSurvival = v;
  updateNodes([this] (Node *n) {
    n->setVisible(Node::MIN_SURVIVAL, n->survival() >= _config.minSurvival);
  });
  updateLayout();
}

void PhylogenyViewer_base::updateMinEnveloppe(int v) {
  _config.minEnveloppe = v / 100.;
  updateNodes([this] (Node *n) {
    n->setVisible(Node::MIN_FULLNESS, n->fullness() >= _config.minEnveloppe);
  });
  updateLayout();
}

void PhylogenyViewer_base::updateClippingRange (uint t) {
  _config.clippingRange = t;
  updateNodes([this] (Node *n) {
    n->updateNode(_config.clippingRange <= n->disappearance());
    n->setVisible(Node::CLIP_RANGE, n->appearance() <= _config.clippingRange);
  });

  if (_items.initialized) {
    QPainterPath dimPath;
    float r = radius();
    if (t != r) {
      dimPath.addEllipse(QPointF(0, 0), r, r);
      dimPath.addEllipse(QPointF(0, 0), t, t);
    }
    _items.dimmer->setDimmingPath(dimPath);
  }

  updateLayout();
}


void PhylogenyViewer_base::toggleShowNames(void) {
  _config.showNames = !_config.showNames;
  updateNodes([this] (Node *n) {
    n->setVisible(Node::SHOW_NAME, _config.showNames);
  });
}

void PhylogenyViewer_base::makeFit(bool autofit) {
  _config.autofit = autofit;
  if (_config.autofit)  _view->fitInView(_items.scene->sceneRect(), Qt::KeepAspectRatio);
}


// ============================================================================
// == Phylogeny update
// ============================================================================

void PhylogenyViewer_base::treeStepped (uint step, const LivingSet &living) {
  updatePens();

  for (Node *n: _items.nodes) n->updateNode(living.find(n->id) != living.end());
  _items.border->setRadius(step);
  _items.scene->setSceneRect(_items.border->boundingRect());
  makeFit(_config.autofit);

  if (_config.screenshots) {
    QString filename;
    QTextStream qss(&filename);
    qss << "snapshots/ptree_step" << step << ".png";
    renderTo(filename);
  }
}

void PhylogenyViewer_base::genomeEntersEnveloppe (SID sid, GID) {
  const uint K = config::PTree::enveloppeSize();
  Node *n = _items.nodes.value(sid);
  n->enveloppe = std::min(n->enveloppe + 1, K);
  n->autoscale();
}

void PhylogenyViewer_base::genomeLeavesEnveloppe (SID, GID) {}

void PhylogenyViewer_base::majorContributorChanged(SID sid, SID oldMC, SID newMC) {
  Node *n = _items.nodes.value(sid),
       *oldP = _items.nodes.value(oldMC),
       *newP = _items.nodes.value(newMC);

  assert(n && oldP && newP);
  assert(n->id == sid);
  assert(oldP->id == oldMC);
  assert(newP->id == newMC);

  n->parent = newP;
  oldP->subnodes.removeAll(n);
  newP->subnodes.append(n);
  std::sort(newP->subnodes.begin(), newP->subnodes.end(),
            [] (const Node *lhs, const Node *rhs) {
    return lhs->id > rhs->id;
  });

  if (n->path)  n->path->start = newP;
  n->setVisible(Node::PARENT, newP->subtreeVisible());

  PTGraphBuilder::updateLayout(_items);
  _view->update();

  qDebug() << "Major contributor for species" << uint(sid)
           << "changed from" << uint(oldMC) << "to" << uint(newMC);
}


// ============================================================================
// == User requests
// ============================================================================

void PhylogenyViewer_base::speciesDetailPopup (SID id, QStringList data,
                                               QGraphicsSceneMouseEvent *e) {
  QDialog *dialog = new QDialog (this);
    QVBoxLayout *layout = new QVBoxLayout;
      QLabel *generalLabel = new QLabel (data.takeFirst());
      QListWidget *listLabel = new QListWidget;

    layout->addWidget(generalLabel);
      listLabel->setFlow(QListWidget::LeftToRight);
      listLabel->insertItems(0, data);
    layout->addWidget(listLabel);
  dialog->setLayout(layout);

  dialog->setWindowTitle(QString("Details of species ")
                       + QString::number(std::underlying_type<SID>::type(id)));
  dialog->move(e->screenPos());
  dialog->show();
}

void PhylogenyViewer_base::renderTo (QString filename) {
  Node *hovered = nullptr;
  if (filename.isEmpty()) {
    hovered = _items.contributors->species;
    filename = QFileDialog::getSaveFileName(this, "Save to", ".",
                                            "PDF (*.pdf);; Images (*svg,*.png)");
    if (hovered)  hovered->hoverEnterEvent(nullptr);
  }

  if (filename.isEmpty()) return;

  QString ext = filename.mid(filename.lastIndexOf('.')+1);
  if (ext == "pdf")
    renderToPDF(filename);

  else if (ext == "svg")
    renderToSVG(filename);

  else {
    if (ext != "png")
      std::cout << "Unkown extension type '" << ext.toStdString()
                << "'. Defaulting to png" << std::endl;
    QPixmap pixmap = renderToPixmap();
    pixmap.save(filename);
  }

  if (hovered)  hovered->hoverLeaveEvent(nullptr);
  std::cout << "Saved to " << filename.toStdString() << std::endl;
}

QPixmap PhylogenyViewer_base::renderToPixmap (const QSize &requestedSize) const {
  QPixmap pixmap (requestedSize);
  pixmap.fill(Qt::transparent);
  QRectF bounds = gui::centeredInto(QRectF({0,0}, requestedSize),
                                    _items.scene->sceneRect());

  QPainter painter(&pixmap);
  painter.setRenderHint(QPainter::Antialiasing);
  _items.scene->render(&painter, bounds.toRect());
  painter.end();

  return pixmap;
}

#ifndef NO_PRINTER
void PhylogenyViewer_base::renderToPDF(const QString &filename) const {
  QPrinter printer(QPrinter::HighResolution);
  printer.setPageSizeMM(_items.scene->sceneRect().size());
  printer.setOrientation(QPrinter::Portrait);
  printer.setOutputFormat(QPrinter::PdfFormat);
  printer.setOutputFileName(filename);

  QPainter p;
  if(!p.begin( &printer ))
    std::cerr << "Failed to creating painting device for file "
              << filename.toStdString() << std::endl;

  else
    _items.scene->render(&p);
}
#endif

#ifndef NO_SVG
void PhylogenyViewer_base::renderToSVG(const QString &filename) const {
  static const int S = 500;
  QSvgGenerator generator;
  generator.setFileName(filename);
  generator.setSize(QSize(S, S));
  generator.setViewBox(QRect(0, 0, S, S));
  generator.setTitle("Phylogenic tree");
  generator.setDescription("Idem");

  QPainter p;
  if(!p.begin(&generator)) {
    std::cerr << "Failed to creating painting device for file "
              << filename.toStdString() << std::endl;
    return;
  }

  _items.scene->render(&p);
  p.end();
}
#endif

}
