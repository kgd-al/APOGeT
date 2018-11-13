#include <QGraphicsItem>
#include <QStack>

#include <QToolBar>
#include <QLabel>
#include <QSlider>
#include <QCheckBox>

#include <QStyle>
#include <QStylePainter>
#include <QStyleOptionSlider>
#include <QToolTip>

#include <QPixmap>
#include <QFileDialog>

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

/// Wrapper for a slider with a floating current value indicator
struct FancySlider : public QSlider {

  /// The prefix used for displaying the current value
  const QString label;

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
    new QGraphicsScene(this), nullptr, nullptr, nullptr, {},
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
  Qt::Orientation orientation;
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
  QSlider *mSSlider = makeSlider(orientation, "Min. survival",
                                 0, steps, this, &PhylogenyViewer_base::updateMinSurvival);
  connect(this, &PhylogenyViewer_base::onTreeStepped, [mSSlider] (uint step, const auto &) {
    mSSlider->setMaximum(step);
  });

  // Enveloppe slider
  QSlider *mESlider = makeSlider(orientation, "Min. enveloppe", 0, 100, this, &PhylogenyViewer_base::updateMinEnveloppe);

  // Print action
  QAction *print = new QAction(style()->standardPixmap(QStyle::SP_DialogSaveButton), "Print", this);
  print->setShortcut(Qt::ControlModifier + Qt::Key_P);
  connect(print, &QAction::triggered, this, &PhylogenyViewer_base::render);

  // Show names checkbox
  QCheckBox *showNames = new QCheckBox("Show names");
  showNames->setChecked(_config.showNames);
  connect(showNames, &QCheckBox::toggled, this, &PhylogenyViewer_base::toggleShowNames);

  // Autofit checkbox
  QCheckBox *autofit = new QCheckBox("AutoFit");
  autofit->setChecked(_config.autofit);
  connect(autofit, &QCheckBox::toggled, this, &PhylogenyViewer_base::makeFit);

  // Position widgets with respect to the orientation
  if (orientation == Qt::Horizontal) {
    toolbar->addWidget(mSSlider);
    toolbar->addWidget(mESlider);
  } else {
    auto *holder = new QWidget;
    auto *layout = new QHBoxLayout;
    layout->addWidget(mSSlider);
    layout->addWidget(mESlider);
    holder->setLayout(layout);
    toolbar->addWidget(holder);
  }

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
}

void PhylogenyViewer_base::genomeEntersEnveloppe (SID sid, GID) {
  const uint K = config::PTree::enveloppeSize();
  Node *n = _items.nodes.value(sid);
  n->enveloppe = std::min(n->enveloppe + 1, K);
  n->autoscale();
}

void PhylogenyViewer_base::genomeLeavesEnveloppe (SID, GID) {}

void PhylogenyViewer_base::majorContributorChanged(SID sid, SID oldMC, SID newMC) {
  assert(false);
}


// ============================================================================
// == User requests
// ============================================================================

void PhylogenyViewer_base::renderTo (QString filename) {
  if (filename.isEmpty())
    filename = QFileDialog::getSaveFileName(this, "Save to", ".", "Images (*.png)");

  if (filename.isEmpty()) return;

  QPixmap pixmap = renderToPixmap();
  pixmap.save(filename);
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

}
