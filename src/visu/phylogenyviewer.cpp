#include <QVBoxLayout>
#include <QGraphicsItem>
#include <QStack>

#include <QToolBar>
#include <QSlider>
#include <QCheckBox>

#include <QStyle>
#include <QStyleOptionSlider>
#include <QToolTip>

#include <QPixmap>
#include <QFileDialog>

#include <QDebug>

#include "phylogenyviewer.h"
#include "graphicsviewzoom.h"

namespace gui {
struct FancySlider : public QSlider {
  FancySlider(Qt::Orientation orientation, QWidget * parent = 0) : QSlider(orientation, parent) {}

  void sliderChange(QAbstractSlider::SliderChange change) {
    QSlider::sliderChange(change);

    if (change == QAbstractSlider::SliderValueChange) {
      QStyleOptionSlider opt;
      initStyleOption(&opt);

      QRect sr = style()->subControlRect(QStyle::CC_Slider, &opt, QStyle::SC_SliderHandle, this);
      QPoint bottomRightCorner = sr.bottomLeft();

      QToolTip::showText(mapToGlobal( QPoint( bottomRightCorner.x(), bottomRightCorner.y() ) ), QString::number(value()), this);
    }
  }
};

template <typename O, typename F>
auto makeSlider (int min, int max, O *object, F callback) {
  QSlider *slider = new FancySlider(Qt::Horizontal);
  slider->setMinimum(min);
  slider->setMaximum(max);
  QObject::connect(slider, &QSlider::valueChanged, object, callback);

  return slider;
}

void PhylogenyViewer_base::constructorDelegate(uint steps) {
  _items = {new QGraphicsScene(this), nullptr, nullptr, {}};

  _view = new QGraphicsView(_items.scene, this);
  _items.scene->setBackgroundBrush(Qt::transparent);
  _view->setRenderHint(QPainter::Antialiasing, true);
//  _view->viewport()->setMinimumSize(512, 512);
  _view->setDragMode(QGraphicsView::ScrollHandDrag);
  _view->setBackgroundBrush(Qt::white);

  QVBoxLayout *layout = new QVBoxLayout;

  QToolBar *toolbar = new QToolBar;

  QSlider *mSSlider = makeSlider(0, steps, this, &PhylogenyViewer_base::updateMinSurvival);
  connect(this, &PhylogenyViewer_base::onTreeStepped, [mSSlider] (uint step, const auto &) {
    mSSlider->setMaximum(step);
  });

  QSlider *mESlider = makeSlider(0, 100, this, &PhylogenyViewer_base::updateMinEnveloppe);

  QAction *print = new QAction(style()->standardPixmap(QStyle::SP_DialogSaveButton), "Print", this);
  print->setShortcut(Qt::ControlModifier + Qt::Key_P);
  connect(print, &QAction::triggered, this, &PhylogenyViewer_base::print);

  QCheckBox *showNames = new QCheckBox("Show names");
  showNames->setChecked(_config.showNames);
  connect(showNames, &QCheckBox::toggled, this, &PhylogenyViewer_base::toggleShowNames);

  QCheckBox *autofit = new QCheckBox("AutoFit");
  autofit->setChecked(_config.autofit);
  connect(autofit, &QCheckBox::toggled, this, &PhylogenyViewer_base::makeFit);

  toolbar->addWidget(mSSlider);
  toolbar->addWidget(mESlider);
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
  printTo(filename);
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
void PhylogenyViewer_base::treeStepped (uint step, const std::set<uint> &living) {
  for (Node *n: _items.nodes) {
    bool isAlive = (living.find(n->id) != living.end());
    n->updateNode(isAlive, step);
  }
  _items.border->setHeight(step);
  _items.scene->setSceneRect(_items.border->boundingRect());
  makeFit(_config.autofit);
}

void PhylogenyViewer_base::genomeEntersEnveloppe (uint sid, uint) {
  const uint K = PTreeConfig::enveloppeSize();
  Node *n = _items.nodes.value(sid);
  n->enveloppe = std::min(n->enveloppe + 1, K);
  n->autoscale();
}

void PhylogenyViewer_base::genomeLeavesEnveloppe (uint, uint) {}


void PhylogenyViewer_base::printTo (QString filename) {
  if (filename.isEmpty())
    filename = QFileDialog::getSaveFileName(this, "Save to", ".", "Images (*.png)");

  if (filename.isEmpty()) return;

  int Scale = 2;

  QPixmap pixmap (Scale * _items.scene->sceneRect().size().toSize());
  pixmap.fill(Qt::transparent);

  QPainter painter(&pixmap);
  painter.setRenderHint(QPainter::Antialiasing);
  _items.scene->render(&painter);
  painter.end();

  pixmap.save(filename);
  std::cout << "Saved to " << filename.toStdString() << std::endl;
}

}
