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
  _scene = new QGraphicsScene(this);
  _view = new QGraphicsView(_scene, this);
  _scene->setBackgroundBrush(Qt::transparent);
  _view->setRenderHint(QPainter::Antialiasing, true);
//  _view->viewport()->setMinimumSize(512, 512);
  _view->setDragMode(QGraphicsView::ScrollHandDrag);
  _view->setBackgroundBrush(Qt::white);

  QVBoxLayout *layout = new QVBoxLayout;

  QToolBar *toolbar = new QToolBar;

  QSlider *mSSlider = makeSlider(0, steps, this, &PhylogenyViewer_base::updateMinSurvival);
  connect(this, &PhylogenyViewer_base::updatedMaxSurvival, mSSlider, &QSlider::setMaximum);

  QSlider *mESlider = makeSlider(0, 100, this, &PhylogenyViewer_base::updateMinEnveloppe);

  QAction *print = new QAction(style()->standardPixmap(QStyle::SP_DialogSaveButton), "Print", this);
  print->setShortcut(Qt::ControlModifier + Qt::Key_P);
  connect(print, &QAction::triggered, this, &PhylogenyViewer_base::print);

  QCheckBox *showNames = new QCheckBox("Show names");
  showNames->setChecked(_config.showNames);
  connect(showNames, &QCheckBox::toggled, this, &PhylogenyViewer_base::toggleShowNames);

  QCheckBox *circular = new QCheckBox("Circular");
  circular->setChecked(_config.circular);
  connect(circular, &QCheckBox::toggled, this, &PhylogenyViewer_base::toggleCircular);

  QCheckBox *autofit = new QCheckBox("AutoFit");
  autofit->setChecked(_config.autofit);
  connect(autofit, &QCheckBox::toggled, this, &PhylogenyViewer_base::makeFit);

  toolbar->addWidget(mSSlider);
  toolbar->addWidget(mESlider);
  toolbar->addWidget(showNames);
  toolbar->addWidget(circular);
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

void PhylogenyViewer_base::render(uint step, QString filename) {
  if (filename.isEmpty()) {
    QTextStream qss(&file);
    qss << "snapshots/ptree_step" << step << ".png";
    static int i=0;
    qDebug() << "[" << i++ << "] saved " << filename;
  }
  printTo(filename);
}

void PhylogenyViewer_base::resizeEvent(QResizeEvent*) {
  makeFit(_config.autofit);
}

void PhylogenyViewer_base::updateMinSurvival(int v) {
  _config.minSurvival = v;
  update([v] (QGraphicsItem *item) {
    if (Node *n = dynamic_cast<Node*>(item))
      n->setVisible(n->survival() >= v);
  });
}

void PhylogenyViewer_base::updateMinEnveloppe(int v) {
  _config.minEnveloppe = v / 100.;
  update();
}

void PhylogenyViewer_base::toggleShowNames(void) {
  _config.showNames = !_config.showNames;
  update();
}

void PhylogenyViewer_base::toggleCircular(void) {
  _config.circular = !_config.circular;
  update();
}

void PhylogenyViewer_base::makeFit(bool autofit) {
  _config.autofit = autofit;
  if (_config.autofit)  _view->fitInView(_scene->sceneRect(), Qt::KeepAspectRatio);
}

void PhylogenyViewer_base::printTo (QString filename) {
  if (filename.isEmpty())
    filename = QFileDialog::getSaveFileName(this, "Save to", ".", "Images (*.png)");

  if (filename.isEmpty()) return;

  int Scale = 2;

  QPixmap pixmap (Scale * _scene->sceneRect().size().toSize());
  pixmap.fill(Qt::transparent);

  QPainter painter(&pixmap);
  painter.setRenderHint(QPainter::Antialiasing);
  _scene->render(&painter);
  painter.end();

  pixmap.save(filename);
  std::cout << "Saved to " << filename.toStdString() << std::endl;
}

}
