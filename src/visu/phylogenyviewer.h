#ifndef PHYLOGENYVIEWER_H
#define PHYLOGENYVIEWER_H

#include <QDialog>
#include <QGraphicsView>

#include "ptreeintrospecter.hpp"

namespace gui {

class PhylogenyViewer_base : public QDialog {
  Q_OBJECT
public:
  struct Config {
    uint minSurvival;
    float minEnveloppe;
    bool showNames;
    bool circular;
    bool autofit;
  };

  PhylogenyViewer_base (QWidget *parent, Config config) : QDialog(parent), _config(config) {}

  virtual void update (void) = 0;
  void update (uint step);

  void resizeEvent(QResizeEvent *event);

  static Config defaultConfig (void) {
    return { 0, 0, true, true, true };
  }

public slots:
  void updateMinSurvival (int v);
  void updateMinEnveloppe (int v);
  void toggleShowNames (void);
  void toggleCircular (void);
  void makeFit (bool autofit);

  void print (void) {
    printTo("");
  }
  void printTo (QString filename);

protected:
  Config _config;

  QGraphicsScene *_scene;
  QGraphicsView *_view;

  void constructorDelegate (uint steps);
};

template <typename GENOME>
class PhylogenyViewer : public PhylogenyViewer_base {
public:
  using PTree = phylogeny::PhylogenicTree<GENOME>;
  using PTI = phylogeny::PTreeIntrospecter<GENOME>;

  PhylogenyViewer(QWidget *parent, PTree &ptree) : PhylogenyViewer(parent, ptree, defaultConfig()) {}
  PhylogenyViewer(QWidget *parent, PTree &ptree, Config config)
    : PhylogenyViewer_base(parent, config), _ptree(ptree) {
    constructorDelegate(PTI::totalSteps(_ptree));
  }

  void update (void) override {
    _scene->clear();

    if (PTI::totalSteps(_ptree) >= 0) {
      PTI::fillScene(_ptree, _scene, _config);
      if (_config.circular)
        PTI::toCircular(_scene);

      makeFit(_config.autofit);
    }

    PhylogenyViewer_base::update(PTI::totalSteps(_ptree));
  }

private:
  const PTree &_ptree;
};

} // end namespace gui

#endif // PHYLOGENYVIEWER_H
