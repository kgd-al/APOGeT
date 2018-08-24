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

  virtual void update (bool save = false) = 0;
  void update (uint step, bool save);

  void resizeEvent(QResizeEvent *event);

  static Config defaultConfig (void) {
    return { 0, 0, true, true, true };
  }

signals:
  void updatedMaxSurvival (int v);

  void newSpecies (uint sid);
  void genomeEntersEnveloppe (uint sid, uint gid);
  void genomeLeavesEnveloppe (uint sid, uint gid);

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
  using PTCallbacks = typename PTree::Callbacks;
  friend PTCallbacks;

  using PTI = phylogeny::PTreeIntrospecter<GENOME>;

  PhylogenyViewer(QWidget *parent, PTree &ptree) : PhylogenyViewer(parent, ptree, defaultConfig()) {}
  PhylogenyViewer(QWidget *parent, PTree &ptree, Config config)
    : PhylogenyViewer_base(parent, config), _ptree(ptree), callbacks(this) {
    _ptree.setCallbacks(&callbacks);
    constructorDelegate(PTI::totalSteps(_ptree));
  }

  void update (bool save = false) override {
    uint step = PTI::totalSteps(_ptree);
    _scene->clear();

    if (step >= 0) {
      PTI::fillScene(_ptree, _scene, _config);
      if (_config.circular)
        PTI::toCircular(_scene);

      makeFit(_config.autofit);
    }

    PhylogenyViewer_base::update(step, save);
  }

private:
  const PTree &_ptree;
  PTCallbacks callbacks;
};

} // end namespace gui

template <typename GENOME>
struct phylogeny::Callbacks_t<phylogeny::PhylogenicTree<GENOME>> {
  Callbacks_t (gui::PhylogenyViewer_base *v) : viewer(v) {}

  void onNewSpecies (uint sid) {
    emit viewer->newSpecies(sid);
  }

  void onGenomeEntersEnveloppe (uint sid, uint gid) {
    emit viewer->genomeEntersEnveloppe(sid, gid);
  }

  void onGenomeLeavesEnveloppe (uint sid, uint gid) {
    emit viewer->genomeLeavesEnveloppe(sid, gid);
  }

private:
  gui::PhylogenyViewer_base *viewer;
};

#endif // PHYLOGENYVIEWER_H
