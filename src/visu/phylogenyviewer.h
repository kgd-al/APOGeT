#ifndef PHYLOGENYVIEWER_H
#define PHYLOGENYVIEWER_H

#include <QDialog>
#include <QGraphicsView>

#include "ptgraphbuilder.h"

namespace gui {

class PhylogenyViewer_base : public QDialog {
  Q_OBJECT
protected:
  using Config = gui::ViewerConfig;
public:
  PhylogenyViewer_base (QWidget *parent, Config config) : QDialog(parent), _config(config) {}

  virtual void update (void) = 0;
  void render (uint step, QString filename);

  void resizeEvent(QResizeEvent *event);

  static Config defaultConfig (void) {
    return { 0, 0, true, true, true };
  }

signals:
  void treeStepped (uint step, const std::set<uint> &living);
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

    uint step = _ptree.step();
    constructorDelegate(step);
    build();
  }

  void build (void) {
    PTI::fillScene(_ptree, _scene, _config);
    makeFit(_config.autofit);
  }

  void render (QString filename = "") const {
    uint step = _ptree.step();
    PhylogenyViewer_base::render(step, filename);
  }

private:
  const PTree &_ptree;
  PTCallbacks callbacks;
};

} // end namespace gui

template <typename GENOME>
struct phylogeny::Callbacks_t<phylogeny::PhylogenicTree<GENOME>> {
  using PT = phylogeny::PhylogenicTree<GENOME>;
  using GID = typename PT::GID;
  using SID = typename PT::SID;
  using LivingSet = typename PT::LivingSet;

  Callbacks_t (gui::PhylogenyViewer_base *v) : viewer(v) {}

  void onStepped (uint step, const LivingSet &living) {
    emit viewer->treeStepped(step, living);
  }

  void onNewSpecies (SID sid) {
    emit viewer->newSpecies(sid);
  }

  void onGenomeEntersEnveloppe (SID sid, GID gid) {
    emit viewer->genomeEntersEnveloppe(sid, gid);
  }

  void onGenomeLeavesEnveloppe (SID sid, GID gid) {
    emit viewer->genomeLeavesEnveloppe(sid, gid);
  }

private:
  gui::PhylogenyViewer_base *viewer;
};

#endif // PHYLOGENYVIEWER_H
