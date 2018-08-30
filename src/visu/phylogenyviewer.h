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

  void render (uint step);

  void resizeEvent(QResizeEvent *event);

  static Config defaultConfig (void) {
    return { 0, 0, true, true, true };
  }

signals:
  void onTreeStepped (uint step, const std::set<uint> &living);
  void onNewSpecies (uint pid, uint sid);
  void onGenomeEntersEnveloppe (uint sid, uint gid);
  void onGenomeLeavesEnveloppe (uint sid, uint gid);

public slots:
  // Callbacks from PTree
  void treeStepped (uint step, const std::set<uint> &living);
//  void newSpecies (uint pid, uint sid);
  void genomeEntersEnveloppe (uint sid, uint gid);
  void genomeLeavesEnveloppe (uint sid, uint gid);

  // Config update
  void updateMinSurvival (uint v);
  void updateMinEnveloppe (int v);
  void toggleShowNames (void);
  void makeFit (bool autofit);

  void print (void) {
    printTo("");
  }
  void printTo (QString filename);

protected:
  Config _config;

  GUIItems _items;
  QGraphicsView *_view;

  void constructorDelegate (uint steps);

  template <typename F>
  void updateNodes (F f) {
    for (Node *n: _items.nodes) f(n);
  }

  template <typename F>
  void updateRoot (F f) {
    f(_items.root);
  }

  virtual void updateLayout (void) = 0;
};

template <typename GENOME>
class PhylogenyViewer : public PhylogenyViewer_base {
public:
  using PTree = phylogeny::PhylogenicTree<GENOME>;
  using PTCallbacks = typename PTree::Callbacks;
  friend PTCallbacks;

  using Builder = PTGraphBuilder;

  PhylogenyViewer(QWidget *parent, PTree &ptree) : PhylogenyViewer(parent, ptree, defaultConfig()) {}
  PhylogenyViewer(QWidget *parent, PTree &ptree, Config config)
    : PhylogenyViewer_base(parent, config), _ptree(ptree), callbacks(this) {
    _ptree.setCallbacks(&callbacks);

    uint step = _ptree.step();
    constructorDelegate(step);

    build();
  }

  auto cache (void) {
    return PTreeBuildingCache { _config, _ptree.step(), _items };
  }

  void build (void) {
    auto c = cache();
    Builder::fillScene(_ptree, c);
    Builder::updateLayout(_items);
    makeFit(_config.autofit);
  }

  void render (void) {
    PhylogenyViewer_base::render(_ptree.step());
  }

public /*slots*/:
  void newSpecies (uint pid, uint sid) {
    auto c = cache();
    Node *parent = (pid != PTree::NoID) ? _items.nodes[pid] : nullptr;
    const auto &pn = *_ptree.nodeAt(sid);
    Builder::addSpecies(parent, pn, c);
    Builder::updateLayout(_items);
    _items.border->setEmpty(false);
    _view->update();
  }

private:
  const PTree &_ptree;
  PTCallbacks callbacks;

  void updateLayout (void) override {
    Builder::updateLayout(_items);
  }
};

} // end namespace gui

template <typename GENOME>
struct phylogeny::Callbacks_t<phylogeny::PhylogenicTree<GENOME>> {
  using PT = phylogeny::PhylogenicTree<GENOME>;
  using PV = gui::PhylogenyViewer<GENOME>;

  using GID = typename PT::GID;
  using SID = typename PT::SID;
  using LivingSet = typename PT::LivingSet;

  Callbacks_t (PV *v) : viewer(v) {}

  void onStepped (uint step, const LivingSet &living) {
    viewer->treeStepped(step, living);
    emit viewer->onTreeStepped(step, living);
  }

  void onNewSpecies (SID pid, SID sid) {
    viewer->newSpecies(pid, sid);
    emit viewer->onNewSpecies(pid, sid);
  }

  void onGenomeEntersEnveloppe (SID sid, GID gid) {
    viewer->genomeEntersEnveloppe(sid, gid);
    emit viewer->onGenomeEntersEnveloppe(sid, gid);
  }

  void onGenomeLeavesEnveloppe (SID sid, GID gid) {
    viewer->genomeLeavesEnveloppe(sid, gid);
    emit viewer->onGenomeLeavesEnveloppe(sid, gid);
  }

private:
  PV *viewer;
};

#endif // PHYLOGENYVIEWER_H
