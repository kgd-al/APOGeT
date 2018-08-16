#ifndef PHYLOGENYVIEWER_H
#define PHYLOGENYVIEWER_H

#include <QDialog>
#include <QGraphicsView>

#include "../core/phylogenictree.hpp"

namespace gui {

template <typename GENOME>
class PhylogenyViewer : public QDialog {
//  Q_OBJECT
public:
  using PTree = genotype::PhylogenicTree<GENOME>;

  struct Config {
    uint minSurvival;
    float minEnveloppe;
    bool showNames;
    bool circular;
    bool autofit;
  };

  PhylogenyViewer(QWidget *parent, PTree &ptree) : PhylogenyViewer(parent, ptree, defaultConfig()) {}
  PhylogenyViewer(QWidget *parent, PTree &ptree, Config config);

  void update (void);

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

private:
  const PTree &_ptree;
  Config _config;

  QGraphicsScene *_scene;
  QGraphicsView *_view;
};

} // end namespace visu

#endif // PHYLOGENYVIEWER_H
