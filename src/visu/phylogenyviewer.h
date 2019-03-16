#ifndef PHYLOGENYVIEWER_H
#define PHYLOGENYVIEWER_H

#include <QDialog>
#include <QGraphicsView>
#include <QBoxLayout>

#include "../core/tree/treetypes.h"

#include "ptgraphbuilder.h"
#include "layer.hpp"

/*!
 * \file phylogenyviewer.h
 *
 * Contains definition for the phylogeny top-level viewer
 */

namespace gui {

/// Base class for the phylogeny viewer. No template just the common functions,
/// signals and slots.
class PhylogenyViewer_base : public QDialog {
  Q_OBJECT
protected:
  /// \copydoc genotype::BOCData::GID
  using GID = genotype::BOCData::GID;

  /// Helper alias to the species identificator used in the phylogenic tree
  using SID = phylogeny::SID;

  /// Helper alias to the phylogenic tree's collection of living individuals
  using LivingSet = phylogeny::LivingSet;

  /// Configuration data controlling what to draw and how
  using Config = gui::ViewerConfig;

public:
  /// The direction in which to layout components
  using Direction = QBoxLayout::Direction;

  /// Create a phylogeny viewer with given \p parent and using \p config as
  /// its initial configuration
  PhylogenyViewer_base (QWidget *parent, Config config)
    : QDialog(parent), _config(config) {}

  /// Render the tree to file
  void render (uint step);

  /// Intercepted to allow for tree autofitting
  void resizeEvent(QResizeEvent *event);

  /// \returns the tree radius
  auto radius (void) const {
    return _items.border->radius;
  }

  /// \returns the tree bounding rectangle
  auto boundingRect (void) const {
    return _items.border->boundingRect();
  }

  /// \returns the pen of the requested type
  QPen pathPen (details::PenType t) {
    return _items.pens.value(t);
  }

  /// \return the default configuration
  static auto defaultConfig (void) {
    return Config{};
  }

signals:
  /// \brief Emitted when the tree is stepped
  /// \copydetails phylogeny::Callbacks_t::onStepped
  void onTreeStepped (uint step, const LivingSet &living);

  /// Emitted when a species has been added to the tree
  /// \copydetails phylogeny::Callbacks_t::onNewSpecies
  void onNewSpecies (SID pid, SID sid);

  /// Emitted when a species' enveloppe has changed
  /// \copydetails phylogeny::Callbacks_t::onGenomeEntersEnveloppe
  void onGenomeEntersEnveloppe (SID sid, GID gid);

  /// Emitted when a species has been added to the tree
  /// \copydetails phylogeny::Callbacks_t::onGenomeLeavesEnveloppe
  void onGenomeLeavesEnveloppe (SID sid, GID gid);

  /// Emitted when a species starts/stops being hovered
  void onSpeciesHoverEvent (SID sid, bool entered);

  /// Emitted when a species has changed its rooting point
  /// \copydetails phylogeny::Callbacks_t::onMajorContributorChanged
  void onMajorContributorChanged (SID sid, SID oldMC, SID newMC);


protected slots:
  // ===========================================================================
  // Callbacks from PTree

  /// Process a step event (new timestamp/living species)
  void treeStepped (uint step, const LivingSet &living);

  /// Process a enveloppe change event
  void genomeEntersEnveloppe (SID sid, GID gid);

  /// \copydoc genomeEntersEnveloppe
  void genomeLeavesEnveloppe (SID sid, GID gid);

  /// Process a rooting change event
  void majorContributorChanged(SID sid, SID oldMC, SID newMC);

  // ===========================================================================
  // Config update

  /// Called when the corresponding checkbox has been changed
  void toggleShowOnlySurvivors (void);

  /// Called when the corresponding slider's value has been changed
  void updateMinSurvival (uint v);

  /// \copydoc updateMinSurvival
  void updateMinEnveloppe (int v);

  /// \copydoc updateMinSurvival
  void updateClippingRange (uint t);

  /// \copydoc toggleShowOnlySurvivors
  void toggleShowNames (void);

public slots:
  /// Requests the scale of the view to be adapted to the size of the scene
  void makeFit (bool autofit);

  /// Pops a detailed view a species node contents up
  void speciesDetailPopup (SID id, QStringList data, const QString &summary,
                           QGraphicsSceneMouseEvent *e);

  /// Prints the current scene to the image file \p filename
  void renderTo (QString filename = "");

  /// Prints the current scene into a pixmap of size \p requestedSize
  QPixmap renderToPixmap (QSize requestedSize) const;

#ifndef NO_PRINTER
  /// Prints the current scene into a (vector?) portable document format file
  void renderToPDF (const QString &filename) const;
#endif

#ifndef NO_SVG
  /// Prints the current scene into scalable vector graphics file
  void renderToSVG (const QString &filename) const;
#endif

  /// Process a hover event
  virtual void hoverEvent (SID sid, bool entered) = 0;

  /// Process a double click event
  virtual void doubleClickEvent (const Node &n, QGraphicsSceneMouseEvent *e) = 0;

protected:
  /// The graphics config
  Config _config;

  /// Cache data for all graphics items managed by this viewer
  GUIItems _items;

  /// The view in which the graphics items reside
  QGraphicsView *_view;

  /// Constructor delegate called by template instantiations
  void constructorDelegate (uint steps,
                            Direction direction = Direction::LeftToRight);

  /// Ensure pens are consistent with the current state of the ptree
  void updatePens (void) {
    PTGraphBuilder::updatePenSet(radius(), _items.pens);
  }

  /// Apply function \p f to all of the scene's current nodes
  template <typename F>
  void updateNodes (F f) {
    for (Node *n: _items.nodes) f(n);
  }

  /// Apply function \p f starting from the scene's root node.
  /// \attention The recursivity (or lack thereof) is left to the caller's
  /// discretion
  template <typename F>
  void updateRoot (F f) {
    f(_items.root);
  }

  /// Called when the nodes/paths position must be changed radically (e.g.
  /// new species inserted or visibility criterion modified)
  virtual void updateLayout (void) = 0;
};

/// \brief Instantiable phylogeny viewer for type \p GENOME.
///
/// This class can be used as a standalone (top-level) viewer or be embedded in
/// any kind of layout. It is possible to control the orientation of the controls
/// with respect to the main view.
template <typename GENOME, typename UDATA>
class PhylogenyViewer : public PhylogenyViewer_base {
public:

  /// Helper alias to the PTree visualized
  using PTree = phylogeny::PhylogeneticTree<GENOME, UDATA>;

  /// Helper alias to the callbacks type used by the PTree
  using PTCallbacks = typename PTree::Callbacks;

  /// These need priviliged access
  friend PTCallbacks;

  /// Helper alias for the graph builder
  using Builder = PTGraphBuilder;

  /// Build a phylogeny viewer with given parameters
  ///
  /// \param parent The dialog's parent (can be null)
  /// \param ptree The tree of which this object is a view
  /// \param direction The layout direction
  /// \param config The initial configuration
  PhylogenyViewer(QWidget *parent, PTree &ptree,
                  Direction direction = Direction::TopToBottom,
                  Config config = defaultConfig())
    : PhylogenyViewer_base(parent, config), _ptree(ptree), callbacks(this) {
    _ptree.setCallbacks(&callbacks);

    uint step = _ptree.step();
    constructorDelegate(step, direction);

    build();
  }

  /// Helper function for getting a tree building cache
  auto cache (void) {
    return PTreeBuildingCache {
      this, _config, _ptree.step(), _items
    };
  }

  /// Request full parsing of the associated PTree for a complete graph generation
  void build (void) {
    auto c = cache();
    Builder::fillScene(_ptree, c);
    Builder::updateLayout(_items);

    updatePens();
    makeFit(_config.autofit);
  }

  /// Requests the base class to render the current view
  void render (void) {
    PhylogenyViewer_base::render(_ptree.step());
  }

  /// Process a new species event (add a new node to the graph)
  void newSpecies (SID pid, SID sid) {
    auto c = cache();
    Node *parent = (pid != SID::INVALID) ? _items.nodes[pid] : nullptr;
    const auto &pn = *_ptree.nodeAt(sid);
    Builder::addSpecies(parent, pn, c);
    Builder::updateLayout(_items);
    _items.border->setEmpty(false);
    _view->update();
  }

protected:
  void hoverEvent (SID sid, bool entered) override {
    if (entered)
          _items.contributors->show(sid, _items,
                                    _ptree.nodeAt(sid)->contributors);
    else  _items.contributors->hide();

    emit onSpeciesHoverEvent(sid, entered);
  }

  void doubleClickEvent (const Node &gn, QGraphicsSceneMouseEvent *e) override {
    const typename PTree::Node &n = *_ptree.nodeAt(gn.id);

    QStringList data;
    std::vector<GENOME> genomes;

    data.append(gn.computeTooltip());
    for (const auto &ep: n.enveloppe) {
      data.append(dumpEnveloppePoint(ep));
      genomes.push_back(ep.genome);
    }

    std::ostringstream oss;
    GENOME::aggregate(oss, genomes, config::PTree::speciesDetailVerbosity());

    speciesDetailPopup(gn.id, data, QString::fromStdString(oss.str()), e);
  }

private:
  /// The PTree associated to this view
  const PTree &_ptree;

  /// The callbacks used by #_ptree
  PTCallbacks callbacks;

  void updateLayout (void) override {
    Builder::updateLayout(_items);
    update();
  }

  /// \returns a description of the data contained by this enveloppe point
  QString dumpEnveloppePoint (const typename PTree::Node::EnvPoint &ep) {
    QString s;
    s += "Genome: ";
    s += QString::fromStdString(nlohmann::json(ep.genome).dump(2));
    s += "\nUser data: ";
    s += QString::fromStdString(nlohmann::json(*ep.userData).dump(2));
    s += "\n";
    return s;
  }
};



} // end namespace gui

/// \brief Specialization of the callbacks for a PTree with parameter \p GENOME.
///
/// Forwards events to the viewer.
template <typename GENOME, typename UDATA>
struct phylogeny::Callbacks_t<phylogeny::PhylogeneticTree<GENOME, UDATA>> {
  /// The PTree type at the source of these callbacks
  using PT = phylogeny::PhylogeneticTree<GENOME, UDATA>;

  /// The PViewer type receiving these forwarded events
  using PV = gui::PhylogenyViewer<GENOME, UDATA>;

  /// Helper alias to a genomic identificator
  using GID = phylogeny::GID;

  /// Helper alias to a species identificator
  using SID = phylogeny::SID;

  /// Helper alias to a collection of still-alive species
  using LivingSet = phylogeny::LivingSet;

  /// Creates a callback object associated with a specific viewer
  Callbacks_t (PV *v) : viewer(v) {}

  /// Notify both the viewer and the outside world that the associated tree has
  /// been stepped
  ///
  /// \copydetails phylogeny::Callbacks_t::onStepped
  void onStepped (uint step, const LivingSet &living) {
    viewer->treeStepped(step, living);
    emit viewer->onTreeStepped(step, living);
  }

  /// Notify both the viewer and the outside world that a new species has been
  /// added to the associated tree
  ///
  /// \copydetails phylogeny::Callbacks_t::onNewSpecies
  void onNewSpecies (SID pid, SID sid) {
    viewer->newSpecies(pid, sid);
    emit viewer->onNewSpecies(pid, sid);
  }

  /// Notify both the viewer and the outside world that a species of the associated
  /// tree has changed
  ///
  /// \copydetails phylogeny::Callbacks_t::onGenomeEntersEnveloppe
  void onGenomeEntersEnveloppe (SID sid, GID gid) {
    viewer->genomeEntersEnveloppe(sid, gid);
    emit viewer->onGenomeEntersEnveloppe(sid, gid);
  }

  /// Notify both the viewer and the outside world that a species of the associated
  /// tree has changed
  ///
  /// \copydetails phylogeny::Callbacks_t::onGenomeLeavesEnveloppe
  void onGenomeLeavesEnveloppe (SID sid, GID gid) {
    viewer->genomeLeavesEnveloppe(sid, gid);
    emit viewer->onGenomeLeavesEnveloppe(sid, gid);
  }

  /// Notify both the viewer and the outside world that a species of the
  /// associated tree has changed its rooting point
  ///
  /// \copydetails phylogeny::Callbacks_t::onMajorContributorChanged
  void onMajorContributorChanged (SID sid, SID oldMC, SID newMC) {
    viewer->majorContributorChanged(sid, oldMC, newMC);
    emit viewer->onMajorContributorChanged(sid, oldMC, newMC);
  }

private:
  /// The associated viewer
  PV *viewer;
};

#endif // PHYLOGENYVIEWER_H
