#ifndef PHYLOGENYVIEWER_H
#define PHYLOGENYVIEWER_H

#include <QDialog>
#include <QGraphicsView>
#include <QBoxLayout>

#include "ptgraphbuilder.h"

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

  /// \copydoc phylogeny::TreeTypes::SID
  using SID = phylogeny::TreeTypes::SID;

  /// \copydoc phylogeny::TreeTypes::LivingSet
  using LivingSet = phylogeny::TreeTypes::LivingSet;

  /// Configuration data controlling what to draw and how
  using Config = gui::ViewerConfig;

  /// The direction in which to layout components
  using Direction = QBoxLayout::Direction;

public:
  /// Create a phylogeny viewer with given \p parent and using \p config as
  /// its initial configuration
  PhylogenyViewer_base (QWidget *parent, Config config)
    : QDialog(parent), _config(config) {}

  /// Render the tree to file
  void render (uint step);

  /// Intercepted to allow for tree autofitting
  void resizeEvent(QResizeEvent *event);

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

protected slots:
  // ===========================================================================
  // Callbacks from PTree

  /// Process a step event (new timestamp/living species)
  void treeStepped (uint step, const LivingSet &living);

  /// Move to template version (need to provide the node)
//  void newSpecies (uint pid, uint sid);

  /// Process a enveloppe change event
  void genomeEntersEnveloppe (SID sid, GID gid);

  /// \copydoc genomeEntersEnveloppe
  void genomeLeavesEnveloppe (SID sid, GID gid);

  // ===========================================================================
  // Config update

  /// Called when the corresponding slider's value has been changed
  void updateMinSurvival (uint v);

  /// \copydoc updateMinSurvival
  void updateMinEnveloppe (int v);

  /// Called when the corresponding checkbox has been changed
  void toggleShowNames (void);

public slots:
  /// Requests the scale of the view to be adapted to the size of the scene
  void makeFit (bool autofit);

  /// Prints the current scene to the image file \p filename
  void renderTo (QString filename = "");

  /// Prints the current scene into a pixmap
  QPixmap renderToPixmap (void) const {
    return renderToPixmap(_items.scene->sceneRect().size().toSize());
  }

  /// Prints the current scene into a pixmap of size \p requestedSize
  QPixmap renderToPixmap (const QSize &requestedSize) const;

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
template <typename GENOME>
class PhylogenyViewer : public PhylogenyViewer_base {
public:

  /// Helper alias to the PTree visualized
  using PTree = phylogeny::PhylogenicTree<GENOME>;

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
      _config, _ptree.step(), _items,
      [this] (auto sid, auto entered) {
        emit onSpeciesHoverEvent(sid, entered);
      }
    };
  }

  /// Request full parsing of the associated PTree for a complete graph generation
  void build (void) {
    auto c = cache();
    Builder::fillScene(_ptree, c);
    Builder::updateLayout(_items);
    makeFit(_config.autofit);
  }

  /// Requests the base class to render the current view
  void render (void) {
    PhylogenyViewer_base::render(_ptree.step());
  }

public:
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

private:
  /// The PTree associated to this view
  const PTree &_ptree;

  /// The callbacks used by #_ptree
  PTCallbacks callbacks;

  void updateLayout (void) override {
    Builder::updateLayout(_items);
  }
};

} // end namespace gui

/// \brief Specialization of the callbacks for a PTree with parameter \p GENOME.
///
/// Forwards events to the viewer.
template <typename GENOME>
struct phylogeny::Callbacks_t<phylogeny::PhylogenicTree<GENOME>> {
  /// The PTree type at the source of these callbacks
  using PT = phylogeny::PhylogenicTree<GENOME>;

  /// The PViewer type receiving these forwarded events
  using PV = gui::PhylogenyViewer<GENOME>;

  /// Helper alias to a genomic identificator
  using GID = typename PT::GID;

  /// Helper alias to a species identificator
  using SID = typename PT::SID;

  /// Helper alias to a collection of still-alive species
  using LivingSet = phylogeny::TreeTypes::LivingSet;

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

private:
  /// The associated viewer
  PV *viewer;
};

#endif // PHYLOGENYVIEWER_H
