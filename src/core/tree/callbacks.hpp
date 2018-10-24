#ifndef CALLBACKS_HPP
#define CALLBACKS_HPP

/*!
 * \file callbacks.hpp
 *
 * Contains the definition for the callbacks sent by the phylogenic tree
 */

namespace phylogeny {

/// \cond internal

/// Contains a set of functions called by the phylogenic tree when each of the
/// corresponding events occur
template <typename PT>
struct Callbacks_t {
  /// Alias to the genetic ID
  using GID = typename PT::GID;

  /// Alias to the species ID
  using SID = typename PT::SID;

  /// Alias to the collection of still-alive species id
  using LivingSet = typename PT::LivingSet;

  /// \brief Called when the PTree has been stepped.
  ///
  /// Provides the current step and the set of still-alive species
  void onStepped (uint step, const LivingSet &living);

  /// \brief Called to notify of a newly created species.
  ///
  /// Provides the identificators of both parent (if any) and new species
  void onNewSpecies (SID pid, SID sid);

  /// \brief Called when a genome is added to an enveloppe.
  ///
  /// Provides the id of the species whose enveloppe just changed and the id
  /// of the newly inserted genome
  void onGenomeEntersEnveloppe (SID sid, GID gid);

  /// \brief Called when a genome is removed from an enveloppe.
  ///
  /// Provides the id of the species whose enveloppe just changed and the id
  /// of the newly removed genome
  void onGenomeLeavesEnveloppe (SID sid, GID gid);

  /// \brief Called when a species' major contributor has changed
  ///
  /// Provides the id of the species whose major contributor just changed
  /// as well as the id of the previous and new MC
  void onMajorContributorChanged (SID sid, SID oldMC, SID newMC);
};

/// \endcond

} // end of namespace phylogeny

#endif // CALLBACKS_HPP
