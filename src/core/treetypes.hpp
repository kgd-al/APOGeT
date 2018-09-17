#ifndef _PTREE_TYPES_HPP_
#define _PTREE_TYPES_HPP_

#include <type_traits>
#include <set>

namespace phylogeny {

/// Helper container for the phylogeny-related types
struct TreeTypes {

  /// Alias for the species identificator
  enum class SID : uint {
    INVALID = uint(-1)  ///< Value indicating an unspecified species
  };

  /// Auto convert outstream operator
  friend std::ostream& operator<< (std::ostream &os, SID sid) {
    return os << std::underlying_type<SID>::type(sid);
  }

  /// Collections of still-alive species identificators
  using LivingSet = std::set<SID>;
};

} // end of namespace phylogeny

#endif // _PTREE_TYPES_HPP_
