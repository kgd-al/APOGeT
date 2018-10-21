#ifndef _PTREE_TYPES_HPP_
#define _PTREE_TYPES_HPP_

#include <type_traits>
#include <vector>
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

  /// std::vector extension for managing collections indexed by species
  /// identificators
  template <typename T>
  class sidvector {
    /// Helper alias to the internal vector type
    using V = std::vector<T>;

    /// Helper alias to the integral equivalent of SID
    using SID_t = std::underlying_type<SID>::type;

    V vec;  ///< Internal vector

  public:
    /// Build from nothing (empty vector)
    sidvector(void) {}

    /// Build from initializer list
    sidvector(std::initializer_list<T> &&list) : vec(list) {}

    /// Build from static C array
    template <size_t N>
    sidvector(const T (&a) [N]) : vec(a, std::end(a)) {}

    /// Access an element
    auto& operator[] (SID sid) {
      return vec[SID_t(sid)];
    }

    /// Access an immutable element
    const auto& operator[] (SID sid) const {
      return vec[SID_t(sid)];
    }

    /// Access an element (with bounds checking)
    auto& at (SID sid) {
      return vec.at(SID_t(sid));
    }

    /// Access an immutable element (with bounds checking)
    const auto& at (SID sid) const {
      return vec.at(SID_t(sid));
    }

    /// \overload
    void push_back (const T& val) {
      vec.push_back(val);
    }

    /// \overload
    void push_back (T&& val) {
      vec.push_back(val);
    }

    /// \overload
    auto resize (typename V::size_type n) {
      vec.resize(n);
    }
  };
};

} // end of namespace phylogeny

#endif // _PTREE_TYPES_HPP_
