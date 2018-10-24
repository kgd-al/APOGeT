#ifndef ENUM_VECTOR_HPP
#define ENUM_VECTOR_HPP

/*!
 * \file enumvector.hpp
 *
 * Contains the definition for a specialized vector
 */

#include <vector>

namespace phylogeny {

/// std::vector extension for managing collections indexed by an enumeration
template <typename ENUM, typename T>
class enumvector {
  /// Helper alias to the internal vector type
  using V = std::vector<T>;

  /// Helper alias to the integral equivalent of SID
  using ENUM_t = typename std::underlying_type<ENUM>::type;

  V vec;  ///< Internal vector

public:
  /// Build from nothing (empty vector)
  enumvector(void) {}

  /// Build from initializer list
  enumvector(std::initializer_list<T> &&list) : vec(list) {}

  /// Build from static C array
  template <size_t N>
  enumvector(const T (&a) [N]) : vec(a, std::end(a)) {}

  /// Access an element
  auto& operator[] (ENUM sid) {
    return vec[ENUM_t(sid)];
  }

  /// Access an immutable element
  const auto& operator[] (ENUM sid) const {
    return vec[ENUM_t(sid)];
  }

  /// Access an element (with bounds checking)
  auto& at (ENUM sid) {
    return vec.at(ENUM_t(sid));
  }

  /// Access an immutable element (with bounds checking)
  const auto& at (ENUM sid) const {
    return vec.at(ENUM_t(sid));
  }

  /// \return the size of the underlying buffer
  auto size (void) const {
    return vec.size();
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

  /// \overload
  auto begin (void) {
    return vec.begin();
  }

  /// \overload
  auto rbegin (void) {
    return vec.rbegin();
  }

  /// \overload
  auto end (void) {
    return vec.end();
  }

  /// \overload
  auto rend (void) {
    return vec.rend();
  }
};

} // end namespace phylogeny

#endif // ENUM_VECTOR_HPP
