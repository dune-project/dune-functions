// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_FUNCTIONS_COMMON_OVERFLOWARRAY_HH
#define DUNE_FUNCTIONS_COMMON_OVERFLOWARRAY_HH

#include <algorithm>
#include <iostream>
#include <cstddef>
#include <array>
#include <initializer_list>

#include <dune/common/genericiterator.hh>



namespace Dune::Functions {


/**
 * \brief A dynamically sized array-like class with overflow
 *
 * \tparam BA Type of base array
 * \tparam maxSize Maximal size of OverflowArray
 *
 * This class publicly inherits from a statically sized array-like
 * base class BA and extends it by an overflow such that
 * a total capacity of maxSize is available. Within this bound
 * the size is managed dynamically.
 *
 * Potential usecase: If you want to construct a statically
 * sized array but need dynamic resizing while building it,
 * you can use an OverflowArray<std::array<T,finalSize>, ...>
 * and cast the result to the base class type.
 *
 * Similar to Dune::ReservedVector this uses a std::array
 * internally with the following implications: Entries must be
 * default-constructible. The whole capacity will always be filled
 * with entries, even if size<capacity. Entries are only destructed
 * when the OverflowArray is destructed - not when shrinking
 * or clearing it.
 */
template<class BA, std::size_t maxSize = std::tuple_size_v<BA>>
class OverflowArray :
  public BA
{
  static constexpr std::size_t baseSize = std::tuple_size_v<BA>;

public:
  using BaseArray = BA;

  using value_type = typename BaseArray::value_type;
  using reference = value_type&;
  using const_reference = const value_type&;
  using pointer = value_type*;
  using difference_type = std::ptrdiff_t;
  using size_type = std::size_t;
  using iterator = Dune::GenericIterator<OverflowArray, value_type>;
  using const_iterator = Dune::GenericIterator<const OverflowArray, const value_type>;

private:
  using OverflowBuffer = std::array<value_type, maxSize-baseSize>;

public:

  OverflowArray() = default;

  OverflowArray(const std::initializer_list<value_type>& l) {
    assert(l.size() <= capacity());
    size_ = l.size();
    std::copy_n(l.begin(), size_, begin());
  }

  bool operator == (const OverflowArray& other) const {
    if (size() != other.size())
      return false;
    for (size_type i=0; i<size(); ++i)
      if ((*this)[i] != other[i])
        return false;
    return true;
  }

  //! Erases all elements.
  void clear() {
    size_ = 0;
  }

  /**
   * \brief Specifies a new size for the OverflowArray.
   *
   * The new size must not exceed max_size().
   * This is an O(1) operation.
   */
  void resize(size_type n) {
    assert(n <= capacity());
    size_ = n;
  }

  /**
   * \brief Appends an element to the end of the OverflowArray,
   *
   * The new size must not exceed max_size().
   * This is an O(1) operation.
   */
  void push_back(const value_type& t) {
    assert(size() < capacity());
    (*this)[size_++] = t;
  }

  //! Erases the last element of the OverflowArray, O(1) time.
  void pop_back() {
    assert(size() > 0);
    if (! empty())
      size_--;
  }

  /**
   * \brief Inserts an element to the begin of the OverflowArray,
   *
   * The new size must not exceed max_size().
   * This is an O(size()) operation.
   */
  void push_front(const value_type& t) {
    assert(size() < capacity());
    for (size_type i=0; i<size(); i++)
      (*this)[i+1] = (*this)[i];
    (*this)[0] = t;
  }

  //! Returns a iterator pointing to the beginning of the OverflowArray.
  iterator begin() {
    return iterator(*this, 0);
  }

  //! Returns a const_iterator pointing to the beginning of the OverflowArray.
  const_iterator begin() const {
    return const_iterator(*this, 0);
  }

  //! Returns an iterator pointing to the end of the OverflowArray.
  iterator end() {
    return iterator(*this, size());
  }

  //! Returns a const_iterator pointing to the end of the OverflowArray.
  const_iterator end() const {
    return const_iterator(*this, size());
  }

  //! Returns reference to the i'th element.
  reference operator[] (size_type i) {
    assert(i < size());
    // If there's no padding between the base class and the overflow_ member,
    // the compiler should be able to optimize this to
    // return *(&BaseArray::operator[](0) + i);
    if (i<baseSize)
      return BaseArray::operator[](i);
    return overflow_[i-baseSize];
  }

  //! Returns a const reference to the i'th element.
  const_reference operator[] (size_type i) const {
    assert(i < size());
    // If there's no padding between the base class and the overflow_ member,
    // the compiler should be able to optimize this to
    // return *(&BaseArray::operator[](0) + i);
    if (i<baseSize)
      return BaseArray::operator[](i);
    return overflow_[i-baseSize];
  }

  //! Returns reference to first element of OverflowArray.
  reference front() {
    assert(size() > 0);
    return (*this)[0];
  }

  //! Returns const reference to first element of OverflowArray.
  const_reference front() const {
    assert(size() > 0);
    return (*this)[0];
  }

  //! Returns reference to last element of OverflowArray.
  reference back() {
    assert(size() > 0);
    return (*this)[size()-1];
  }

  //! Returns const reference to last element of OverflowArray.
  const_reference back() const {
    assert(size() > 0);
    return (*this)[size()-1];
  }

  //! Returns number of elements in the OverflowArray.
  size_type size () const {
    return size_;
  }

  //! Returns true if OverflowArray has no elements.
  bool empty() const {
    return size() == 0;
  }

  //! Returns the capacity of the  OverflowArray.
  static constexpr size_type capacity() {
    return maxSize;
  }

  //! Returns the maximum length of the OverflowArray.
  static constexpr size_type max_size() {
    return maxSize;
  }

  //! Compute hash value
  inline friend std::size_t hash_value(const OverflowArray& v) noexcept {
    return hash_range(v.begin(), v.end());
  }

  //! Write container to an output stream
  friend std::ostream& operator<< (std::ostream& s, const OverflowArray& c) {
    for (const auto& ci : c)
      s << ci << "  ";
    return s;
  }

private:
  OverflowBuffer overflow_;
  size_type size_ = 0;
};



} // namespace Dune::Functions



#endif // DUNE_FUNCTIONS_COMMON_OVERFLOWARRAY_HH
