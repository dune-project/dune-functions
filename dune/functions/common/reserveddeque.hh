// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_FUNCTIONS_COMMON_RESERVEDDEQUE_HH
#define DUNE_FUNCTIONS_COMMON_RESERVEDDEQUE_HH

/** \file
 * \brief An stl-compliant double-ended queue which stores everything on the stack
 */

#include <algorithm>
#include <iostream>
#include <cstddef>
#include <initializer_list>

#include <dune/common/genericiterator.hh>

#ifdef CHECK_RESERVEDDEQUE
#define CHECKSIZE(X) assert(X)
#else
#define CHECKSIZE(X) {}
#endif

namespace Dune {
namespace Functions {

  /**
     \brief A Vector class with statically reserved memory.

     \ingroup Utility

     ReservedDeque is something between Dune::array and std::deque.
     You have a double ended queue which can be extended and shrunk using methods like
     push_back and pop_back at the end, or via push_front and pop_front,
     but reserved memory is predefined.

     This implies that the vector can not grow bigger than the predefined
     maximum size.

     \tparam T The data type ReservedVector stores.
     \tparam n The maximum number of objects the ReservedVector can store.

   */
  template<class T, int n>
  class ReservedDeque
  {
  public:

    /** @{ Typedefs */

    //! The type of object, T, stored in the vector.
    typedef T value_type;
    //! Pointer to T.
    typedef T* pointer;
    //! Reference to T
    typedef T& reference;
    //! Const reference to T
    typedef const T& const_reference;
    //! An unsigned integral type.
    typedef size_t size_type;
    //! A signed integral type.
    typedef std::ptrdiff_t difference_type;
    //! Iterator used to iterate through a vector.
    typedef Dune::GenericIterator<ReservedDeque, value_type> iterator;
    //! Const iterator used to iterate through a vector.
    typedef Dune::GenericIterator<const ReservedDeque, const value_type> const_iterator;

    /** @} */

    /** @{ Constructors */

    //! Constructor
    ReservedDeque() :
      size_(0),
      first_(0)
    {}

    ReservedDeque(std::initializer_list<T> const &l)
    {
      assert(l.size() <= n);// Actually, this is not needed any more!
      size_ = l.size();
      std::copy_n(l.begin(), size_, data_);
    }

    /** @} */

    /** @{ Data access operations */

    //! Erases all elements.
    void clear()
    {
      first_ = 0;
      size_ = 0;
    }

    //! Specifies a new size for the vector.
    void resize(size_t s)
    {
      CHECKSIZE(s<=n);
      size_ = s;
    }

    //! Appends an element to the end of a vector, up to the maximum size n, O(1) time.
    void push_back(const T& t)
    {
      CHECKSIZE(size_<n);
      data_[size_++ % n] = t;
    }

    //! Prepends an element to the begin of a deque if size<capacity, O(1) time.
    void push_front(const T& t)
    {
      CHECKSIZE(size_<n);
      ++size_;
      first_ = (first_+(n-1)) % n;
      data_[first_] = t;
    }

    //! Erases the last element of the vector, O(1) time.
    void pop_back()
    {
      if (! empty())
        size_--;
    }

    //! Erases the first element of the vector, O(1) time.
    void pop_front()
    {
      if (! empty())
      {
        size_--;
        first_ = (++first_) % n;
      }
    }

    //! Returns a iterator pointing to the beginning of the vector.
    iterator begin(){
      return iterator(*this, 0);
    }

    //! Returns a const_iterator pointing to the beginning of the vector.
    const_iterator begin() const {
      return const_iterator(*this, 0);
    }

    //! Returns an iterator pointing to the end of the vector.
    iterator end(){
      return iterator(*this, size_);
    }

    //! Returns a const_iterator pointing to the end of the vector.
    const_iterator end() const {
      return const_iterator(*this, size_);
    }

    //! Returns reference to the i'th element.
    reference operator[] (size_type i)
    {
      CHECKSIZE(size_>i);
      return data_[(first_ + i) % n];
    }

    //! Returns a const reference to the i'th element.
    const_reference operator[] (size_type i) const
    {
      CHECKSIZE(size_>i);
      return data_[(first_ + i) % n];
    }

    //! Returns reference to first element of vector.
    reference front()
    {
      CHECKSIZE(size_>0);
      return data_[first_];
    }

    //! Returns const reference to first element of vector.
    const_reference front() const
    {
      CHECKSIZE(size_>0);
      return data_[first_];
    }

    //! Returns reference to last element of vector.
    reference back()
    {
      CHECKSIZE(size_>0);
      return data_[(first_ + size_-1) % n];
    }

    //! Returns const reference to last element of vector.
    const_reference back() const
    {
      CHECKSIZE(size_>0);
      return data_[(first_ + size_-1) % n];
    }

    /** @} */

    /** @{ Informative Methods */

    //! Returns number of elements in the vector.
    size_type size () const
    {
      return size_;
    }

    //! Returns true if vector has no elements.
    bool empty() const
    {
      return size_==0;
    }

    //! Returns current capacity (allocated memory) of the vector.
    static DUNE_CONSTEXPR size_type capacity()
    {
      return n;
    }

    //! Returns the maximum length of the vector.
    static DUNE_CONSTEXPR size_type max_size()
    {
      return n;
    }

    /** @} */

    //! Send ReservedVector to an output stream
    friend std::ostream& operator<< (std::ostream& s, const ReservedDeque& v)
    {
      for (size_t i=0; i<v.size(); i++)
        s << v[i] << "  ";
      return s;
    }

  private:
    T data_[n];
    size_type first_;
    size_type size_;
  };

} // end namespace Functions
} // end namespace Dune

#undef CHECKSIZE

#endif // DUNE_FUNCTIONS_COMMON_RESERVEDDEQUE_HH
