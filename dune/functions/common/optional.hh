// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_FUNCTIONS_COMMON_OPTIONAL_HH
#define DUNE_FUNCTIONS_COMMON_OPTIONAL_HH

#include <utility>

#include <dune/common/typeutilities.hh>

namespace Dune {
namespace Functions {


/**
 * \brief A wrapper that can either contain an object of T or be empty
 *
 * \ingroup Utility
 *
 * \tparam T Type of wrapped objects
 */
template<class T>
class Optional
{
public:

  Optional() :
    p_(nullptr)
  {}

  template<class TT, disableCopyMove<Optional, TT> = 0>
  Optional(TT&& t) :
    p_(nullptr)
  {
    emplace(std::forward<TT>(t));
  }

  Optional(Optional&& other)
  {
    if (other)
      p_ = new (buffer_) T(std::move(other.value()));
    else
      p_ = nullptr;
  }

  Optional(const Optional& other)
  {
    if (other)
      p_ = new (buffer_) T(other.value());
    else
      p_ = nullptr;
  }

  ~Optional()
  {
    if (operator bool())
      p_->~T();
  }

  template<class TT, disableCopyMove<Optional, TT> = 0 >
  Optional& operator=(TT&& t)
  {
    if (operator bool())
      *p_ = std::forward<T>(t);
    else
      p_ = new (buffer_) T(std::forward<T>(t));
    return *this;
  }

  Optional& operator=(const Optional& other)
  {
    if (other)
      *this = other.value();
    else if (operator bool())
    {
      p_->~T();
      p_ = nullptr;
    }
    return *this;
  }

  Optional& operator=(Optional&& other)
  {
    if (other)
      *this = std::move(other.value());
    else if (operator bool())
    {
      p_->~T();
      p_ = nullptr;
    }
    return *this;
  }

  explicit operator bool() const
  {
    return p_;
  }

  const T& value() const
  {
    return *p_;
  }

  T& value()
  {
    return *p_;
  }

  template< class... Args >
  void emplace(Args&&... args)
  {
    if (operator bool())
      p_->~T();
    p_ = new (buffer_) T(std::forward<Args>(args)...);
  }

  void release()
  {
    if (operator bool())
    {
      p_->~T();
      p_ = nullptr;
    }
  }

private:

  alignas(T) char buffer_[sizeof(T)];
  T* p_;
};


} // namespace Functions
} // namespace Dune

#endif // DUNE_FUNCTIONS_COMMON_OPTIONAL_HH
