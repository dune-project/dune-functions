// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

// SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file AUTHORS.md
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception OR LGPL-3.0-or-later

#ifndef DUNE_FUNCTIONS_COMMON_DENSEVECTORVIEW_HH
#define DUNE_FUNCTIONS_COMMON_DENSEVECTORVIEW_HH

#include <type_traits>

#include <dune/common/densevector.hh>


namespace Dune::Functions::Impl {

  /**
   * \brief Wrapper providing the DenseVector interface for a random access container
   */
  template<class R>
  class DenseVectorView
    : public DenseVector<DenseVectorView<R>>
  {
    R& data_;
    using Base = DenseVector<DenseVectorView<R>>;

    using mutable_reference = typename std::decay_t<R>::reference;

  public:

    //! The type used for array indices and sizes
    using size_type = typename std::decay_t<R>::size_type;

    //! The type of values
    using value_type = typename std::decay_t<R>::value_type;

    //! The type used for const references to the vector entry
    using const_reference = typename std::decay_t<R>::const_reference;

    //! The type used for references to the vector entry
    using reference = std::conditional_t<std::is_const_v<R>,
      const_reference,
      mutable_reference>;

    //! Construct from a pointer to a scalar
    DenseVectorView (R& data)
      : data_(data)
    {}

    //! Copy constructor
    DenseVectorView (const DenseVectorView &other) :
      Base(),
      data_(other.data_)
    {}

    //! Move constructor
    DenseVectorView (DenseVectorView &&other) :
      Base(),
      data_( other.data_ )
    {}

    //! Copy assignment operator
    DenseVectorView& operator= (const DenseVectorView& other)
    {
      data_ = other.data_;
      return *this;
    }

    //! Copy assignment operator
    template<class RR>
    DenseVectorView& operator= (const DenseVectorView<RR>& other)
    {
      data_ = other.data_;
      return *this;
    }

    //! Container size
    size_type size () const
    {
      return data_.size();
    }

    //! Random access operator
    reference operator[] (size_type i)
    {
      return data_[i];
    }

    //! Random access operator
    const_reference operator[] (size_type i) const
    {
      return data_[i];
    }

  }; // class DenseVectorView

}

namespace Dune {

  template< class R>
  struct DenseMatVecTraits< Dune::Functions::Impl::DenseVectorView<R> >
  {
    using derived_type = Dune::Functions::Impl::DenseVectorView<R>;
    using value_type = typename R::value_type;
    using size_type = typename R::size_type;
  };

  template< class R >
  struct FieldTraits< Dune::Functions::Impl::DenseVectorView<R> >
    : public FieldTraits<std::remove_const_t<typename Dune::Functions::Impl::DenseVectorView<R>::value_type>>
  {};

}


#endif // DUNE_FUNCTIONS_COMMON_DENSEVECTORVIEW_HH
