// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_FUNCTIONS_COMMON_TUPLEVECTOR_HH
#define DUNE_FUNCTIONS_COMMON_TUPLEVECTOR_HH

#include <tuple>

#include <dune/common/indices.hh>

namespace Dune
{
namespace Functions
{

  /**
   * \brief A class augmenting std::tuple by element access via operator[]
   *
   * \ingroup Utility
   */
  template<class... T>
  class TupleVector : public std::tuple<T...>
  {
    using Base = std::tuple<T...>;

  public:

    /** \brief Construct from a set of arguments */
    template<class... TT>
    constexpr TupleVector(TT&&... tt) :
      Base(std::forward<TT>(tt)...)
    {}

    /** \brief Default constructor */
    constexpr TupleVector()
    {}

    /** \brief Const access to the tuple elements */
    template<std::size_t i>
    auto operator[](const Dune::index_constant<i>&) const
      ->decltype(std::get<i>(*this))
    {
      return std::get<i>(*this);
    }

    /** \brief Non-const access to the tuple elements */
    template<std::size_t i>
    auto operator[](const Dune::index_constant<i>&)
      ->decltype(std::get<i>(*this))
    {
      return std::get<i>(*this);
    }

    /** \brief Number of elements of the tuple */
    static constexpr std::size_t size()
    {
      return std::tuple_size<Base>::value;
    }

  };

}  // namespace Functions

}  // namespace Dune

#endif
