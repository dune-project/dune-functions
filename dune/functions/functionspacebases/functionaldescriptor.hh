// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

// SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file AUTHORS.md
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception OR LGPL-3.0-or-later
#ifndef DUNE_FUNCTIONS_FUNCTIONSPACEBASES_FUNCTIONALDESCRIPTOR_HH
#define DUNE_FUNCTIONS_FUNCTIONSPACEBASES_FUNCTIONALDESCRIPTOR_HH

#include <array>


namespace Dune::Functions::Impl
{



/**
 * \brief A class describing a dual functional
 *
 * This class allows to describe which type of evaluations
 * is performed by the dual basis implemented by a
 * LocalFiniteElement::LocalInterpolation.
 * For a LocalInterpolation object lint
 * lint.functionalDescriptor(k) should return
 * a FunctionalDescriptor which describes the
 * evaluation type performed by the k-th
 * dual-basis function and thus specifies
 * the nature of the k-th DOF.
 *
 * \warning: This interface is highly experimental
 * and subject to change.
 */
template<std::size_t dim>
class FunctionalDescriptor
{
public:

  using Order = std::array<unsigned int, dim>;

  FunctionalDescriptor()
    : partialDerivativeOrder_{}
    , normalDerivativeOrder_(0)
  {}

  explicit FunctionalDescriptor(const Order& partialDerivativeOrder)
    : partialDerivativeOrder_{partialDerivativeOrder}
    , normalDerivativeOrder_(0)
  {}

  explicit FunctionalDescriptor(unsigned int normalDerivativeOrder)
    : partialDerivativeOrder_{}
    , normalDerivativeOrder_(normalDerivativeOrder)
  {}

  bool isNormalDerivative() const
  {
    return normalDerivativeOrder_>0;
  }

  bool isPartialDerivative() const
  {
    for(auto i: Dune::range(dim))
    {
      if (partialDerivativeOrder(i)>0)
        return true;
    }
    return false;
  }

  unsigned int normalDerivativeOrder() const
  {
    return normalDerivativeOrder_;
  }

  const Order& partialDerivativeOrder() const
  {
    return partialDerivativeOrder_;
  }

private:
  Order partialDerivativeOrder_;
  unsigned int normalDerivativeOrder_;
};



} // namespace Dune::Functions::Impl

#endif // DUNE_FUNCTIONS_FUNCTIONSPACEBASES_FUNCTIONALDESCRIPTOR_HH
