// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

// SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file AUTHORS.md
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception OR LGPL-3.0-or-later
#ifndef DUNE_FUNCTIONS_COMMON_SQUEEZETENSOR_HH
#define DUNE_FUNCTIONS_COMMON_SQUEEZETENSOR_HH

namespace Dune::Functions
{
  namespace Impl
  {
    /** \brief Remove empty axes from tensors
     *
     * This default overload returns the unchanged object.
     */
    template<class Object>
    constexpr Object& squeezeTensor(Object& o){
      return o;
    }

    /** \brief Remove empty axes from tensors
     *
     * Default overload for `const` objects. It returns the unchanged object.
     */
    template<class Object>
    constexpr Object const& squeezeTensor(Object const& o){
      return o;
    }

    /** \brief Const overload of \ref squeeze for 1-dimensional vectors
     *
     * Overload for vectors; returns a scalar
     */
    template<class K>
    constexpr K const& squeezeTensor(Dune::FieldVector<K,1> const& v){
      return v[0];
    }

    /** \brief Mutable overload of \ref squeeze for 1-dimensional vectors
     */
    template<class K>
    constexpr K& squeezeTensor(Dune::FieldVector<K,1>& v){
      return v[0];
    }

    /** \brief Const overload of \ref squeeze for 1 x N matrices
     */
    template<class K, int N>
    constexpr Dune::FieldVector<K, N> const& squeezeTensor(Dune::FieldMatrix<K,1,N> const& m){
      return m[0];
    }

    /** \brief Mutable overload of \ref squeeze for 1 x N-dimensional matrices
     */
    template<class K, int N>
    constexpr Dune::FieldVector<K, N>& squeezeTensor(Dune::FieldMatrix<K,1,N>& m){
      return m[0];
    }

    /** \brief Const overload of \ref squeeze for 1,N,M-dimensional tensors
     */
    template<class K, int N, int M>
    constexpr Dune::FieldMatrix<K, N, M> const& squeezeTensor(std::array<Dune::FieldMatrix<K,N,M>, 1> const& m){
      return m[0];
    }

    /** \brief Mutable overload of \ref squeeze for 1,N,M-dimensional tensors
     */
    template<class K, int N, int M>
    constexpr Dune::FieldMatrix<K, N, M>& squeezeTensor(std::array<Dune::FieldMatrix<K,N,M>, 1>& m){
      return m[0];
    }

  }  // namespace Impl

}  // namespace Dune::Functions

#endif  // DUNE_FUNCTIONS_COMMON_SQUEEZETENSOR_HH
