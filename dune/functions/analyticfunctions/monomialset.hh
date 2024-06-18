// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_FUNCTIONS_ANALYTICFUNCTIONS_MONOMIALSET_HH
#define DUNE_FUNCTIONS_ANALYTICFUNCTIONS_MONOMIALSET_HH

#include <array>

#include <dune/common/fvector.hh>
#include <dune/common/fmatrix.hh>
#include <dune/common/math.hh>
#include <dune/common/rangeutilities.hh>



namespace Dune::Functions {

  namespace Impl {

    // Computes `y = [1, x, x^2, x^3, ..., x^maxDegree]` with `x` a single coordinate
    template<int maxDegree, class DomainFieldType, class RangeType>
    void computePowers(const DomainFieldType& x, RangeType& y)
    {
      if constexpr(maxDegree >= 0)
      {
        y[0] = 1;
        for(auto k : Dune::range(maxDegree))
          y[k+1] = y[k]*x;
      }
    }

  } // namespace Impl


  /**
   * \brief Function, which evaluates all monomials up to degree \p maxDegree in
   * a given coordinate.
   *
   * \ingroup FunctionImplementations
   *
   * \tparam RangeFieldType scalar type.
   * \tparam dim Domain dimension.
   * \tparam maxDegree Maximal monomial degree.
   *
   * The `Range` of this (differentiable) function is a vector of monomial evaluations
   * `[1,x,y,z,xx,xy,yy,xz,yz,zz,...]` in the coordinate vector `[x,y,z]`. The
   * maximal degree of the provided monomial evaluations is given by the parameter
   * \p maxDegree. The number of coordinate components is given by the \p dimension.
   * From \p maxDegree and \p dimension the total number of monomials can be computed
   * as `binomial(maxDegree+dimension, dimension)` and is provided as static
   * constant `::size` in the class.

   * The function models the \ref Concept::DifferentiableFunction<Range(Domain)>
   * concept.
   *
   * The MonomialSet function is specialized for \p dim = 1, \p dim = 2, and
   * \p dim = 3 only.
   **/
  template<class RangeFieldType, int dimension, int maxDegree>
  struct MonomialSet
  {
    static constexpr int dim = dimension;
    static constexpr int size = Dune::binomial(maxDegree + dim, dim);
    static_assert(1 <= dim && dim <= 3, "MonomialSet is specialized for dimensions 1,2,3 only.");

    /**
     * \brief Return array of monomial exponents with shape `size x dim`
     *
     * The k-the entry of the returned array is the exponent
     * multi-index of the monomial corresponding the the k-th
     * component of the function. Note that the ordering is tensor based,
     * e.g., for the set of 3d monomials `1,x,y,z,xx,xy,yy,xz,yz,zz, ...` we get
     * the exponents:
     * `[ [0,0,0], [1,0,0], [0,1,0], [0,0,1], [2,0,0], [1,1,0], [0,2,0], [1,0,1],
     *    [0,1,1], [0,0,2], ...]`
     **/
    static constexpr std::array<std::array<std::size_t,dim>,size> exponents();

    /**
     * \brief Return array of monomial evaluations
     *
     * The k-the entry of the returned array is the value
     * of the monomial corresponding the the k-th
     * entry of the return value of exponents().
     *
     * \tparam DomainFieldType The scalar type of the domain.
     **/
    template<class DomainFieldType>
    constexpr Dune::FieldVector<RangeFieldType,size> operator()(const Dune::FieldVector<DomainFieldType,dim>& x) const;

    /**
     * \brief Set of all first order derivatives of monomials up to degree
     * \p maxDegree as vector of vector valued functions.
     **/
    struct Derivative
    {
      /**
       * \brief Return array of arrays of monomial derivatives
       *
       * The ith component of the k-the entry of the returned structure
       * is the derivative in direction i of the monomial corresponding
       * the the k-th entry of the return value of exponents().
       *
       * \tparam DomainFieldType The scalar type of the domain.
       **/
      template<class DomainFieldType>
      constexpr Dune::FieldMatrix<RangeFieldType,size, dim> operator()(const Dune::FieldVector<DomainFieldType,dim>& x) const;

      /**
       * \brief Set of all second order derivatives of monomials up to degree
       * \p maxDegree as vector of matrix valued functions.
       **/
      struct Hessian
      {
        /**
         * \brief Return array of Matrices of monomial derivatives
         *
         * The (i,j)th component of the k-the entry of the returned structure
         * is the derivative in direction (i,j) of the monomial corresponding
         * the the k-th entry of the return value of exponents().
         *
         * \tparam DomainFieldType The scalar type of the domain.
         **/
        template<class DomainFieldType>
        constexpr std::array<FieldMatrix<RangeFieldType,dim, dim>, size> operator()(const Dune::FieldVector<DomainFieldType,dim>& x) const;

      };

      /**
       * \brief Construct the Hessian object from a Derivative
       **/
      constexpr friend auto derivative(const Derivative & d)
      {
        return Hessian{};
      }

    };

    /**
     * \brief Construct the Derivative object from a MonomialSet
     **/
    constexpr friend auto derivative(const MonomialSet& m)
    {
      return Derivative{};
    }

  };


#ifndef DOXYGEN
  // Specialization for dim = 1
  template<class RangeFieldType, int maxDegree>
  struct MonomialSet<RangeFieldType, 1, maxDegree>
  {
    static constexpr int dim = 1;
    static constexpr int size = maxDegree+1;

    static constexpr auto exponents()
    {
      auto p = std::array<std::array<std::size_t,1>,size>{};
      for(auto k : Dune::range(size))
        p[k][0] = k;
      return p;
    }

    template<class DomainFieldType>
    constexpr auto operator()(const Dune::FieldVector<DomainFieldType,1>& x) const
    {
      auto y = Dune::FieldVector<RangeFieldType,size>{};
      Impl::computePowers<maxDegree>(x[0], y);
      return y;
    }

    struct Derivative
    {
      template<class DomainFieldType>
      constexpr auto operator()(const Dune::FieldVector<DomainFieldType,1>& x) const
      {
        auto xPowers = Dune::FieldVector<RangeFieldType,size>{};
        Impl::computePowers<maxDegree-1>(x[0], xPowers);

        auto y = Dune::FieldMatrix<RangeFieldType,size,1>{};
        for(auto degree : Dune::range(1, maxDegree+1))
          y[degree][0] = degree*xPowers[degree-1];
        return y;
      }

      struct Hessian
      {
        template<class DomainFieldType>
        constexpr auto operator()(const Dune::FieldVector<DomainFieldType,1>& x) const
        {
          auto xPowers = std::array<RangeFieldType,maxDegree+1>{};
          Impl::computePowers<maxDegree-2>(x[0],xPowers);

          auto y = std::array<Dune::FieldMatrix<RangeFieldType,1,1>,size>{};
          for(auto degree : Dune::range(maxDegree+1))
            if (degree-1 > 0)
              y[degree][0][0] = degree*(degree-1)*xPowers[degree-2];
          return y;
        }

      };

      constexpr friend auto derivative(const Derivative& d)
      {
        return Hessian{};
      }

    };

    constexpr friend auto derivative(const MonomialSet& m)
    {
      return Derivative{};
    }

  };


  // Specialization for dim = 2
  template<class RangeFieldType, int maxDegree>
  struct MonomialSet<RangeFieldType, 2, maxDegree>
  {
    static constexpr int dim = 2;
    static constexpr int size = (maxDegree+1)*(maxDegree+2)/2;

    static constexpr auto exponents()
    {
      auto p = std::array<std::array<std::size_t,2>,size>{};
      std::size_t index = 0;
      for(auto degree : Dune::range(maxDegree+1))
      {
        for(auto k : Dune::range(degree+1))
        {
          p[index][0] = degree-k;
          p[index][1] = k;
          ++index;
        }
      }
      return p;
    }

    template<class DomainFieldType>
    constexpr auto operator()(const Dune::FieldVector<DomainFieldType,2>& x) const
    {
      auto xPowers = std::array<Dune::FieldVector<RangeFieldType,size>,dim>{};
      for(auto j : Dune::range(dim))
        Impl::computePowers<maxDegree>(x[j], xPowers[j]);

      auto y = Dune::FieldVector<RangeFieldType,size>{};
      std::size_t index = 0;
      for(auto degree : Dune::range(maxDegree+1))
      {
        for(auto k : Dune::range(degree+1))
        {
          y[index] = xPowers[0][degree-k]*xPowers[1][k];
          ++index;
        }
      }
      return y;
    }

    struct Derivative
    {
      template<class DomainFieldType>
      constexpr auto operator()(const Dune::FieldVector<DomainFieldType,2>& x) const
      {
        auto xPowers = std::array<Dune::FieldVector<RangeFieldType,size>,dim>{};
        for(auto j : Dune::range(dim))
          Impl::computePowers<maxDegree-1>(x[j], xPowers[j]);

        auto y = Dune::FieldMatrix<RangeFieldType,size,2>{};
        std::size_t index = 0;
        for(auto degree : Dune::range(maxDegree+1))
        {
          for(auto k : Dune::range(degree+1))
          {
            if (degree-k > 0)
              y[index][0] = (degree-k)*xPowers[0][degree-k-1]*xPowers[1][k];
            if (k > 0)
              y[index][1] = k*xPowers[0][degree-k]*xPowers[1][k-1];
            ++index;
          }
        }
        return y;
      }

      struct Hessian
      {
        template<class DomainFieldType>
        constexpr auto operator()(const Dune::FieldVector<DomainFieldType,2>& x) const
        {
          auto xPowers = std::array<std::array<RangeFieldType,maxDegree+1>,dim>{};
          for(auto j : Dune::range(dim))
            Impl::computePowers<maxDegree-2>(x[j], xPowers[j]);

          auto y = std::array<Dune::FieldMatrix<RangeFieldType,2,2>,size>{};
          std::size_t index = 0;
          for(auto degree : Dune::range(maxDegree+1))
          {
            for(auto k : Dune::range(degree+1))
            {
              if (degree-k > 1)
                y[index][0][0] = (degree-k-1)*(degree-k)*xPowers[0][degree-k-2]*xPowers[1][k];
              if (k > 0 and degree-k > 0){
                auto mixed = k*(degree-k)*xPowers[0][degree-k-1]*xPowers[1][k-1];
                y[index][0][1]= mixed;
                y[index][1][0]= mixed;
              }
              if (k > 1)
                y[index][1][1] = k*(k-1)*xPowers[0][degree-k]*xPowers[1][k-2];

              ++index;
            }
          }
          return y;
        }

      };

      constexpr friend auto derivative(const Derivative & d)
      {
        return Hessian{};
      }

    };

    constexpr friend auto derivative(const MonomialSet& m)
    {
      return Derivative{};
    }

  };


  // Specialization for dim = 3
  template<class RangeFieldType, int maxDegree>
  struct MonomialSet<RangeFieldType, 3, maxDegree>
  {
    static constexpr int dim = 3;
    static constexpr int size = Dune::binomial(maxDegree + dim, dim);

    static constexpr auto exponents()
    {
      auto p = std::array<std::array<std::size_t,3>,size>{};
      std::size_t index = 0;
      for(auto degree : Dune::range(maxDegree+1))
      {
        for(auto k : Dune::range(degree+1))
        {
          for (auto l : Dune::range(degree-k+1))
          {
            p[index][0] = degree-k-l;
            p[index][1] = l;
            p[index][2] = k;
            ++index;
          }

        }
      }
      return p;
    }

    template<class DomainFieldType>
    constexpr auto operator()(const Dune::FieldVector<DomainFieldType,3>& x) const
    {
      auto xPowers = std::array<std::array<RangeFieldType,maxDegree+1>,dim>{};
      for(auto j : Dune::range(dim))
        Impl::computePowers<maxDegree>(x[j], xPowers[j]);

      auto y = Dune::FieldVector<RangeFieldType,size>{};
      std::size_t index = 0;
      for(auto degree : Dune::range(maxDegree+1))
      {
        for(auto k : Dune::range(degree+1))
        {
          for (auto l : Dune::range(degree-k+1))
          {
            y[index] = xPowers[0][degree-k-l]*xPowers[1][l]*xPowers[2][k];
            ++index;
          }
        }
      }
      return y;
    }

    struct Derivative
    {
      template<class DomainFieldType>
      constexpr auto operator()(const Dune::FieldVector<DomainFieldType,3>& x) const
      {
        auto xPowers = std::array<std::array<RangeFieldType,maxDegree+1>,dim>{};
        for(auto j : Dune::range(dim))
        {
          xPowers[j][0] = 1.0;
          for(auto k : Dune::range(maxDegree))
            xPowers[j][k+1] = xPowers[j][k]*x[j];
        }

        auto y = Dune::FieldMatrix<RangeFieldType,size,3>{};
        std::size_t index = 0;
        for(auto degree : Dune::range(maxDegree+1))
        {
          for(auto k : Dune::range(degree+1))
          {
            for (auto l : Dune::range(degree-k+1))
            {
              if (degree-k-l > 0)
                y[index][0] = (degree-k-l)*xPowers[0][degree-k-l-1]*xPowers[1][l]*xPowers[2][k];
              if (l > 0)
                y[index][1] = l*xPowers[0][degree-k-l]*xPowers[1][l-1]*xPowers[2][k];
              if (k > 0)
                y[index][2] = k*xPowers[0][degree-k-l]*xPowers[1][l]*xPowers[2][k-1];
              ++index;
            }
          }
        }
        return y;
      }

      struct Hessian
      {
        template<class DomainFieldType>
        constexpr auto operator()(const Dune::FieldVector<DomainFieldType,3>& x) const
        {
          auto xPowers = std::array<std::array<RangeFieldType,maxDegree+1>,dim>{};
          for(auto j : Dune::range(dim))
            Impl::computePowers<maxDegree>(x[j], xPowers[j]);

          auto y = std::array<Dune::FieldMatrix<RangeFieldType,3,3>,size>{};
          std::size_t index = 0;
          for(auto degree : Dune::range(maxDegree+1))
          {
            for(auto k : Dune::range(degree+1))
            {
              for (auto l : Dune::range(degree-k+1))
              {
                // xx
                if (degree-k-l-1 > 0)
                  y[index][0][0] = (degree-k-l)*(degree-k-l-1)*xPowers[0][degree-k-l-2]*xPowers[1][l]*xPowers[2][k];
                // xy and yx
                if (degree-k-l > 0 and l > 0){
                  y[index][0][1] = (degree-k-l)*l*xPowers[0][degree-k-l-1]*xPowers[1][l-1]*xPowers[2][k];
                  y[index][1][0] = y[index][0][1];
                }
                // yy
                if (l-1 > 0)
                  y[index][1][1] = l*(l-1)*xPowers[0][degree-k-l]*xPowers[1][l-2]*xPowers[2][k];
                // xz and zx
                if (k > 0 and degree-k-l > 0)
                {
                  y[index][0][2] = (degree-k-l)*k*xPowers[0][degree-k-l-1]*xPowers[1][l]*xPowers[2][k-1];
                  y[index][2][0] = y[index][0][2];
                }
                // yz
                if (l > 0 and k > 0)
                {
                  y[index][1][2] = l*k*xPowers[0][degree-k-l]*xPowers[1][l-1]*xPowers[2][k-1];
                  y[index][2][1] = y[index][1][2];
                }
                // zz
                if (k-1 > 0)
                  y[index][2][2] = (k-1)*k*xPowers[0][degree-k-l]*xPowers[1][l]*xPowers[2][k-2];
                ++index;
              }
            }
          }
          return y;
        }

      };

      constexpr friend auto derivative(const Derivative & d)
      {
        return Hessian{};
      }

    };

    constexpr friend auto derivative(const MonomialSet& m)
    {
      return Derivative{};
    }

  };
  #endif
} // namespace Dune::Functions

#endif // DUNE_FUNCTIONS_ANALYTICFUNCTIONS_MONOMIALSET_HH
