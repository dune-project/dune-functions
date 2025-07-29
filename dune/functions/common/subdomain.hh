// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

// SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file AUTHORS.md
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception OR LGPL-3.0-or-later

#ifndef DUNE_FUNCTIONS_COMMON_SUBDOMAIN_HH
#define DUNE_FUNCTIONS_COMMON_SUBDOMAIN_HH

#include <array>
#include <cstddef>
#include <limits>
#include <type_traits>
#include <utility>
#include <vector>

#include <dune/common/exceptions.hh>
#include <dune/common/iteratorfacades.hh>
#include <dune/common/iteratorrange.hh>
#include <dune/common/rangeutilities.hh>

#include <dune/geometry/type.hh>
#include <dune/geometry/typeindex.hh>

#include <dune/grid/common/mcmgmapper.hh>
#include <dune/grid/common/partitionset.hh>

namespace Dune::Functions::Experimental {

  namespace Impl {

    template<class GridView>
    using GlobalIntersectionIteratorTraits = Dune::DefaultIteratorTraits<
                                                std::forward_iterator_tag,
                                                decltype(*std::declval<typename GridView::IntersectionIterator>())>;


    template<class GV, class ContainsCallback>
    class GlobalIntersectionIt
      : public Dune::IteratorFacadeForTraits<GlobalIntersectionIt<GV, ContainsCallback>, GlobalIntersectionIteratorTraits<GV>>
    {
      using Facade = Dune::IteratorFacadeForTraits<GlobalIntersectionIt<GV, ContainsCallback>, GlobalIntersectionIteratorTraits<GV>>;

    public:

      using GridView = GV;
      using Element = typename GridView::template Codim<0>::Entity;
      using ElementIterator = typename GridView::template Codim<0>::Iterator;
      using IntersectionIterator = typename GridView::IntersectionIterator;

      class SentinelIterator
      {};

      GlobalIntersectionIt(const GridView& gridView, const ContainsCallback& contains, ElementIterator elementIt, ElementIterator elementEnd)
        : gridView_(gridView)
        , contains_(contains)
        , elementIt_(std::move(elementIt))
        , elementEnd_(std::move(elementEnd))
      {
        if (elementIt_ != elementEnd_)
        {
          element_ = *elementIt_;
          iIt_ = gridView_.ibegin(element_);
          iEnd_ = gridView_.iend(element_);
          if (not contains_(*iIt_))
            ++(*this);
        }
      }

      using reference = typename Facade::reference;

      reference operator*() const
      {
        return *iIt_;
      }

      GlobalIntersectionIt& operator++()
      {
        while(true)
        {
          ++iIt_;
          if (iIt_ == iEnd_)
          {
            ++elementIt_;
            if (elementIt_ == elementEnd_)
              return *this;
            element_ = *elementIt_;
            iIt_ = gridView_.ibegin(element_);
            iIt_ = gridView_.ibegin(element_);
            iEnd_ = gridView_.iend(element_);
          }
          if (contains_(*iIt_))
            return *this;
        }
        return *this;
      }

      friend bool operator==(const GlobalIntersectionIt& it1, const GlobalIntersectionIt& it2)
      {
        if (it1.elementIt_ != it2.elementIt_)
          return false;
        if (it1.elementIt_ == it1.elementEnd_)
          return true;
        return (it1.iIt_ == it2.iIt_);
      }

      friend bool operator==(const GlobalIntersectionIt& it1, const SentinelIterator& it2)
      {
        return it1.elementIt_ == it1.elementEnd_;
      }

    private:
      GridView gridView_;
      ContainsCallback contains_;
      ElementIterator elementIt_;
      ElementIterator elementEnd_;
      Element element_;
      IntersectionIterator iIt_;
      IntersectionIterator iEnd_;
    };

  }



  /**
   * \brief An IndexSet for a sub-domain
   *
   * \ingroup Utility
   *
   * A SubDomainIndexSet provides indices for a subset of the
   * entities of a grid view. Once created, new entities can be
   * inserted using `subDomainIndexSet.insertElement(element)`
   * which will insert the grid `element` of codim 0 and all
   * its sub entities. This will increment the internally stored
   * size information accordingly. Querying the index of an
   * entity or sub-entity that has not been inserted before
   * will lead to an exception.
   *
   * Internally this uses a `std::vector<std::size_t>` which stores
   * the sub-domain index for each entity in the underlying grid view
   * using a magic value for entities not contained in the sub-domain.
   *
   * \tparam HGV The underlying host grid view type.
   */
  template<class HGV>
  class SubDomainIndexSet
  {
    using HostGridView = HGV;

  public:

    using Grid = typename HostGridView::Grid;
    using Types = std::vector<Dune::GeometryType>;
    using IndexType = std::size_t;

    //! Codim specific typedefs
    template<int codim>
    struct Codim
    {
      using Entity = typename Grid::template Codim<codim>::Entity;
      using EntitySeed = typename Grid::template Codim<codim>::EntitySeed;
      using Geometry = typename Grid::template Codim<codim>::Geometry;
      using LocalGeometry = typename Grid::template Codim<codim>::LocalGeometry;
    };

    enum {dimension = Grid::dimension};

  private:

    using AllEntityMapper = Dune::MultipleCodimMultipleGeomTypeMapper<HostGridView>;

    static auto allCodimLayout()
    {
      return [](Dune::GeometryType, int) { return true; };
    }

    static constexpr auto typeIndexSize = Dune::GlobalGeometryTypeIndex::size(dimension);
    static constexpr auto unusesIndex = std::numeric_limits<IndexType>::max();

  public:

    //! Construct SubDomainIndexSet for underlying host grid view
    SubDomainIndexSet(const HostGridView& hostGridView)
      : hostGridView_(hostGridView)
      , allEntityMapper_(hostGridView_, allCodimLayout())
    {
      clear();
    }

    // *********************************
    // IndexSet interface methods
    // *********************************

    IndexType size(Dune::GeometryType gt) const
    {
      return sizePerGT_[Dune::GlobalGeometryTypeIndex::index(gt)];
    }

    IndexType size(int codim) const
    {
      return sizePerCodim_[codim];
    }

    template<class Entity>
    IndexType index(const Entity& entity) const
    {
      auto index = indices_[allEntityMapper_.index(entity)];
      if (index==unusesIndex)
        DUNE_THROW(Dune::InvalidStateException, "Accessing nonexisting entry using SubDomainIndexSet::index()!");
      return index;
    }

    template<int cc>
    IndexType index(const typename Codim<cc>::Entity& entity) const
    {
      return index<typename Codim<cc>::Entity>(entity);
    }

    template<class Entity>
    IndexType subIndex(const Entity& entity, int subEntity, unsigned int codim) const
    {
      auto index = indices_[allEntityMapper_.subIndex(entity, subEntity, codim)];
      if (index==unusesIndex)
        DUNE_THROW(Dune::InvalidStateException, "Accessing nonexisting entry using SubDomainIndexSet::subIndex()!");
      return index;
    }

    template<int cc>
    IndexType subIndex(const typename Codim<cc>::Entity& entity, int subEntity, unsigned int codim) const
    {
      return subIndex<typename Codim<cc>::Entity>(entity, subEntity, codim);
    }

    template<class Entity >
    bool contains(const Entity& entity) const
    {
      return (indices_[allEntityMapper_.index(entity)] != unusesIndex);
    }

    Types types(int codim) const
    {
      return typesPerCodim_[codim];
    }

    // *********************************
    // Extended methods
    // *********************************

    //! Access underlying host grid view
    const HostGridView& hostGridView() const
    {
      return hostGridView_;
    }

    //! Insert element and all its sub-entities into SubDomainIndexSet
    void insertElement(const typename Codim<0>::Entity& element)
    {
      const auto& re = referenceElement(element);
      for (auto codim : Dune::range(0, dimension+1))
      {
        for (auto subEntity : Dune::range(re.size(codim)))
        {
          auto& index = indices_[allEntityMapper_.subIndex(element, subEntity, codim)];
          if (index==unusesIndex)
          {
            const auto& type = re.type(subEntity, codim);
            const auto typeIndex = Dune::GlobalGeometryTypeIndex::index(type);
            index = sizePerGT_[typeIndex];
            if (sizePerGT_[typeIndex]==0)
              typesPerCodim_[codim].push_back(type);
            sizePerGT_[typeIndex]++;
            sizePerCodim_[codim]++;
          }
        }
      }
    }

  protected:

    // Clear all data
    void clear()
    {
      for(auto& size : sizePerGT_)
        size = 0;
      for(auto& size : sizePerCodim_)
        size = 0;
      for(auto& types : typesPerCodim_)
        types.clear();
      indices_.clear();
      indices_.resize(allEntityMapper_.size(), unusesIndex);
    }

    HostGridView hostGridView_;

    // Global size information
    std::array<std::size_t, typeIndexSize> sizePerGT_;
    std::array<std::size_t, dimension+1> sizePerCodim_;
    std::array<Types, dimension+1> typesPerCodim_;

    AllEntityMapper allEntityMapper_;

    // Index map
    std::vector<IndexType> indices_;
  };



  /**
   * \brief A GridView for a sub-domain
   *
   * \ingroup Utility
   *
   * A SubDomainGridView provides a reduces GridView interface
   * for a subset of the entities of a grid view. Objects of
   * this class store a pointer to a SubDomainIndexSet and
   * can be copied cheaply.
   *
   * \tparam HGV The underlying host grid view type.
   */
  template<class HGV>
  class SubDomainGridView
  {

    template<int codim>
    class NonImplementedIterator
    {
    public:
      NonImplementedIterator()
      {
        static_assert(codim==0, "SubDomainGridView::Codim::Iterator<codim> is only implemented for codim=0");
      }
    };

    template<PartitionIteratorType pit>
    class ElementIterator
    {
      using Element = typename HGV::template Codim<0>::Entity;
    public:

      using HostElementIterator = typename HGV::template Codim<0>::template Partition<pit>::Iterator;

      ElementIterator(const SubDomainIndexSet<HGV>& indexSet, HostElementIterator&& it, HostElementIterator&& endIt)
        : indexSet_(&indexSet)
        , hostIt_(std::move(it))
        , hostEndIt_(std::move(endIt))
      {
        while ((hostIt_!= hostEndIt_) and (not indexSet_->contains(*hostIt_)))
          ++hostIt_;
      }

      ElementIterator& operator++()
      {
        ++hostIt_;
        while ((hostIt_!= hostEndIt_) and (not indexSet_->contains(*hostIt_)))
          ++hostIt_;
        return *this;
      }

      const Element& operator*() const
      {
        return *hostIt_;
      }

      friend bool operator==(const ElementIterator& a, const ElementIterator& b)
      {
        return a.hostIt_==b.hostIt_;
      }

    private:
      HostElementIterator hostIt_;
      HostElementIterator hostEndIt_;
      const SubDomainIndexSet<HGV>* indexSet_;
    };

  public:

    using HostGridView = HGV;

    using Grid = typename HostGridView::Grid;
    using ctype = typename Grid::ctype;
    using IndexSet = SubDomainIndexSet<HostGridView>;
    using Intersection = typename HostGridView::Intersection;
    using IntersectionIterator = typename HostGridView::IntersectionIterator;

    //! Codim specific typedefs
    template<int codim>
    struct Codim
    {
      using Entity = typename Grid::template Codim<codim>::Entity;
      using EntitySeed = typename Grid::template Codim<codim>::EntitySeed;
      using Geometry = typename Grid::template Codim<codim>::Geometry;
      using LocalGeometry = typename Grid::template Codim<codim>::LocalGeometry;
      using Iterator = std::conditional_t<codim==0, ElementIterator<All_Partition>, NonImplementedIterator<codim>>;

      template<PartitionIteratorType pit>
      struct Partition
      {
        using Iterator = std::conditional_t<codim==0, ElementIterator<pit>, NonImplementedIterator<codim>>;
      };
    };

    enum {dimension = Grid::dimension};
    enum {dimensionworld = Grid::dimensionworld};

    SubDomainGridView(const IndexSet& indexSet)
      : indexSet_(&indexSet)
    {}

    SubDomainGridView(const SubDomainGridView& other) = default;

    const Grid& grid() const
    {
      return indexSet_->hostGridView().grid();
    }

    const IndexSet& indexSet() const
    {
      return *indexSet_;
    }

    int size(int codim) const
    {
      return indexSet_->size(codim);
    }

    int size(Dune::GeometryType gt) const
    {
      return indexSet_->size(gt);
    }

    template<class Entity>
    bool contains(const Entity& entity) const
    {
      return indexSet_->contains(entity);
    }

    //! Create an iterator pointing to the begin of the range.
    template<int codim, PartitionIteratorType pit = All_Partition>
    typename Codim<codim>::template Partition<pit>::Iterator begin() const
    {
      static_assert(codim==0, "SubDomainGridView::begin<codim> is only implemented for codim=0");
      return {indexSet(), hostGridView().template begin<codim, pit>(), hostGridView().template end<codim, pit>()};
    }

    //! Create an iterator pointing to the end of the range.
    template<int codim, PartitionIteratorType pit = All_Partition>
    typename Codim<codim>::template Partition<pit>::Iterator end() const
    {
      static_assert(codim==0, "SubDomainGridView::end<codim> is only implemented for codim=0");
      return {indexSet(), hostGridView().template end<codim, pit>(), hostGridView().template end<codim, pit>()};
    }

    decltype(auto) comm() const
    {
      return hostGridView().comm();
    }

    decltype(auto) ibegin(const typename Codim<0>::Entity& element) const
    {
      return hostGridView().ibegin(element);
    }

    decltype(auto) iend(const typename Codim<0>::Entity& element) const
    {
      return hostGridView().iend(element);
    }

    //! Access underlying host grid view
    const HostGridView& hostGridView() const
    {
      return indexSet_->hostGridView();
    }

  protected:
    const IndexSet* indexSet_;
  };



  /**
   * \brief ADL findable access to element range for a SubDomainGridView
   *
   * \ingroup Utility
   */
  template<class HostGridView>
  auto elements(const SubDomainGridView<HostGridView>& subDomainGridView)
  {
    return Dune::IteratorRange(subDomainGridView.template begin<0>(), subDomainGridView.template end<0>());
  }

  /**
   * \brief ADL findable access to element range for a SubDomainGridView
   *
   * \ingroup Utility
   */
  template<class HostGridView, unsigned int partitions>
  auto elements(const SubDomainGridView<HostGridView>& subDomainGridView, Dune::PartitionSet<partitions> partitionSet)
  {
    constexpr auto pit = partitionSet.partitionIterator();
    return Dune::IteratorRange(subDomainGridView.template begin<0, pit>(), subDomainGridView.template end<0, pit>());
  }

  /**
   * \brief ADL findable access to intersection range for an element of a SubDomainGridView
   *
   * \ingroup Utility
   */
  template<class HostGridView, class Element>
  auto intersections(const SubDomainGridView<HostGridView>& subDomainGridView, const Element& element)
  {
    return Dune::IteratorRange(subDomainGridView.ibegin(element), subDomainGridView.iend(element));
  }



  /**
   * \brief Class representing a sub-domain of a GridView
   *
   * \ingroup Utility
   *
   * A SubDomain is a subset of grid elements from a given
   * underlying grid view together with their sub-entities.
   * It allows to create a `SubDomainGridView` which implements
   * a reasonable subset of the grid view interface defined
   * in dune-grid. In particular the `SubDomainGridView`
   * implements an index set, a `contains()` methods, an
   * element iterator and intersection iterators.
   *
   * \tparam HGV The underlying host grid view type.
   */
  template<class HGV>
  class SubDomain
  {
  public:

    using HostGridView = HGV;
    using Grid = typename HostGridView::Grid;
    using IndexSet = SubDomainIndexSet<HostGridView>;
    using GridView = SubDomainGridView<HostGridView>;

    //! Codim specific typedefs
    template<int codim>
    struct Codim
    {
      using Entity = typename Grid::template Codim<codim>::Entity;
      using EntitySeed = typename Grid::template Codim<codim>::EntitySeed;
      using Geometry = typename Grid::template Codim<codim>::Geometry;
      using LocalGeometry = typename Grid::template Codim<codim>::LocalGeometry;
    };

    enum {dimension = Grid::dimension};

    //! Construct SubDomain for underlying host grid view
    SubDomain(const HostGridView& hostGridView)
      : indexSet_(hostGridView)
    {}

    const IndexSet& indexSet() const
    {
      return indexSet_;
    }

    //! Create grid view representing the SubDomain
    GridView gridView() const
    {
      return GridView(indexSet_);
    }

    //! Access underlying host grid view
    HostGridView hostGridView() const
    {
      return indexSet_.hostGridView();
    }

    //! Insert element and all its sub-entities into SubDomain
    void insertElement(const typename Codim<0>::Entity& element)
    {
      indexSet_.insertElement(element);
    }

    //! Check if element is contained in SubDomain
    bool contains(const typename Codim<0>::Entity& element) const
    {
      return indexSet_.contains(element);
    }

  private:
    IndexSet indexSet_;
  };



  /**
   * \brief Class representing the intersection between two subdomains
   *
   * Conceptually this represents a range of intersections.
   */
  template<class SubDomainA, class SubDomainB>
  class SubDomainInterface
  {
    static_assert(
      std::is_same_v<typename SubDomainA::GridView::Intersection, typename SubDomainB::GridView::Intersection>,
      "SubDomainInterface requires that both SubDomain types have the same Intersection type");

  public:

    using Intersection = typename SubDomainA::GridView::Intersection;

    /**
     * \brief Create interface between two subdomains
     *
     * Notice that the order of the passed subdomains does matter,
     * because the intersections visited by the iterator will
     * always be oriented such that A is inside and B outside.
     */
    SubDomainInterface(const SubDomainA& subDomainA, const SubDomainB& subDomainB)
      : subDomainA_(subDomainA)
      , subDomainB_(subDomainB)
    {}

    /**
     * \brief Check if intersection is contained in the interface between the subdomains
     *
     * This evaluates to true for all intersections having element from subdomain A and B
     * as inside and outside or outside and inside, respectively.
     *
     * Notice that this also evaluates to true for intersection with B as inside and
     * A as outside, despite the fact that these are not visited by the iterator.
     */
    bool contains(const Intersection& is) const
    {
      if (is.boundary() or not(is.neighbor()))
        return false;
      return (subDomainA_.contains(is.inside()) && subDomainB_.contains(is.outside()))
          || (subDomainA_.contains(is.outside()) && subDomainB_.contains(is.inside()));
    }

    /**
     * \brief Check if intersection is oriented
     *
     * The intersection is called oriented if its inside and outside
     * element are contained in subDomainA and subDomainB, respectively.
     * The result of this method is undefined if the intersection is not
     * contained in the intersection.
     */
    bool isOriented(const Intersection& is) const
    {
      if (is.boundary() or not(is.neighbor()))
        return false;
      return (subDomainA_.contains(is.inside()) && subDomainB_.contains(is.outside()));
    }

    /**
     * \brief Begin iterator over all intersection between the subdomains
     *
     * The iterator will always have the elements from sub domain A and B
     * as inside and outside, respectively.
     */
    const auto begin() const
    {
      return Impl::GlobalIntersectionIt(subDomainA_.gridView(), [&](const auto& is) {
        if (is.boundary() or not(is.neighbor()))
          return false;
        return subDomainB_.indexSet().contains(is.outside());
      }, subDomainA_.gridView().template begin<0>(), subDomainA_.gridView().template end<0>());
    }

    //! End iterator (sentinel)
    const auto end() const
    {
      return typename decltype(begin())::SentinelIterator();
    }

  private:
    const SubDomainA& subDomainA_;
    const SubDomainB& subDomainB_;
  };



  /**
   * \brief Class representing the skeleton of a subdomain
   *
   * Conceptually this represents a range of intersections.
   */
  template<class SubDomain>
  class SubDomainSkeleton
  {
  public:

    using Intersection = typename SubDomain::GridView::Intersection;

    //! Create skeleton of a subdomain
    SubDomainSkeleton(const SubDomain& subDomain)
      : subDomain_(subDomain)
    {}

    //! Check if intersection is contained in the skeleton of the subdomain
    bool contains(const Intersection& is) const
    {
      if (is.boundary() or not(is.neighbor()))
        return false;
      return subDomain_.contains(is.inside()) and subDomain_.contains(is.outside());
    }

  private:
    const SubDomain& subDomain_;
  };



} // namespace Dune::Functions::Experimental

#endif// DUNE_FUNCTIONS_COMMON_SUBDOMAIN_HH
