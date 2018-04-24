// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_FUNCTIONS_COMMON_PERSISTENT_GRIDVIEW_HH
#define DUNE_FUNCTIONS_COMMON_PERSISTENT_GRIDVIEW_HH


#include <cstddef>
#include <memory>
#include <utility>
#include <vector>
#include <array>
#include <unordered_map>

#include <dune/common/exceptions.hh>
#include <dune/geometry/type.hh>
#include <dune/geometry/typeindex.hh>

namespace Dune {
namespace Functions {
namespace Experimental {



  /**
   * \brief A persistent index set
   *
   * \warning This is experimental and may be removed or
   * modified in a non-compatible way. When using this
   * functionality take care to follow the dune-functions
   * development to be aware of possible changes.
   *
   * A PersistentIndexSet caches the grid indices of
   * the grid view passed during creation in a persistent
   * way. The cache is filled lazily:
   *
   * After creation any accessed index is obtained from
   * the index set of the passed grid view and then stored
   * in an internal cache. Once the cache is filled the
   * PersistentIndexSet can be set to non-mutable and
   * all cached indices can be queried without
   * touching the passed grid view.
   */
  template<class GV>
  class PersistentIndexSet
  {
    using OriginalGridView = GV;
    using Grid = typename GV::Grid;

  public:

    using Types = std::vector<Dune::GeometryType>;
    using IndexType = typename GV::IndexSet::IndexType;

    template<int i>
    using Codim = typename GV::template Codim<i>;

    enum {dimension = GV::dimension};

    PersistentIndexSet(const OriginalGridView& gv) :
      originalGridView_(gv),
      grid_(&gv.grid())
    {
      setup();
    }

    // *********************************
    // size() methods
    // *********************************
    IndexType size(Dune::GeometryType gt) const
    {
      return sizePerGT_[Dune::GlobalGeometryTypeIndex::index(gt)];
    }

    IndexType size(int codim) const
    {
      return sizePerCodim_[codim];
    }

    // *********************************
    // Index methods
    // *********************************
    template<class Entity>
    IndexType index(const Entity& e) const
    {
      auto id = grid_->globalIdSet().id(e);
      auto it = idToIndex_.find(id);
      if (it!=idToIndex_.end())
        return it->second;

      if(not(isMutable_))
        DUNE_THROW(Dune::InvalidStateException, "Accessing nonexisting entry using PersistentIndexSet::index()!");

      auto index = originalGridView_.indexSet().index(e);
      idToIndex_.insert({id, index});
      return index;
    }

    template<int cc>
    IndexType index(const typename Codim<cc>::Entity& e) const
    {
      return index<typename Codim<cc>::Entity>(e);
    }

    template<class Entity>
    IndexType subIndex(const Entity& e, int subEntity, unsigned int codim) const
    {
      auto id = grid_->globalIdSet().subId(e, subEntity, codim);
      auto it = idToIndex_.find(id);
      if (it!=idToIndex_.end())
        return it->second;

      if(not(isMutable_))
        DUNE_THROW(Dune::InvalidStateException, "Accessing nonexisting entry using PersistentIndexSet::subIndex()!");

      auto index = originalGridView_.indexSet().subIndex(e, subEntity, codim);
      idToIndex_.insert({id, index});
      return index;
    }

    template<int cc>
    IndexType subIndex(const typename Codim<cc>::Entity& e, int subEntity, unsigned int codim) const
    {
      return subIndex<typename Codim<cc>::Entity>(e, subEntity, codim);
    }

    // *********************************
    // other methods
    // *********************************

    /**
     * \brief Set mutability
     *
     * In any case accessing already cached indices will
     * return the cached value.
     *
     * If the PersistentIndexSet is set mutable, any access
     * to a non-cached index will query the index form the
     * index set of the grid view passed during creation
     * and put this index into the cache.
     *
     * If the PersistentIndexSet is set non-mutable, any access
     * to a non-cached index throws an exception.
     *
     * After creation the PersistentIndexSet is set to mutable.
     */
    void setMutable(bool m) const
    {
      isMutable_ = m;
    }

    template<class Entity >
    bool contains(const Entity& e) const
    {
      auto it = idToIndex_.find(grid_->globalIdSet().id(e));
      return it != std::end(idToIndex_);
    }

    Types types(int codim) const
    {
      return typesPerCodim_[codim];
    }

  private:

    void setup()
    {
      // store global information
      const auto& originalIndexSet = originalGridView_.indexSet();
      sizePerGT_.clear();
      sizePerGT_.resize(Dune::GlobalGeometryTypeIndex::size(dimension), 0);
      for (std::size_t codim = 0; codim <= dimension; ++codim)
      {
        sizePerCodim_[codim] = originalIndexSet.size(codim);
        typesPerCodim_[codim].clear();
        for (const auto& gt : originalIndexSet.types(codim))
        {
          typesPerCodim_[codim].push_back(gt);
          auto gtIndex = Dune::GlobalGeometryTypeIndex::index(gt);
          sizePerGT_[gtIndex] = originalIndexSet.size(gt);
        }
      }

      // prepare mutable index cache
      idToIndex_.clear();
      isMutable_ = true;

      // Add all (codim-0) elements to the index map. These should be known
      // to the PersistentGridView for the contains() method.
      for(const auto& element : elements(originalGridView_)) {
        auto index = originalGridView_.indexSet().index(element);
        auto id = grid_->globalIdSet().id(element);
        idToIndex_.insert({id, index});
      }
    }

    // global information
    const OriginalGridView originalGridView_;
    const Grid* grid_;

    std::vector<std::size_t> sizePerGT_;
    std::array<std::size_t, dimension+1> sizePerCodim_;
    std::array<Types, dimension+1> typesPerCodim_;

    // index cache
    mutable std::unordered_map<typename Grid::GlobalIdSet::IdType, IndexType> idToIndex_;
    mutable bool isMutable_;
  };



  /**
   * \brief A pseudo-grid view provinding a persistent index set
   *
   * The grid view uses a PersistentIndexSet that builds
   * a cache of any accessed index for the passed grid view.
   * Those indices can later be retrieved from the internal
   * cache without accessing the passed grid view.
   *
   * Once the index set in the PersistentGridView is set
   * to non-mutable, only the cached indices that have been
   * acccessed in mutable state are accessible.
   *
   * Notice that the PersistentGridView does not provide
   * any iterators, but only the persistent index set.
   */
  template<class GV>
  class PersistentGridView
  {
    using OriginalGridView = GV;

  public:

    using Grid = typename GV::Grid;
    using ctype = typename GV::ctype;

    template<int i>
    using Codim = typename GV::template Codim<i>;

    enum {dimension = GV::dimension};
    enum {dimensionworld = GV::dimensionworld};

    using IndexSet = PersistentIndexSet<OriginalGridView>;


    PersistentGridView(OriginalGridView&& originalGridView) :
      indexSet_(std::make_shared<IndexSet>(std::move(originalGridView))),
      grid_(&originalGridView.grid())
    {}

    PersistentGridView(const OriginalGridView& originalGridView) :
      indexSet_(std::make_shared<IndexSet>(originalGridView)),
      grid_(&originalGridView.grid())
    {}

    PersistentGridView(const PersistentGridView& other) = default;

    const Grid& grid() const
    {
      return *grid_;
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
    bool contains(const Entity& e) const
    {
      return indexSet_->contains(e);
    }

  private:
    std::shared_ptr<IndexSet> indexSet_;
    const Grid* grid_;
  };



  /**
   * \brief Create a PersistentGridView
   *
   * \warning This is experimental and may be removed or
   * modified in a non-compatible way. When using this
   * functionality take care to follow the dune-functions
   * development to be aware of possible changes.
   *
   * The result provides access to the indices
   * of the passed grid view.
   */
  template<class GridView>
  auto persistentGridView(const GridView& gridView)
  {
    return PersistentGridView<GridView>(gridView);
  }



}}} // namespace Dune::Functions::Experimental

#endif// DUNE_FUNCTIONS_COMMON_PERSISTENT_GRIDVIEW_HH
