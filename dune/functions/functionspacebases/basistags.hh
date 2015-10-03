// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_FUNCTIONS_FUNCTIONSPACEBASES_BASISTAGS_HH
#define DUNE_FUNCTIONS_FUNCTIONSPACEBASES_BASISTAGS_HH



namespace Dune {
namespace Functions {
namespace BasisTags {



struct IndexTag {};
struct FlatIndex : public IndexTag {};
struct InterleafedIndex : public IndexTag {};
struct BlockedIndex : public IndexTag {};
struct LeafBlockedIndex : public IndexTag {};



} // end namespace BasisTags
} // end namespace Functions
} // end namespace Dune


#endif // DUNE_FUNCTIONS_FUNCTIONSPACEBASES_BASISTAGS_HH
