// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_FUNCTIONS_COMMON_TEST_COLLECTORSTREAM_HH
#define DUNE_FUNCTIONS_COMMON_TEST_COLLECTORSTREAM_HH

#include <sstream>
#include <string>
#include <functional>


#include <dune/functions/common/type_traits.hh>


namespace Dune {
namespace Functions {



/**
 * \brief Data collector stream
 *
 * A class derived from std::ostringstream that allows to
 * collect data via a temporary returned object. To facilitate
 * this it stores a callback that is used to pass the collected
 * data to its creator on destruction.
 *
 * In order to avoid passing the same data twice, copy construction
 * is forbidden and only move construction is allowed.
 */
class CollectorStream : public std::ostringstream
{
public:

  template<class CallBack,
    Dune::Functions::disableCopyMove<CollectorStream, CallBack> = 0>
  CollectorStream(CallBack&& callBack) :
    callBack_(callBack)
  {}

  CollectorStream(const CollectorStream& other) = delete;

  CollectorStream(CollectorStream&& other) :
    callBack_(other.callBack_)
  {
    (*this) << other.str();
    other.callBack_ = [](std::string){};
  }

  ~CollectorStream()
  {
    callBack_(this->str());
  }

private:
  std::function<void(std::string)> callBack_;
};



}} // namespace Dune::Functions


#endif // DUNE_FUNCTIONS_COMMON_TEST_COLLECTORSTREAM_HH
