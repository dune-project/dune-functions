// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_FUNCTIONS_VTK_ADAPTER_HH
#define DUNE_FUNCTIONS_VTK_ADAPTER_HH

#include <dune/common/version.hh>

namespace Dune {
namespace Functions {

template<typename F>
struct VTKLocalFunction
{
private:
  using Fnkt = typename std::decay<F>::type;
  mutable Fnkt _f;
public:
  VTKLocalFunction(F && f) : _f(f) {}

  template<typename Entity>
  void bind(const Entity& e) const
  {
    _f.bind(e);
  }

  void unbind() const
  {
    _f.unbind();
  }

  template<typename X>
  auto operator()(const X & x) const -> decltype(_f(x))
  {
    return _f(x);
  }

};

/**
   \brief generate an adapter class to pass a localizable function to a VTKWriter

   \note This adapter is only needed when using dune-functions with the Dune 2.4 release.
 */
#if DUNE_VERSION_NEWER_REV(DUNE_GRID,2,4,1)
template<typename F>
auto vtkFunction(F && f) -> decltype(std::forward<F>(f))
{
  return std::forward<F>(f);
}
#else
template<typename F>
auto vtkFunction(F && f) -> VTKLocalFunction<decltype(localFunction(std::forward<F>(f)))>
{
  using LF = decltype(localFunction(std::forward<F>(f)));
  return VTKLocalFunction<LF>(localFunction(std::forward<F>(f)));
}
#endif

}
}

#endif // DUNE_FUNCTIONS_VTK_ADAPTER_HH
