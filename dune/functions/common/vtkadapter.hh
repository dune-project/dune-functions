// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_FUNCTIONS_VTK_ADAPTER_HH
#define DUNE_FUNCTIONS_VTK_ADAPTER_HH

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
   \brief generate an adapter class to pass a localozable function to a VTKWriter

   \note this adapter is only needed when using dune-functions with then Dune 2.4 release
 */
template<typename F>
auto vtkFunction(F && f) -> VTKLocalFunction<decltype(localFunction(std::forward<F>(f)))>
{
  using LF = decltype(localFunction(std::forward<F>(f)));
  return VTKLocalFunction<LF>(localFunction(std::forward<F>(f)));
}

}
}

#endif // DUNE_FUNCTIONS_VTK_ADAPTER_HH
