#include <iostream>

#include <dune/grid/yaspgrid.hh>

#include <dune/functions/common/access.hh>
#include <dune/functions/functionspacebases/blocking.hh>
#include <dune/functions/functionspacebases/powerbasis.hh>
#include <dune/functions/functionspacebases/compositebasis.hh>
#include <dune/functions/functionspacebases/lagrangebasis.hh>
#include <dune/functions/functionspacebases/subspacebasis.hh>

#include "taylorhoodbases.hh"

using namespace Dune::Functions;
using namespace Dune::Functions::Blocking;

template <class Blocking>
void printBlocking()
{
  std::cout << __PRETTY_FUNCTION__ << "\n";
}

template <class GlobalBasis>
void printType(GlobalBasis const&)
{
  using Blocking = Blocking_t<GlobalBasis>;
  printBlocking<Blocking>();
}


int main(int argc, char** argv)
{
  Dune::MPIHelper::instance(argc, argv);

  TaylorHoodBases tester{};

  static const std::size_t K = 1; // pressure order for Taylor-Hood

  { // Root: blockedLexicographic, Velocity: flatLexicographic
    auto taylorHoodBasis = tester.basis1();
    using GlobalBasis = decltype(taylorHoodBasis);
    using Blocking = Blocking_t<GlobalBasis>;

    static_assert(std::is_same<Blocking, Blocked<Flat,Flat>>::value, "");
    printType(taylorHoodBasis);
  }


  { // Root: blockedLexicographic, Velocity: flatInterleaved
    auto taylorHoodBasis = tester.basis2();
    using GlobalBasis = decltype(taylorHoodBasis);
    using Blocking = Blocking_t<GlobalBasis>;

    static_assert(std::is_same<Blocking, Blocked<Flat,Flat>>::value, "");
    printType(taylorHoodBasis);
  }

  { // Root: blockedLexicographic, Velocity: blockedLexicographic
    auto taylorHoodBasis = tester.basis3();
    using GlobalBasis = decltype(taylorHoodBasis);
    using Blocking = Blocking_t<GlobalBasis>;

    static_assert(std::is_same<Blocking, Blocked<Blocked<Flat,Flat>,Flat>>::value, "");
    printType(taylorHoodBasis);
  }

  { // Root: blockedLexicographic, Velocity: leafBlockedInterleaved
    auto taylorHoodBasis = tester.basis4();
    using GlobalBasis = decltype(taylorHoodBasis);
    using Blocking = Blocking_t<GlobalBasis>;

    static_assert(std::is_same<Blocking, Blocked<LeafBlocked<2>,Flat>>::value, "");
    printType(taylorHoodBasis);
  }

  { // Root: flatLexicographic, Velocity/Pressure: blockedLexicographic
    auto taylorHoodBasis = tester.basis5();
    using GlobalBasis = decltype(taylorHoodBasis);
    using Blocking = Blocking_t<GlobalBasis>;

    static_assert(std::is_same<Blocking, Unknown>::value, "");
    printType(taylorHoodBasis);
  }


  { // Root: flatLexicographic, Velocity/Pressure: flatLexicographic
    auto taylorHoodBasis = tester.basis6();
    using GlobalBasis = decltype(taylorHoodBasis);
    using Blocking = Blocking_t<GlobalBasis>;

    static_assert(std::is_same<Blocking, Flat>::value, "");
    printType(taylorHoodBasis);
  }
}
