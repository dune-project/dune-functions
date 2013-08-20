#include <config.h>

#include <iostream>

#include <dune/common/exceptions.hh>
#include <dune/common/parallel/mpihelper.hh>

int main ( int argc, char **argv )
try
{
  Dune::MPIHelper &helper = Dune::MPIHelper::instance( argc, argv );
  std::cout << "Hello World! This is dune-functions." << std::endl;

  if( Dune::MPIHelper::isFake )
    std::cout<< "This is a sequential program." << std::endl;
  else
    std::cout << "I am rank " << helper.rank() << " of " << helper.size()
              << " processes!" << std::endl;
  return 0;
}
catch( Dune::Exception &e )
{
  std::cerr << "Dune reported error: " << e << std::endl;
}
catch(...)
{
  std::cerr << "Unknown exception thrown!" << std::endl;
}
