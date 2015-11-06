// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_FUNCTIONS_COMMON_TEST_TESTSUITE_HH
#define DUNE_FUNCTIONS_COMMON_TEST_TESTSUITE_HH

#include <iostream>
#include <sstream>
#include <string>

#include <dune/common/exceptions.hh>

#include <dune/functions/common/type_traits.hh>

#include <dune/functions/common/test/collectorstream.hh>


namespace Dune {
namespace Functions {



/**
 * \brief A Simple helper class to organize your test suite
 *
 * Usage: Construct a TestSuite and call check() or require()
 * with the condition to check and probably a name for this check.
 * These methods return a stream such that you can pipe in an
 * explanantion accomponied by respective data to give a reason
 * for a test failure.
 */
class TestSuite
{
public:

  TestSuite(bool alwaysThrow, std::string name="") :
    name_(name),
    result_(true),
    alwaysThrow_(alwaysThrow)
  {}

  TestSuite(std::string name="", bool alwaysThrow=false) :
    name_(name),
    result_(true),
    alwaysThrow_(alwaysThrow)
  {}

  CollectorStream check(bool condition, std::string name="")
  {
    return CollectorStream([=](std::string reason) {
        this->announceCheckResult(condition, false, "CHECK  ", name, reason);
      });
  }

  CollectorStream require(bool condition, std::string name="")
  {
    return CollectorStream([=](std::string reason) {
        this->announceCheckResult(condition, true, "REQUIRED CHECK", name, reason);
      });
  }

  void subTest(const TestSuite& subTest)
  {
    announceCheckResult(subTest, false, "SUBTEST", subTest.name(), "Some checks or subtests of this subtest failed.");
  }

  operator const bool () const
  {
    return result_;
  }

  std::string name() const
  {
    return name_;
  }

  bool report() const
  {
    if (not result_)
      std::cout << composeMessage(false, "TEST   ", name_, "Some checks or subtests of the test failed.") << std::endl;
    return result_;
  }

  int exit() const
  {
    return (report() ? 0: 1);
  }



protected:

  std::string composeMessage(bool required, std::string type, std::string name, std::string reason) const
  {
    std::ostringstream s;
    if (required)
      s << "REQUIRED ";
    s << type << " FAILED : ";
    if (name!="")
      s << name << " ";
    if (reason!="")
      s << "REASON : " << reason;
    return s.str();
  }

  void announceCheckResult(bool condition, bool required, std::string type, std::string name, std::string reason)
  {
    result_ &= condition;
    if (not condition)
    {
      std::string message = composeMessage(required, type, name, reason);
      std::cout << message << std::endl;
      if (required or alwaysThrow_)
        DUNE_THROW(Dune::Exception, message);
    }
  }

  std::string name_;
  bool result_;
  bool alwaysThrow_;
};



}} // namespace Dune::Functions


#endif // DUNE_FUNCTIONS_COMMON_TEST_TESTSUITE_HH
