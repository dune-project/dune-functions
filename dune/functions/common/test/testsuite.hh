// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_FUNCTIONS_COMMON_TEST_TESTSUITE_HH
#define DUNE_FUNCTIONS_COMMON_TEST_TESTSUITE_HH

#include <iostream>
#include <sstream>
#include <string>

#include <dune/common/exceptions.hh>

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
  enum ThrowPolicy
  {
    AlwaysThrow,
    ThrowOnRequired
  };

  TestSuite(ThrowPolicy policy, std::string name="") :
    name_(name),
    checks_(0),
    failedChecks_(0),
    throwPolicy_(policy==AlwaysThrow)
  {}

  TestSuite(std::string name="", ThrowPolicy policy=ThrowOnRequired) :
    name_(name),
    checks_(0),
    failedChecks_(0),
    throwPolicy_(policy==AlwaysThrow)
  {}

  CollectorStream check(bool condition, std::string name="")
  {
    ++checks_;
    if (not condition)
      ++failedChecks_;

    return CollectorStream([=](std::string reason) {
        if (not condition)
          this->announceCheckResult(throwPolicy_, "CHECK  ", name, reason);
      });
  }

  CollectorStream require(bool condition, std::string name="")
  {
    ++checks_;
    if (not condition)
      ++failedChecks_;

    return CollectorStream([=](std::string reason) {
        if (not condition)
          this->announceCheckResult(true, "REQUIRED CHECK", name, reason);
      });
  }

  void subTest(const TestSuite& subTest)
  {
    checks_ += subTest.checks_;
    failedChecks_ += subTest.failedChecks_;

    if (not subTest)
      announceCheckResult(throwPolicy_, "SUBTEST", subTest.name(), std::to_string(subTest.failedChecks_)+"/"+std::to_string(subTest.checks_) + " checks failed in this subtest.");
  }

  operator const bool () const
  {
    return (failedChecks_==0);
  }

  std::string name() const
  {
    return name_;
  }

  bool report() const
  {
    if (failedChecks_>0)
      std::cout << composeMessage("TEST   ", name(), std::to_string(failedChecks_)+"/"+std::to_string(checks_) + " checks failed in this test.") << std::endl;
    return (failedChecks_==0);
  }

  int exit() const
  {
    return (report() ? 0: 1);
  }



protected:

  static std::string composeMessage(std::string type, std::string name, std::string reason)
  {
    std::ostringstream s;
    s << type << " FAILED";
    if (name!="")
      s << "(" << name << ")";
    s << ": ";
    if (reason!="")
      s << reason;
    return s.str();
  }

  static void announceCheckResult(bool throwException, std::string type, std::string name, std::string reason)
  {
    std::string message = composeMessage(type, name, reason);
    std::cout << message << std::endl;
    if (throwException)
    {
      Dune::Exception ex;
      ex.message(message);
      throw ex;
    }
  }

  std::string name_;
  std::size_t checks_;
  std::size_t failedChecks_;
  bool throwPolicy_;
};



}} // namespace Dune::Functions


#endif // DUNE_FUNCTIONS_COMMON_TEST_TESTSUITE_HH
