# .. cmake_module::
#
#    Module that checks for supported C++14, C++11 and non-standard features.
#
#    The behaviour of this module can be modified by the following variable:
#
#    :ref:`DISABLE_CXX_VERSION_CHECK`
#       Disable checking for std=c++11 (c++14, c++1y)
#
#    This module internally sets the following variables, which are then
#    exported into the config.h of the current dune module.
#
#    :code:`HAS_ATTRIBUTE_UNUSED`
#       True if attribute unused is supported
#
#    :code:`HAS_ATTRIBUTE_DEPRECATED`
#       True if attribute deprecated is supported
#
#    :code:`HAS_ATTRIBUTE_DEPRECATED_MSG`
#       True if attribute deprecated("msg") is supported
#
# .. cmake_variable:: DISABLE_CXX_VERSION_CHECK
#
#    You may set this variable to TRUE to disable checking for
#    std=c++11 (c++14, c++1y). For more details, check :ref:`CheckCXXFeatures`.
#

# COPYRIGHT AND LICENSE:
#
# This file has been backported from the dune-common library. It is not
# licensed under the same conditions as the rest of dune-functions. For
# the list of copyright holders and license of this file see below.
#
# Copyright holders of this file:
# ===============================
#
# Markus Blatt
# Christian Engwer
# Jorrit Fahlke
# Christoph Grüninger
# Dominic Kempf
#
# License of this file:
# =====================
#
# The DUNE library and headers are licensed under version 2 of the GNU
# General Public License (see below), with a special exception for
# linking and compiling against DUNE, the so-called "runtime exception."
# The license is intended to be similar to the GNU Lesser General
# Public License, which by itself isn't suitable for a template library.
#
# The exact wording of the exception reads as follows:
#
#    As a special exception, you may use the DUNE source files as part
#    of a software library or application without restriction.
#    Specifically, if other files instantiate templates or use macros or
#    inline functions from one or more of the DUNE source files, or you
#    compile one or more of the DUNE source files and link them with
#    other files to produce an executable, this does not by itself cause
#    the resulting executable to be covered by the GNU General Public
#    License.  This exception does not however invalidate any other
#    reasons why the executable file might be covered by the GNU General
#    Public License.
#
#
#         GNU GENERAL PUBLIC LICENSE
#            Version 2, June 1991
#
#  Copyright (C) 1989, 1991 Free Software Foundation, Inc.
#  51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA
#  Everyone is permitted to copy and distribute verbatim copies
#  of this license document, but changing it is not allowed.
#
#           Preamble
#
#   The licenses for most software are designed to take away your
# freedom to share and change it.  By contrast, the GNU General Public
# License is intended to guarantee your freedom to share and change free
# software--to make sure the software is free for all its users.  This
# General Public License applies to most of the Free Software
# Foundation's software and to any other program whose authors commit to
# using it.  (Some other Free Software Foundation software is covered by
# the GNU Library General Public License instead.)  You can apply it to
# your programs, too.
#
#   When we speak of free software, we are referring to freedom, not
# price.  Our General Public Licenses are designed to make sure that you
# have the freedom to distribute copies of free software (and charge for
# this service if you wish), that you receive source code or can get it
# if you want it, that you can change the software or use pieces of it
# in new free programs; and that you know you can do these things.
#
#   To protect your rights, we need to make restrictions that forbid
# anyone to deny you these rights or to ask you to surrender the rights.
# These restrictions translate to certain responsibilities for you if you
# distribute copies of the software, or if you modify it.
#
#   For example, if you distribute copies of such a program, whether
# gratis or for a fee, you must give the recipients all the rights that
# you have.  You must make sure that they, too, receive or can get the
# source code.  And you must show them these terms so they know their
# rights.
#
#   We protect your rights with two steps: (1) copyright the software, and
# (2) offer you this license which gives you legal permission to copy,
# distribute and/or modify the software.
#
#   Also, for each author's protection and ours, we want to make certain
# that everyone understands that there is no warranty for this free
# software.  If the software is modified by someone else and passed on, we
# want its recipients to know that what they have is not the original, so
# that any problems introduced by others will not reflect on the original
# authors' reputations.
#
#   Finally, any free program is threatened constantly by software
# patents.  We wish to avoid the danger that redistributors of a free
# program will individually obtain patent licenses, in effect making the
# program proprietary.  To prevent this, we have made it clear that any
# patent must be licensed for everyone's free use or not licensed at all.
#
#   The precise terms and conditions for copying, distribution and
# modification follow.
#
#         GNU GENERAL PUBLIC LICENSE
#    TERMS AND CONDITIONS FOR COPYING, DISTRIBUTION AND MODIFICATION
#
#   0. This License applies to any program or other work which contains
# a notice placed by the copyright holder saying it may be distributed
# under the terms of this General Public License.  The "Program", below,
# refers to any such program or work, and a "work based on the Program"
# means either the Program or any derivative work under copyright law:
# that is to say, a work containing the Program or a portion of it,
# either verbatim or with modifications and/or translated into another
# language.  (Hereinafter, translation is included without limitation in
# the term "modification".)  Each licensee is addressed as "you".
#
# Activities other than copying, distribution and modification are not
# covered by this License; they are outside its scope.  The act of
# running the Program is not restricted, and the output from the Program
# is covered only if its contents constitute a work based on the
# Program (independent of having been made by running the Program).
# Whether that is true depends on what the Program does.
#
#   1. You may copy and distribute verbatim copies of the Program's
# source code as you receive it, in any medium, provided that you
# conspicuously and appropriately publish on each copy an appropriate
# copyright notice and disclaimer of warranty; keep intact all the
# notices that refer to this License and to the absence of any warranty;
# and give any other recipients of the Program a copy of this License
# along with the Program.
#
# You may charge a fee for the physical act of transferring a copy, and
# you may at your option offer warranty protection in exchange for a fee.
#
#   2. You may modify your copy or copies of the Program or any portion
# of it, thus forming a work based on the Program, and copy and
# distribute such modifications or work under the terms of Section 1
# above, provided that you also meet all of these conditions:
#
#     a) You must cause the modified files to carry prominent notices
#     stating that you changed the files and the date of any change.
#
#     b) You must cause any work that you distribute or publish, that in
#     whole or in part contains or is derived from the Program or any
#     part thereof, to be licensed as a whole at no charge to all third
#     parties under the terms of this License.
#
#     c) If the modified program normally reads commands interactively
#     when run, you must cause it, when started running for such
#     interactive use in the most ordinary way, to print or display an
#     announcement including an appropriate copyright notice and a
#     notice that there is no warranty (or else, saying that you provide
#     a warranty) and that users may redistribute the program under
#     these conditions, and telling the user how to view a copy of this
#     License.  (Exception: if the Program itself is interactive but
#     does not normally print such an announcement, your work based on
#     the Program is not required to print an announcement.)
#
# These requirements apply to the modified work as a whole.  If
# identifiable sections of that work are not derived from the Program,
# and can be reasonably considered independent and separate works in
# themselves, then this License, and its terms, do not apply to those
# sections when you distribute them as separate works.  But when you
# distribute the same sections as part of a whole which is a work based
# on the Program, the distribution of the whole must be on the terms of
# this License, whose permissions for other licensees extend to the
# entire whole, and thus to each and every part regardless of who wrote it.
#
# Thus, it is not the intent of this section to claim rights or contest
# your rights to work written entirely by you; rather, the intent is to
# exercise the right to control the distribution of derivative or
# collective works based on the Program.
#
# In addition, mere aggregation of another work not based on the Program
# with the Program (or with a work based on the Program) on a volume of
# a storage or distribution medium does not bring the other work under
# the scope of this License.
#
#   3. You may copy and distribute the Program (or a work based on it,
# under Section 2) in object code or executable form under the terms of
# Sections 1 and 2 above provided that you also do one of the following:
#
#     a) Accompany it with the complete corresponding machine-readable
#     source code, which must be distributed under the terms of Sections
#     1 and 2 above on a medium customarily used for software interchange; or,
#
#     b) Accompany it with a written offer, valid for at least three
#     years, to give any third party, for a charge no more than your
#     cost of physically performing source distribution, a complete
#     machine-readable copy of the corresponding source code, to be
#     distributed under the terms of Sections 1 and 2 above on a medium
#     customarily used for software interchange; or,
#
#     c) Accompany it with the information you received as to the offer
#     to distribute corresponding source code.  (This alternative is
#     allowed only for noncommercial distribution and only if you
#     received the program in object code or executable form with such
#     an offer, in accord with Subsection b above.)
#
# The source code for a work means the preferred form of the work for
# making modifications to it.  For an executable work, complete source
# code means all the source code for all modules it contains, plus any
# associated interface definition files, plus the scripts used to
# control compilation and installation of the executable.  However, as a
# special exception, the source code distributed need not include
# anything that is normally distributed (in either source or binary
# form) with the major components (compiler, kernel, and so on) of the
# operating system on which the executable runs, unless that component
# itself accompanies the executable.
#
# If distribution of executable or object code is made by offering
# access to copy from a designated place, then offering equivalent
# access to copy the source code from the same place counts as
# distribution of the source code, even though third parties are not
# compelled to copy the source along with the object code.
#
#   4. You may not copy, modify, sublicense, or distribute the Program
# except as expressly provided under this License.  Any attempt
# otherwise to copy, modify, sublicense or distribute the Program is
# void, and will automatically terminate your rights under this License.
# However, parties who have received copies, or rights, from you under
# this License will not have their licenses terminated so long as such
# parties remain in full compliance.
#
#   5. You are not required to accept this License, since you have not
# signed it.  However, nothing else grants you permission to modify or
# distribute the Program or its derivative works.  These actions are
# prohibited by law if you do not accept this License.  Therefore, by
# modifying or distributing the Program (or any work based on the
# Program), you indicate your acceptance of this License to do so, and
# all its terms and conditions for copying, distributing or modifying
# the Program or works based on it.
#
#   6. Each time you redistribute the Program (or any work based on the
# Program), the recipient automatically receives a license from the
# original licensor to copy, distribute or modify the Program subject to
# these terms and conditions.  You may not impose any further
# restrictions on the recipients' exercise of the rights granted herein.
# You are not responsible for enforcing compliance by third parties to
# this License.
#
#   7. If, as a consequence of a court judgment or allegation of patent
# infringement or for any other reason (not limited to patent issues),
# conditions are imposed on you (whether by court order, agreement or
# otherwise) that contradict the conditions of this License, they do not
# excuse you from the conditions of this License.  If you cannot
# distribute so as to satisfy simultaneously your obligations under this
# License and any other pertinent obligations, then as a consequence you
# may not distribute the Program at all.  For example, if a patent
# license would not permit royalty-free redistribution of the Program by
# all those who receive copies directly or indirectly through you, then
# the only way you could satisfy both it and this License would be to
# refrain entirely from distribution of the Program.
#
# If any portion of this section is held invalid or unenforceable under
# any particular circumstance, the balance of the section is intended to
# apply and the section as a whole is intended to apply in other
# circumstances.
#
# It is not the purpose of this section to induce you to infringe any
# patents or other property right claims or to contest validity of any
# such claims; this section has the sole purpose of protecting the
# integrity of the free software distribution system, which is
# implemented by public license practices.  Many people have made
# generous contributions to the wide range of software distributed
# through that system in reliance on consistent application of that
# system; it is up to the author/donor to decide if he or she is willing
# to distribute software through any other system and a licensee cannot
# impose that choice.
#
# This section is intended to make thoroughly clear what is believed to
# be a consequence of the rest of this License.
#
#   8. If the distribution and/or use of the Program is restricted in
# certain countries either by patents or by copyrighted interfaces, the
# original copyright holder who places the Program under this License
# may add an explicit geographical distribution limitation excluding
# those countries, so that distribution is permitted only in or among
# countries not thus excluded.  In such case, this License incorporates
# the limitation as if written in the body of this License.
#
#   9. The Free Software Foundation may publish revised and/or new versions
# of the General Public License from time to time.  Such new versions will
# be similar in spirit to the present version, but may differ in detail to
# address new problems or concerns.
#
# Each version is given a distinguishing version number.  If the Program
# specifies a version number of this License which applies to it and "any
# later version", you have the option of following the terms and conditions
# either of that version or of any later version published by the Free
# Software Foundation.  If the Program does not specify a version number of
# this License, you may choose any version ever published by the Free Software
# Foundation.
#
#   10. If you wish to incorporate parts of the Program into other free
# programs whose distribution conditions are different, write to the author
# to ask for permission.  For software which is copyrighted by the Free
# Software Foundation, write to the Free Software Foundation; we sometimes
# make exceptions for this.  Our decision will be guided by the two goals
# of preserving the free status of all derivatives of our free software and
# of promoting the sharing and reuse of software generally.
#
#           NO WARRANTY
#
#   11. BECAUSE THE PROGRAM IS LICENSED FREE OF CHARGE, THERE IS NO WARRANTY
# FOR THE PROGRAM, TO THE EXTENT PERMITTED BY APPLICABLE LAW.  EXCEPT WHEN
# OTHERWISE STATED IN WRITING THE COPYRIGHT HOLDERS AND/OR OTHER PARTIES
# PROVIDE THE PROGRAM "AS IS" WITHOUT WARRANTY OF ANY KIND, EITHER EXPRESSED
# OR IMPLIED, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF
# MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.  THE ENTIRE RISK AS
# TO THE QUALITY AND PERFORMANCE OF THE PROGRAM IS WITH YOU.  SHOULD THE
# PROGRAM PROVE DEFECTIVE, YOU ASSUME THE COST OF ALL NECESSARY SERVICING,
# REPAIR OR CORRECTION.
#
#   12. IN NO EVENT UNLESS REQUIRED BY APPLICABLE LAW OR AGREED TO IN WRITING
# WILL ANY COPYRIGHT HOLDER, OR ANY OTHER PARTY WHO MAY MODIFY AND/OR
# REDISTRIBUTE THE PROGRAM AS PERMITTED ABOVE, BE LIABLE TO YOU FOR DAMAGES,
# INCLUDING ANY GENERAL, SPECIAL, INCIDENTAL OR CONSEQUENTIAL DAMAGES ARISING
# OUT OF THE USE OR INABILITY TO USE THE PROGRAM (INCLUDING BUT NOT LIMITED
# TO LOSS OF DATA OR DATA BEING RENDERED INACCURATE OR LOSSES SUSTAINED BY
# YOU OR THIRD PARTIES OR A FAILURE OF THE PROGRAM TO OPERATE WITH ANY OTHER
# PROGRAMS), EVEN IF SUCH HOLDER OR OTHER PARTY HAS BEEN ADVISED OF THE
# POSSIBILITY OF SUCH DAMAGES.
#
#          END OF TERMS AND CONDITIONS
#
#       How to Apply These Terms to Your New Programs
#
#   If you develop a new program, and you want it to be of the greatest
# possible use to the public, the best way to achieve this is to make it
# free software which everyone can redistribute and change under these terms.
#
#   To do so, attach the following notices to the program.  It is safest
# to attach them to the start of each source file to most effectively
# convey the exclusion of warranty; and each file should have at least
# the "copyright" line and a pointer to where the full notice is found.
#
#     <one line to give the program's name and a brief idea of what it does.>
#     Copyright (C) <year>  <name of author>
#
#     This program is free software; you can redistribute it and/or modify
#     it under the terms of the GNU General Public License as published by
#     the Free Software Foundation; either version 2 of the License, or
#     (at your option) any later version.
#
#     This program is distributed in the hope that it will be useful,
#     but WITHOUT ANY WARRANTY; without even the implied warranty of
#     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#     GNU General Public License for more details.
#
#     You should have received a copy of the GNU General Public License
#     along with this program; if not, write to the Free Software
#     Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA
#
#
# Also add information on how to contact you by electronic and paper mail.
#
# If the program is interactive, make it output a short notice like this
# when it starts in an interactive mode:
#
#     Gnomovision version 69, Copyright (C) year name of author
#     Gnomovision comes with ABSOLUTELY NO WARRANTY; for details type `show w'.
#     This is free software, and you are welcome to redistribute it
#     under certain conditions; type `show c' for details.
#
# The hypothetical commands `show w' and `show c' should show the appropriate
# parts of the General Public License.  Of course, the commands you use may
# be called something other than `show w' and `show c'; they could even be
# mouse-clicks or menu items--whatever suits your program.
#
# You should also get your employer (if you work as a programmer) or your
# school, if any, to sign a "copyright disclaimer" for the program, if
# necessary.  Here is a sample; alter the names:
#
#   Yoyodyne, Inc., hereby disclaims all copyright interest in the program
#   `Gnomovision' (which makes passes at compilers) written by James Hacker.
#
#   <signature of Ty Coon>, 1 April 1989
#   Ty Coon, President of Vice
#
# This General Public License does not permit incorporating your program into
# proprietary programs.  If your program is a subroutine library, you may
# consider it more useful to permit linking proprietary applications with the
# library.  If this is what you want to do, use the GNU Lesser General
# Public License instead of this License.


include(CMakePushCheckState)
cmake_push_check_state()

# test for C++14 flags
if(NOT DISABLE_CXX_VERSION_CHECK)
  # try to use compiler flag -std=c++14
  include(TestCXXAcceptsFlag)
  check_cxx_accepts_flag("-std=c++14" CXX_FLAG_CXX14)

  include(CheckCXXSourceCompiles)
  cmake_push_check_state()
  set(CMAKE_REQUIRED_FLAGS "${CMAKE_REQUIRED_FLAGS} -std=c++14")
  check_cxx_source_compiles("
      #include <memory>

      int main() {
        std::make_unique<int>();
      }
    " CXX_LIB_SUPPORTS_CXX14)
  cmake_pop_check_state()
endif()

if(CXX_FLAG_CXX14 AND CXX_LIB_SUPPORTS_CXX14)
  set(CMAKE_REQUIRED_FLAGS "${CMAKE_REQUIRED_FLAGS} -std=c++14")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++14 ")
  set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -std=c++14 ")
  set(CMAKE_CXX_FLAGS_MINSIZEREL "${CMAKE_CXX_FLAGS_MINSIZEREL} -std=c++14 ")
  set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -std=c++14 ")
  set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "${CMAKE_CXX_FLAGS_RELWITHDEBINFO} -std=c++14 ")
  set(CXX_STD14_FLAGS "-std=c++14")
else()
  if(NOT DISABLE_CXX_VERSION_CHECK)
    # try to use compiler flag -std=c++1y for older compilers
    check_cxx_accepts_flag("-std=c++1y" CXX_FLAG_CXX1Y)

    include(CheckCXXSourceCompiles)
    cmake_push_check_state()
    set(CMAKE_REQUIRED_FLAGS "${CMAKE_REQUIRED_FLAGS} -std=c++1y")
    check_cxx_source_compiles("
        #include <memory>

        int main() {
          std::make_unique<int>();
        }
      " CXX_LIB_SUPPORTS_CXX1Y)
    cmake_pop_check_state()
  endif()
  if(CXX_FLAG_CXX1Y AND CXX_LIB_SUPPORTS_CXX1Y)
    set(CMAKE_REQUIRED_FLAGS "${CMAKE_REQUIRED_FLAGS} -std=c++1y" )
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++1y ")
    set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -std=c++1y ")
    set(CMAKE_CXX_FLAGS_MINSIZEREL "${CMAKE_CXX_FLAGS_MINSIZEREL} -std=c++1y ")
    set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -std=c++1y ")
    set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "${CMAKE_CXX_FLAGS_RELWITHDEBINFO} -std=c++1y ")
  set(CXX_STD14_FLAGS "-std=c++1y")
  endif()
endif()

# test for C++11 flags
if(NOT DISABLE_CXX_VERSION_CHECK
   AND NOT ((CXX_FLAG_CXX14 AND CXX_LIB_SUPPORTS_CXX14)
            OR (CXX_FLAG_CXX1Y AND CXX_LIB_SUPPORTS_CXX1Y)))
  # try to use compiler flag -std=c++11
  check_cxx_accepts_flag("-std=c++11" CXX_FLAG_CXX11)

  if(CXX_FLAG_CXX11)
    set(CMAKE_REQUIRED_FLAGS "${CMAKE_REQUIRED_FLAGS} -std=c++11")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 ")
    set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -std=c++11 ")
    set(CMAKE_CXX_FLAGS_MINSIZEREL "${CMAKE_CXX_FLAGS_MINSIZEREL} -std=c++11 ")
    set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -std=c++11 ")
    set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "${CMAKE_CXX_FLAGS_RELWITHDEBINFO} -std=c++11 ")
    set(CXX_STD11_FLAGS "-std=c++11")
  else()
    message(FATAL_ERROR "Your compiler does not seem to support C++11. If it does, please add any required flags to your CMAKE_CXX_FLAGS or set DISABLE_CXX_VERSION_CHECK.")
  endif()
endif()

# perform tests
include(CheckCXXSourceCompiles)

# nullptr
check_cxx_source_compiles("
    int main(void)
    {
      char* ch = nullptr;
      return 0;
    }
"  HAVE_NULLPTR
  )

if(HAVE_NULLPTR)
  check_cxx_source_compiles("
    #include <cstddef>

    int main(void)
    {
      std::nullptr_t npt = nullptr;
      return 0;
    }
"  HAVE_NULLPTR_T
    )
  if (NOT HAVE_NULLPTR_T)
    if(${CMAKE_CXX_COMPILER_ID} STREQUAL "Intel")
      message(FATAL_ERROR "Your compiler supports the 'nullptr' keyword, but not the type std::nullptr_t.
You are using an Intel compiler, where this error is typically caused by an outdated underlying system GCC.
Check if 'gcc --version' gives you at least GCC 4.6, otherwise please install a newer GCC and point your Intel compiler at it using the '-gcc-name' or '-gxx-name' options (see 'man icc' for further details)."
        )
    else()
      message(FATAL_ERROR "Your compiler supports the 'nullptr' keyword, but not the type std::nullptr_t.
THIS SHOULD NOT HAPPEN FOR YOUR COMPILER!
Please submit a bug at https://dune-project.org/flyspray/"
        )
    endif()
  endif()
endif()

# __attribute__((unused))
check_cxx_source_compiles("
   int main(void)
   {
     int __attribute__((unused)) foo;
     return 0;
   };
"  HAS_ATTRIBUTE_UNUSED
)

# __attribute__((deprecated))
check_cxx_source_compiles("
#define DEP __attribute__((deprecated))
   class bar
   {
     bar() DEP;
   };

   class peng { } DEP;

   template <class T>
   class t_bar
   {
     t_bar() DEP;
   };

   template <class T>
   class t_peng {
     t_peng() {};
   } DEP;

   void foo() DEP;

   void foo() {}

   int main(void)
   {
     return 0;
   };
"  HAS_ATTRIBUTE_DEPRECATED
)

# __attribute__((deprecated("msg")))
check_cxx_source_compiles("
#define DEP __attribute__((deprecated(\"message\")))
   class bar {
     bar() DEP;
   };

   class peng { } DEP;

   template <class T>
   class t_bar
   {
     t_bar() DEP;
   };

   template <class T>
   class t_peng
   {
     t_peng() {};
   } DEP;

   void foo() DEP;

   void foo() {}

   int main(void)
   {
     return 0;
   };
"  HAS_ATTRIBUTE_DEPRECATED_MSG
)

# constexpr
check_cxx_source_compiles("
  constexpr int foo()
  { return 0; }

  template<int v>
  struct A
  {
    static const int value = v;
  };

  int main(void)
  {
    return A<foo()>::value;
  }
" HAVE_CONSTEXPR
)

# keyword final
check_cxx_source_compiles("
  struct Foo
  {
    virtual void foo() final;
  };

  int main(void)
  {
    return 0;
  }
" HAVE_KEYWORD_FINAL
)

# range-based for
check_cxx_source_compiles("
  int main(void)
  {
    int arr[3];
    for(int &val : arr)
      val = 0;
  }
" HAVE_RANGE_BASED_FOR
)

# nonexcept specifier
check_cxx_source_compiles("
  void func1() noexcept {}

  void func2() noexcept(true) {}

  template <class T>
  void func3() noexcept(noexcept(T())) {}

  int main(void)
  {
    func1();
    func2();
    func3<int>();
  }
" HAVE_NOEXCEPT_SPECIFIER
)

# std::declval()
check_cxx_source_compiles("
  #include <utility>

  template<typename T>
  struct check;

  template<>
  struct check<int&&>
  {
    int pass() { return 0; }
  };

  int main(void)
  {
    return check<decltype(std::declval<int>())>().pass();
  }
" HAVE_STD_DECLVAL
  )

# full support for is_indexable (checking whether a type supports operator[])
check_cxx_source_compiles("
  #include <utility>
  #include <type_traits>
  #include <array>

  template <class T>
  typename std::add_rvalue_reference<T>::type declval();

  namespace detail {

    template<typename T, typename I, typename = int>
    struct _is_indexable
      : public std::false_type
    {};

    template<typename T, typename I>
    struct _is_indexable<T,I,typename std::enable_if<(sizeof(declval<T>()[declval<I>()]) > 0),int>::type>
      : public std::true_type
    {};

  }

  template<typename T, typename I = std::size_t>
  struct is_indexable
    : public detail::_is_indexable<T,I>
  {};

  struct foo_type {};

  int main()
  {
    double x;
    std::array<double,4> y;
    double z[5];
    foo_type f;

    static_assert(not is_indexable<decltype(x)>::value,\"scalar type\");
    static_assert(is_indexable<decltype(y)>::value,\"indexable class\");
    static_assert(is_indexable<decltype(z)>::value,\"array\");
    static_assert(not is_indexable<decltype(f)>::value,\"not indexable class\");
    static_assert(not is_indexable<decltype(y),foo_type>::value,\"custom index type\");

    return 0;
  }
" HAVE_IS_INDEXABLE_SUPPORT
  )


cmake_pop_check_state()

# find the threading library
# Use a copy FindThreads from CMake 3.1 due to its support of pthread
if(NOT DEFINED THREADS_PREFER_PTHREAD_FLAG)
  set(THREADS_PREFER_PTHREAD_FLAG 1)
endif()
if(${CMAKE_VERSION} VERSION_LESS "3.1")
  find_package(ThreadsCMake31)
else()
  find_package(Threads)
endif()

# see whether threading needs -no-as-needed
if(EXISTS /etc/dpkg/origins/ubuntu)
  set(NO_AS_NEEDED "-Wl,-no-as-needed ")
else(EXISTS /etc/dpkg/origins/ubuntu)
  set(NO_AS_NEEDED "")
endif(EXISTS /etc/dpkg/origins/ubuntu)

set(STDTHREAD_LINK_FLAGS "${NO_AS_NEEDED}${CMAKE_THREAD_LIBS_INIT}"
    CACHE STRING "Linker flags needed to get working C++11 threads support.  On Ubuntu it may be necessary to include -Wl,-no-as-needed (see FS#1650).")

# set linker flags
#
# in all implementations I know it is sufficient to set the linker flags when
# linking the final executable, so this should work.  In cmake, this appears
# to only work when building the project however, not for later config tests
# (contrary to CMAKE_CXX_FLAGS).  Luckily, later tests don't seem to use any
# threading...  (except for our own sanity check)
if(NOT STDTHREAD_LINK_FLAGS STREQUAL "")
  #set(vars CMAKE_EXE_LINKER_FLAGS ${CMAKE_CONFIGURATION_TYPES})
  # CMAKE_CONFIGURATION_TYPES seems to be empty.  Use the configurations from
  # adding -std=c++11 above instead.
  set(vars CMAKE_EXE_LINKER_FLAGS DEBUG MINSIZEREL RELEASE RELWITHDEBINFO)
  string(REPLACE ";" ";CMAKE_EXE_LINKER_FLAGS_" vars "${vars}")
  string(TOUPPER "${vars}" vars)
  foreach(var ${vars})
    if(NOT var STREQUAL "")
      set(${var} "${${var}} ${STDTHREAD_LINK_FLAGS}")
    endif()
  endforeach(var ${vars})
endif(NOT STDTHREAD_LINK_FLAGS STREQUAL "")

include(CheckCXXSourceRuns)
# check that the found configuration works
if(CMAKE_CROSSCOMPILING)
  message(WARNING "Crosscompiling, cannot run test program to see whether "
    "std::thread works.  Assuming that the found configuration does indeed "
    "work.")
endif(CMAKE_CROSSCOMPILING)

if(NOT DEFINED STDTHREAD_WORKS)
  if(NOT CMAKE_CROSSCOMPILING)
    # The value is not in the cache, so run check
    cmake_push_check_state()
    # tests seem to ignore CMAKE_EXE_LINKER_FLAGS
    set(CMAKE_REQUIRED_LIBRARIES "${STDTHREAD_LINK_FLAGS} ${CMAKE_REQUIRED_LIBRARIES}")
    check_cxx_source_runs("
      #include <thread>

      void dummy() {}

      int main() {
        std::thread t(dummy);
        t.join();
      }
    " STDTHREAD_WORKS)
    cmake_pop_check_state()
  endif(NOT CMAKE_CROSSCOMPILING)
  # put the found value into the cache.  Put it there even if we're
  # cross-compiling, so the user can find it.  Use FORCE:
  # check_cxx_source_runs() already puts the value in the cache but without
  # documentation; also the "if(NOT DEFINED STDTHREAD_WORKS)" will prevent us
  # from overwriting a value set by the user.
  set(STDTHREAD_WORKS "${STDTHREAD_WORKS}"
      CACHE BOOL "Whether std::thread works." FORCE)
endif(NOT DEFINED STDTHREAD_WORKS)

if(NOT STDTHREAD_WORKS)
  # Working C++11 threading support is required for dune.  In particular to
  # make things like lazyly initialized caches thread safe
  # (e.g. QuadratureRules::rule(), which needs std::call_once()).  If we don't
  # include the correct options during linking, there will be very funny
  # errors at runtime, ranging from segfaults to
  #
  #  terminate called after throwing an instance of 'std::system_error'
  #    what():  Unknown error 18446744073709551615
  message(FATAL_ERROR "Your system does not seem to have a working "
    "implementation of std::thread.  If it does, please set the linker flags "
    "required to get std::thread working in the cache variable "
    "STDTHREAD_LINK_FLAGS.  If you think this test is wrong, set the cache "
    "variable STDTHREAD_WORKS.")
endif(NOT STDTHREAD_WORKS)

cmake_pop_check_state()
