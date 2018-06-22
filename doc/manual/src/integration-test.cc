// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
//
// Copyright (c) 2015, Carsten Gräser
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//
// 1. Redistributions of source code must retain the above copyright notice,
//    this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright notice,
//    this list of conditions and the following disclaimer in the documentation
//    and/or other materials provided with the distribution.
//
// 3. Neither the name of the copyright holder nor the names of its
//    contributors may be used to endorse or promote products derived from this
//    software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
// LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
// CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
// SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
// INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
// CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
// ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
// POSSIBILITY OF SUCH DAMAGE.

#include <memory>
#include <iostream>
#include <iomanip>
#include <array>
#include <chrono>
#include <functional>


template<class K, std::size_t N>
class StaticVector : public std::array<K, N>
{
public:

  StaticVector()
  {
    auto& self = *this;
    for (std::size_t i=0; i<N; ++i)
      self[i] = 0;
  }

  StaticVector(const K& k)
  {
    auto& self = *this;
    for (std::size_t i=0; i<N; ++i)
      self[i] = k;
  }

  void axpy(const K& k, const StaticVector& other)
  {
    auto& self = *this;
    for (std::size_t i=0; i<N; ++i)
      self[i] += k*other[i];
  }

  void operator*=(const K& k)
  {
    auto& self = *this;
    for (std::size_t i=0; i<N; ++i)
      self[i] *= k;
  }
};

template <class K, std::size_t N>
inline std::ostream& operator<< (std::ostream& s, const StaticVector<K,N>& v)
{
  s << "[";
  for (std::size_t i=0; i<N-1; i++)
    s << v[i] << ",";
    s << v[N-1] << "]";
  return s;
}

template<class Domain, class Range, class F>
Range integrateOverInterval(const F& f, Domain a, Domain b, int n)
{
  Domain h = (b-a)/n;
  Range result = 0.0;

  std::array<std::pair<double,double>,1> quadRule = {{ {0.5, 1} }};

  Range localResult = 0.0;
  for (std::size_t i=0; i<n; ++i)
  {
    for (auto&& quadPoint : quadRule)
//      localResult.axpy(quadPoint.second, f(a + (i+quadPoint.first)*h));
      result.axpy(quadPoint.second, f(a + (i+quadPoint.first)*h));
//    result.axpy(h, localResult);
  }
  result *= h;
  return result;
}



template<class Domain, class Range, class F>
Range integrateOverInterval_Evaluate(const F& f, Domain a, Domain b, int n)
{
  Domain h = (b-a)/n;
  Range result = 0.0;

  std::array<std::pair<double,double>,1> quadRule = {{ {0.5, 1} }};

  Range temp;

  Range localResult = 0.0;
  for (std::size_t i=0; i<n; ++i)
  {
    for (auto&& quadPoint : quadRule)
    {
      f.evaluate(a + (i+quadPoint.first)*h, temp);
//      localResult.axpy(quadPoint.second, temp);
      result.axpy(quadPoint.second, temp);
    }
//    result.axpy(h, localResult);
  }
  result *= h;
  return result;
}




class Timer
{
  using TimeStamp = decltype(std::chrono::high_resolution_clock::now());

  TimeStamp start_;
public:

  Timer () :
    start_(std::chrono::high_resolution_clock::now())
  {}

  double msElapsed() const
  {
    std::chrono::duration<double, std::milli> elapsed
      = std::chrono::high_resolution_clock::now()-start_;
    return elapsed.count();
  }
};



template<class D, class R>
class VirtualFunction
{
public:
  using Domain = D;
  using Range = R;
  virtual void evaluate(const Domain&, Range&) const = 0;
};

template<class D, class R>
class Function
{
public:
  using Domain = D;
  using Range = R;
};






template<class Base>
class F : public Base
{
public:
  using Range = typename Base::Range;
  using Domain = typename Base::Domain;

  void imp(const Domain& x, Range& y) const
  {
    for(std::size_t i=0; i<y.size(); ++i)
    {
//      y[i] = x*x;
      y[i] = x+i;
//      y[i] = x*x + i;
    }
  }
  void evaluate(const Domain& x, Range& y) const
  {
    imp(x,y);
  }

  Range operator()(const Domain& x) const
  {
    Range y;
    imp(x, y);
    return y;
  }
};



template<std::size_t blockSize>
void testWithVectorSize(std::size_t n)
{
  double a = 0;
  double b = 1;

  using Range = StaticVector<double,blockSize>;
  using Domain = double;

  std::array<double, 4> timing;
  std::array<Range, 4> result;

  std::size_t N = n/blockSize;

  // For each test do one warump run with N-1 before.
  {
    using Interface = Function<Domain, Range>;
    auto f = F<Interface>();
    for (int k=N-1; k<=N; ++k)
    {
      Timer t;
      result[0] = integrateOverInterval<Domain, Range>(f, a, b, k);
      timing[0] = t.msElapsed();
    }
  }

  {
    using Interface = Function<Domain, Range>;
    auto f = F<Interface>();
    for (int k=N-1; k<=N; ++k)
    {
      Timer t;
      result[1] = integrateOverInterval_Evaluate<Domain, Range>(f, a, b, k);
      timing[1] = t.msElapsed();
    }
  }


  {
    using Interface = Function<Domain, Range>;
    auto f = F<Interface>();
    for (int k=N-1; k<=N; ++k)
    {
      Timer t;
      result[2] = integrateOverInterval<Domain, Range>(std::function<Range(Domain)>(f), a, b, k);
      timing[2] = t.msElapsed();
    }
  }

  {
    using Interface = VirtualFunction<Domain, Range>;
    std::shared_ptr<Interface> f = std::make_shared<F<Interface>>();
    Interface& fBase = *f;
    for (int k=N-1; k<=N; ++k)
    {
      Timer t;
      result[3] = integrateOverInterval_Evaluate<Domain, Range, Interface>(fBase, a, b, k);
      timing[3] = t.msElapsed();
    }
  }

  std::cout << std::setw(3) << blockSize << "   ";
  for(int i=0; i<4; ++i)
    std::cout << std::setw(14) << (int)(timing[i]);
  for(int i=1; i<4; ++i)
    if (result[0]!=result[i])
      std::cout << " error";
  std::cout << std::endl;
}


int main(int argc, char** argv)
{

  std::size_t n = 100000000;

  std::cout << "size      operator()    evaluate() std::function  v evaluate()" << std::endl;
  testWithVectorSize<1>(n);
  testWithVectorSize<2>(n);
  testWithVectorSize<3>(n);
  testWithVectorSize<4>(n);
  testWithVectorSize<5>(n);
  testWithVectorSize<6>(n);
  testWithVectorSize<7>(n);
  testWithVectorSize<8>(n);
  testWithVectorSize<9>(n);
  testWithVectorSize<10>(n);
  testWithVectorSize<11>(n);
  testWithVectorSize<12>(n);
  testWithVectorSize<13>(n);
  testWithVectorSize<14>(n);
  testWithVectorSize<15>(n);
  testWithVectorSize<16>(n);

  return 0;
}
