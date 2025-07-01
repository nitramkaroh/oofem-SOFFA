/* A helper class that automatically returns the internal contraction
   of a rank 5 tensor. */

#pragma once

namespace FTensor
{

  template <class A, class T, int Dim> class Tensor5_contracted_34
{
  A iterA;

  public:
  T operator()( const int N0, const int N1, const int N2 ) const
  {
    T result( 0 );
    for ( int ii = 0; ii < Dim; ++ii ) {
      result += iterA( N0, N1, N2, ii, ii );
    }
    return result;
  }
  Tensor5_contracted_34( const A &a ) : iterA( a ) {}
};

  template <class A, class T, int Dim> class Tensor5_contracted_24
{
  A iterA;

  public:
  T operator()( const int N0, const int N1, const int N2 ) const
  {
    T result( 0 );
    for ( int ii = 0; ii < Dim; ++ii ) {
      result += iterA( N0, N1, ii, N2, ii );
    }
    return result;
  }
  Tensor5_contracted_24( const A &a ) : iterA( a ) {}
};


  template <class A, class T, int Dim> class Tensor5_contracted_23
  {
    A iterA;

  public:
    T operator()(const int N0, const int N1, const int N2) const
    {
      T result(0);
      for(int ii = 0; ii < Dim; ++ii)
        {
          result += iterA(N0, N1, ii, ii, N2);
        }
      return result;
    }
    Tensor5_contracted_23(const A &a) : iterA(a) {}
  };

    template <class A, class T, int Dim> class Tensor5_contracted_14
  {
    A iterA;

    public:
    T operator()( const int N0, const int N1, const int N2 ) const
    {
      T result( 0 );
      for ( int ii = 0; ii < Dim; ++ii ) {
        result += iterA( N0, ii, N1, N2, ii );
      }
      return result;
    }
    Tensor5_contracted_14( const A &a ) : iterA( a ) {}
  };

  template <class A, class T, int Dim> class Tensor5_contracted_13
  {
    A iterA;

  public:
    T operator()(const int N0, const int N1, const int N2) const
    {
      T result(0);
      for(int ii = 0; ii < Dim; ++ii)
        {
          result += iterA(N0, ii, N1, ii, N2);
        }
      return result;
    }
    Tensor5_contracted_13(const A &a) : iterA(a) {}
  };

  template <class A, class T, int Dim> class Tensor5_contracted_12
  {
    A iterA;

  public:
    T operator()(const int N0, const int N1, const int N2) const
    {
      T result(0);
      for(int ii = 0; ii < Dim; ++ii)
        {
          result += iterA(N0, ii, ii, N1, N2);
        }
      return result;
    }
    Tensor5_contracted_12(const A &a) : iterA(a) {}
  };

  template <class A, class T, int Dim> class Tensor5_contracted_04
  {
    A iterA;

  public:
    T operator()( const int N0, const int N1, const int N2 ) const
    {
      T result( 0 );
      for ( int ii = 0; ii < Dim; ++ii ) {
        result += iterA( ii, N0, N1, N2, ii );
      }
      return result;
    }
    Tensor5_contracted_04( const A &a ) : iterA( a ) {}
  };


  template <class A, class T, int Dim> class Tensor5_contracted_03
  {
    A iterA;

  public:
    T operator()(const int N0, const int N1, const int N2) const
    {
      T result(0);
      for(int ii = 0; ii < Dim; ++ii)
        {
          result += iterA(ii, N0, N1, ii, N2);
        }
      return result;
    }
    Tensor5_contracted_03(const A &a) : iterA(a) {}
  };

  template <class A, class T, int Dim> class Tensor5_contracted_02
  {
    A iterA;

  public:
    T operator()(const int N0, const int N1, const int N2) const
    {
      T result(0);
      for(int ii = 0; ii < Dim; ++ii)
        {
          result += iterA(ii, N0, ii, N1, N2);
        }
      return result;
    }
    Tensor5_contracted_02(const A &a) : iterA(a) {}
  };

  template <class A, class T, int Dim> class Tensor5_contracted_01
  {
    A iterA;

  public:
    T operator()(const int N0, const int N1, const int N2) const
    {
      T result(0);
      for(int ii = 0; ii < Dim; ++ii)
        {
          result += iterA(ii, ii, N0, N1, N2);
        }
      return result;
    }
    Tensor5_contracted_01(const A &a) : iterA(a) {}
  };
}
