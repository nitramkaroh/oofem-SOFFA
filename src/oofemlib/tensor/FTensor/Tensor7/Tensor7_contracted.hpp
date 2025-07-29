/* A helper class that automatically returns the internal contraction
   of a rank 7 tensor. */

#pragma once

namespace FTensor {

template <class A, class T, int Dim>
class Tensor7_contracted_56
{
  A iterA;

  public:
  T operator()( const int N0, const int N1, const int N2, const int N3, const int N4 ) const
  {
    T result( 0 );
    for ( int ii = 0; ii < Dim; ++ii ) {
      result += iterA( N0, N1, N2, N3, N4, ii, ii );
    }
    return result;
  }
  Tensor7_contracted_56( const A &a ) :
      iterA( a ) {}
};

template <class A, class T, int Dim>
class Tensor7_contracted_46
{
  A iterA;

  public:
  T operator()( const int N0, const int N1, const int N2, const int N3, const int N4 ) const
  {
    T result( 0 );
    for ( int ii = 0; ii < Dim; ++ii ) {
      result += iterA( N0, N1, N2, N3, ii, N4, ii );
    }
    return result;
  }
  Tensor7_contracted_46( const A &a ) :
      iterA( a ) {}
};

template <class A, class T, int Dim>
class Tensor7_contracted_45
{
  A iterA;

  public:
  T operator()( const int N0, const int N1, const int N2, const int N3, const int N4 ) const
  {
    T result( 0 );
    for ( int ii = 0; ii < Dim; ++ii ) {
      result += iterA( N0, N1, N2, N3, ii, ii, N4 );
    }
    return result;
  }
  Tensor7_contracted_45( const A &a ) :
      iterA( a ) {}
};

template <class A, class T, int Dim>
class Tensor7_contracted_36
{
  A iterA;

  public:
  T operator()( const int N0, const int N1, const int N2, const int N3, const int N4 ) const
  {
    T result( 0 );
    for ( int ii = 0; ii < Dim; ++ii ) {
      result += iterA( N0, N1, N2, ii, N3, N4, ii );
    }
    return result;
  }
  Tensor7_contracted_36( const A &a ) :
      iterA( a ) {}
};

template <class A, class T, int Dim>
class Tensor7_contracted_35
{
  A iterA;

  public:
  T operator()( const int N0, const int N1, const int N2, const int N3, const int N4 ) const
  {
    T result( 0 );
    for ( int ii = 0; ii < Dim; ++ii ) {
      result += iterA( N0, N1, N2, ii, N3, ii, N4 );
    }
    return result;
  }
  Tensor7_contracted_35( const A &a ) :
      iterA( a ) {}
};

template <class A, class T, int Dim>
class Tensor7_contracted_34
{
  A iterA;

  public:
  T operator()( const int N0, const int N1, const int N2, const int N3, const int N4 ) const
  {
    T result( 0 );
    for ( int ii = 0; ii < Dim; ++ii ) {
      result += iterA( N0, N1, N2, ii, ii, N3, N4 );
    }
    return result;
  }
  Tensor7_contracted_34( const A &a ) :
      iterA( a ) {}
};

template <class A, class T, int Dim>
class Tensor7_contracted_26
{
  A iterA;

  public:
  T operator()( const int N0, const int N1, const int N2, const int N3, const int N4 ) const
  {
    T result( 0 );
    for ( int ii = 0; ii < Dim; ++ii ) {
      result += iterA( N0, N1, ii, N2, N3, N4, ii );
    }
    return result;
  }
  Tensor7_contracted_26( const A &a ) :
      iterA( a ) {}
};

template <class A, class T, int Dim>
class Tensor7_contracted_25
{
  A iterA;

  public:
  T operator()( const int N0, const int N1, const int N2, const int N3, const int N4 ) const
  {
    T result( 0 );
    for ( int ii = 0; ii < Dim; ++ii ) {
      result += iterA( N0, N1, ii, N2, N3, ii, N4 );
    }
    return result;
  }
  Tensor7_contracted_25( const A &a ) :
      iterA( a ) {}
};

template <class A, class T, int Dim>
class Tensor7_contracted_24
{
  A iterA;

  public:
  T operator()( const int N0, const int N1, const int N2, const int N3, const int N4 ) const
  {
    T result( 0 );
    for ( int ii = 0; ii < Dim; ++ii ) {
      result += iterA( N0, N1, ii, N2, ii, N3, N4 );
    }
    return result;
  }
  Tensor7_contracted_24( const A &a ) :
      iterA( a ) {}
};


template <class A, class T, int Dim>
class Tensor7_contracted_23
{
  A iterA;

  public:
  T operator()( const int N0, const int N1, const int N2, const int N3, const int N4 ) const
  {
    T result( 0 );
    for ( int ii = 0; ii < Dim; ++ii ) {
      result += iterA( N0, N1, ii, ii, N2, N3, N4 );
    }
    return result;
  }
  Tensor7_contracted_23( const A &a ) :
      iterA( a ) {}
};

template <class A, class T, int Dim>
class Tensor7_contracted_16
{
  A iterA;

  public:
  T operator()( const int N0, const int N1, const int N2, const int N3, const int N4 ) const
  {
    T result( 0 );
    for ( int ii = 0; ii < Dim; ++ii ) {
      result += iterA( N0, ii, N1, N2, N3, N4, ii);
    }
    return result;
  }
  Tensor7_contracted_16( const A &a ) :
      iterA( a ) {}
};

template <class A, class T, int Dim>
class Tensor7_contracted_15
{
  A iterA;

  public:
  T operator()( const int N0, const int N1, const int N2, const int N3, const int N4 ) const
  {
    T result( 0 );
    for ( int ii = 0; ii < Dim; ++ii ) {
      result += iterA( N0, ii, N1, N2, N3, ii, N4 );
    }
    return result;
  }
  Tensor7_contracted_15( const A &a ) :
      iterA( a ) {}
};

template <class A, class T, int Dim>
class Tensor7_contracted_14
{
  A iterA;

  public:
  T operator()( const int N0, const int N1, const int N2, const int N3, const int N4 ) const
  {
    T result( 0 );
    for ( int ii = 0; ii < Dim; ++ii ) {
      result += iterA( N0, ii, N1, N2, ii, N3, N4 );
    }
    return result;
  }
  Tensor7_contracted_14( const A &a ) :
      iterA( a ) {}
};

template <class A, class T, int Dim>
class Tensor7_contracted_13
{
  A iterA;

  public:
  T operator()( const int N0, const int N1, const int N2, const int N3, const int N4 ) const
  {
    T result( 0 );
    for ( int ii = 0; ii < Dim; ++ii ) {
      result += iterA( N0, ii, N1, ii, N2, N3, N4 );
    }
    return result;
  }
  Tensor7_contracted_13( const A &a ) :
      iterA( a ) {}
};

template <class A, class T, int Dim>
class Tensor7_contracted_12
{
  A iterA;

  public:
  T operator()( const int N0, const int N1, const int N2, const int N3, const int N4 ) const
  {
    T result( 0 );
    for ( int ii = 0; ii < Dim; ++ii ) {
      result += iterA( N0, ii, ii, N1, N2, N3, N4 );
    }
    return result;
  }
  Tensor7_contracted_12( const A &a ) :
      iterA( a ) {}
};

template <class A, class T, int Dim>
class Tensor7_contracted_06
{
  A iterA;

  public:
  T operator()( const int N0, const int N1, const int N2, const int N3, const int N4 ) const
  {
    T result( 0 );
    for ( int ii = 0; ii < Dim; ++ii ) {
      result += iterA( ii, N0, N1, N2, N3, N4, ii );
    }
    return result;
  }
  Tensor7_contracted_06( const A &a ) :
      iterA( a ) {}
};

template <class A, class T, int Dim>
class Tensor7_contracted_05
{

  A iterA;

  public:
  T operator()( const int N0, const int N1, const int N2, const int N3, const int N4 ) const
  {
    T result( 0 );
    for ( int ii = 0; ii < Dim; ++ii ) {
      result += iterA( ii, N0, N1, N2, N3, ii, N4 );
    }
    return result;
  }
  Tensor7_contracted_05( const A &a ) :
      iterA( a ) {}
};

template <class A, class T, int Dim>
class Tensor7_contracted_04
{
  A iterA;

  public:
  T operator()( const int N0, const int N1, const int N2, const int N3, const int N4 ) const
  {
    T result( 0 );
    for ( int ii = 0; ii < Dim; ++ii ) {
      result += iterA( ii, N0, N1, N2, ii, N3, N4 );
    }
    return result;
  }
  Tensor7_contracted_04( const A &a ) :
      iterA( a ) {}
};


template <class A, class T, int Dim>
class Tensor7_contracted_03
{
  A iterA;

  public:
  T operator()( const int N0, const int N1, const int N2, const int N3, const int N4 ) const
  {
    T result( 0 );
    for ( int ii = 0; ii < Dim; ++ii ) {
      result += iterA( ii, N0, N1, ii, N2, N3, N4 );
    }
    return result;
  }
  Tensor7_contracted_03( const A &a ) :
      iterA( a ) {}
};

template <class A, class T, int Dim>
class Tensor7_contracted_02
{
  A iterA;

  public:
  T operator()( const int N0, const int N1, const int N2, const int N3, const int N4 ) const
  {
    T result( 0 );
    for ( int ii = 0; ii < Dim; ++ii ) {
      result += iterA( ii, N0, ii, N1, N2, N3, N4 );
    }
    return result;
  }
  Tensor7_contracted_02( const A &a ) :
      iterA( a ) {}
};

template <class A, class T, int Dim>
class Tensor7_contracted_01
{
  A iterA;

  public:
  T operator()( const int N0, const int N1, const int N2, const int N3, const int N4 ) const
  {
    T result( 0 );
    for ( int ii = 0; ii < Dim; ++ii ) {
      result += iterA( ii, ii, N0, N1, N2, N3, N4 );
    }
    return result;
  }
  Tensor7_contracted_01( const A &a ) :
      iterA( a ) {}
};
} // namespace FTensor
