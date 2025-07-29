#pragma once

namespace FTensor
{
////////////////////////// A(i, j, k, l, m, n, )*B(i, ) //////////////////////////
// Base class:
template <class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
class Tensor6_times_Tensor1_single1
{
  Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> iterA;
  Tensor1_Expr<B, U, Dim0, i> iterB;

  template <int Current_Dim>
  typename promote<T, U>::V
  eval( const int N1, const int N2, const int N3, const int N4, const int N5, const Number<Current_Dim> & ) const
  {
    return iterA( Current_Dim - 1, N1, N2, N3, N4, N5 ) * iterB( Current_Dim - 1 ) + eval( N1, N2, N3, N4, N5, Number<Current_Dim - 1>() );
  }
  typename promote<T, U>::V
  eval( const int N1, const int N2, const int N3, const int N4, const int N5, const Number<1> & ) const
  {
    return iterA( 0, N1, N2, N3, N4, N5 ) * iterB( 0 );
  }

  public:
  Tensor6_times_Tensor1_single1(
      const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
      const Tensor1_Expr<B, U, Dim0, i> &b ) :
      iterA( a ), iterB( b )
  {
  }

  typename promote<T, U>::V operator()( const int N1, const int N2, const int N3, const int N4, const int N5 ) const
  {
    return eval( N1, N2, N3, N4, N5, Number<Dim0>() );
  }
};

// operator for A*b:
template <class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
Tensor5_Expr<Tensor6_times_Tensor1_single1<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n>, typename promote<T, U>::V, Dim1, Dim2, Dim3, Dim4, Dim5, j, k, l, m, n>
operator*( const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
    const Tensor1_Expr<B, U, Dim0, i> &b )
{
  using TensorExpr
      = Tensor6_times_Tensor1_single1<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n>;
  return Tensor5_Expr<TensorExpr, typename promote<T, U>::V, Dim1, Dim2, Dim3, Dim4, Dim5, j, k, l, m, n>( TensorExpr( a, b ) );
}

// operator for b*A:
template <class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
Tensor5_Expr<Tensor6_times_Tensor1_single1<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n>, typename promote<T, U>::V, Dim1, Dim2, Dim3, Dim4, Dim5, j, k, l, m, n>
operator*( const Tensor1_Expr<B, U, Dim0, i> &b,
    const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a )
{
  using TensorExpr
      = Tensor6_times_Tensor1_single1<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n>;
  return Tensor5_Expr<TensorExpr, typename promote<T, U>::V, Dim1, Dim2, Dim3, Dim4, Dim5, j, k, l, m, n>( TensorExpr( a, b ) );
}

////////////////////////// A(i, j, k, l, m, n, )*B(j, ) //////////////////////////
// Base class:
template <class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
class Tensor6_times_Tensor1_single2
{
  Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> iterA;
  Tensor1_Expr<B, U, Dim1, j> iterB;

  template <int Current_Dim>
  typename promote<T, U>::V
  eval( const int N1, const int N2, const int N3, const int N4, const int N5, const Number<Current_Dim> & ) const
  {
    return iterA( N1, Current_Dim - 1, N2, N3, N4, N5 ) * iterB( Current_Dim - 1 ) + eval( N1, N2, N3, N4, N5, Number<Current_Dim - 1>() );
  }
  typename promote<T, U>::V
  eval( const int N1, const int N2, const int N3, const int N4, const int N5, const Number<1> & ) const
  {
    return iterA( N1, 0, N2, N3, N4, N5 ) * iterB( 0 );
  }

  public:
  Tensor6_times_Tensor1_single2(
      const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
      const Tensor1_Expr<B, U, Dim1, j> &b ) :
      iterA( a ), iterB( b )
  {
  }

  typename promote<T, U>::V operator()( const int N1, const int N2, const int N3, const int N4, const int N5 ) const
  {
    return eval( N1, N2, N3, N4, N5, Number<Dim1>() );
  }
};

// operator for A*b:
template <class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
Tensor5_Expr<Tensor6_times_Tensor1_single2<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n>, typename promote<T, U>::V, Dim0, Dim2, Dim3, Dim4, Dim5, i, k, l, m, n>
operator*( const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
    const Tensor1_Expr<B, U, Dim1, j> &b )
{
  using TensorExpr
      = Tensor6_times_Tensor1_single2<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n>;
  return Tensor5_Expr<TensorExpr, typename promote<T, U>::V, Dim0, Dim2, Dim3, Dim4, Dim5, i, k, l, m, n>( TensorExpr( a, b ) );
}

// operator for b*A:
template <class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
Tensor5_Expr<Tensor6_times_Tensor1_single2<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n>, typename promote<T, U>::V, Dim0, Dim2, Dim3, Dim4, Dim5, i, k, l, m, n>
operator*( const Tensor1_Expr<B, U, Dim1, j> &b,
    const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a )
{
  using TensorExpr
      = Tensor6_times_Tensor1_single2<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n>;
  return Tensor5_Expr<TensorExpr, typename promote<T, U>::V, Dim0, Dim2, Dim3, Dim4, Dim5, i, k, l, m, n>( TensorExpr( a, b ) );
}

////////////////////////// A(i, j, k, l, m, n, )*B(k, ) //////////////////////////
// Base class:
template <class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
class Tensor6_times_Tensor1_single3
{
  Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> iterA;
  Tensor1_Expr<B, U, Dim2, k> iterB;

  template <int Current_Dim>
  typename promote<T, U>::V
  eval( const int N1, const int N2, const int N3, const int N4, const int N5, const Number<Current_Dim> & ) const
  {
    return iterA( N1, N2, Current_Dim - 1, N3, N4, N5 ) * iterB( Current_Dim - 1 ) + eval( N1, N2, N3, N4, N5, Number<Current_Dim - 1>() );
  }
  typename promote<T, U>::V
  eval( const int N1, const int N2, const int N3, const int N4, const int N5, const Number<1> & ) const
  {
    return iterA( N1, N2, 0, N3, N4, N5 ) * iterB( 0 );
  }

  public:
  Tensor6_times_Tensor1_single3(
      const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
      const Tensor1_Expr<B, U, Dim2, k> &b ) :
      iterA( a ), iterB( b )
  {
  }

  typename promote<T, U>::V operator()( const int N1, const int N2, const int N3, const int N4, const int N5 ) const
  {
    return eval( N1, N2, N3, N4, N5, Number<Dim2>() );
  }
};

// operator for A*b:
template <class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
Tensor5_Expr<Tensor6_times_Tensor1_single3<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n>, typename promote<T, U>::V, Dim0, Dim1, Dim3, Dim4, Dim5, i, j, l, m, n>
operator*( const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
    const Tensor1_Expr<B, U, Dim2, k> &b )
{
  using TensorExpr
      = Tensor6_times_Tensor1_single3<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n>;
  return Tensor5_Expr<TensorExpr, typename promote<T, U>::V, Dim0, Dim1, Dim3, Dim4, Dim5, i, j, l, m, n>( TensorExpr( a, b ) );
}

// operator for b*A:
template <class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
Tensor5_Expr<Tensor6_times_Tensor1_single3<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n>, typename promote<T, U>::V, Dim0, Dim1, Dim3, Dim4, Dim5, i, j, l, m, n>
operator*( const Tensor1_Expr<B, U, Dim2, k> &b,
    const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a )
{
  using TensorExpr
      = Tensor6_times_Tensor1_single3<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n>;
  return Tensor5_Expr<TensorExpr, typename promote<T, U>::V, Dim0, Dim1, Dim3, Dim4, Dim5, i, j, l, m, n>( TensorExpr( a, b ) );
}

////////////////////////// A(i, j, k, l, m, n, )*B(l, ) //////////////////////////
// Base class:
template <class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
class Tensor6_times_Tensor1_single4
{
  Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> iterA;
  Tensor1_Expr<B, U, Dim3, l> iterB;

  template <int Current_Dim>
  typename promote<T, U>::V
  eval( const int N1, const int N2, const int N3, const int N4, const int N5, const Number<Current_Dim> & ) const
  {
    return iterA( N1, N2, N3, Current_Dim - 1, N4, N5 ) * iterB( Current_Dim - 1 ) + eval( N1, N2, N3, N4, N5, Number<Current_Dim - 1>() );
  }
  typename promote<T, U>::V
  eval( const int N1, const int N2, const int N3, const int N4, const int N5, const Number<1> & ) const
  {
    return iterA( N1, N2, N3, 0, N4, N5 ) * iterB( 0 );
  }

  public:
  Tensor6_times_Tensor1_single4(
      const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
      const Tensor1_Expr<B, U, Dim3, l> &b ) :
      iterA( a ), iterB( b )
  {
  }

  typename promote<T, U>::V operator()( const int N1, const int N2, const int N3, const int N4, const int N5 ) const
  {
    return eval( N1, N2, N3, N4, N5, Number<Dim3>() );
  }
};

// operator for A*b:
template <class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
Tensor5_Expr<Tensor6_times_Tensor1_single4<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n>, typename promote<T, U>::V, Dim0, Dim1, Dim2, Dim4, Dim5, i, j, k, m, n>
operator*( const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
    const Tensor1_Expr<B, U, Dim3, l> &b )
{
  using TensorExpr
      = Tensor6_times_Tensor1_single4<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n>;
  return Tensor5_Expr<TensorExpr, typename promote<T, U>::V, Dim0, Dim1, Dim2, Dim4, Dim5, i, j, k, m, n>( TensorExpr( a, b ) );
}

// operator for b*A:
template <class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
Tensor5_Expr<Tensor6_times_Tensor1_single4<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n>, typename promote<T, U>::V, Dim0, Dim1, Dim2, Dim4, Dim5, i, j, k, m, n>
operator*( const Tensor1_Expr<B, U, Dim3, l> &b,
    const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a )
{
  using TensorExpr
      = Tensor6_times_Tensor1_single4<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n>;
  return Tensor5_Expr<TensorExpr, typename promote<T, U>::V, Dim0, Dim1, Dim2, Dim4, Dim5, i, j, k, m, n>( TensorExpr( a, b ) );
}

////////////////////////// A(i, j, k, l, m, n, )*B(m, ) //////////////////////////
// Base class:
template <class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
class Tensor6_times_Tensor1_single5
{
  Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> iterA;
  Tensor1_Expr<B, U, Dim4, m> iterB;

  template <int Current_Dim>
  typename promote<T, U>::V
  eval( const int N1, const int N2, const int N3, const int N4, const int N5, const Number<Current_Dim> & ) const
  {
    return iterA( N1, N2, N3, N4, Current_Dim - 1, N5 ) * iterB( Current_Dim - 1 ) + eval( N1, N2, N3, N4, N5, Number<Current_Dim - 1>() );
  }
  typename promote<T, U>::V
  eval( const int N1, const int N2, const int N3, const int N4, const int N5, const Number<1> & ) const
  {
    return iterA( N1, N2, N3, N4, 0, N5 ) * iterB( 0 );
  }

  public:
  Tensor6_times_Tensor1_single5(
      const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
      const Tensor1_Expr<B, U, Dim4, m> &b ) :
      iterA( a ), iterB( b )
  {
  }

  typename promote<T, U>::V operator()( const int N1, const int N2, const int N3, const int N4, const int N5 ) const
  {
    return eval( N1, N2, N3, N4, N5, Number<Dim4>() );
  }
};

// operator for A*b:
template <class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
Tensor5_Expr<Tensor6_times_Tensor1_single5<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n>, typename promote<T, U>::V, Dim0, Dim1, Dim2, Dim3, Dim5, i, j, k, l, n>
operator*( const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
    const Tensor1_Expr<B, U, Dim4, m> &b )
{
  using TensorExpr
      = Tensor6_times_Tensor1_single5<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n>;
  return Tensor5_Expr<TensorExpr, typename promote<T, U>::V, Dim0, Dim1, Dim2, Dim3, Dim5, i, j, k, l, n>( TensorExpr( a, b ) );
}

// operator for b*A:
template <class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
Tensor5_Expr<Tensor6_times_Tensor1_single5<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n>, typename promote<T, U>::V, Dim0, Dim1, Dim2, Dim3, Dim5, i, j, k, l, n>
operator*( const Tensor1_Expr<B, U, Dim4, m> &b,
    const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a )
{
  using TensorExpr
      = Tensor6_times_Tensor1_single5<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n>;
  return Tensor5_Expr<TensorExpr, typename promote<T, U>::V, Dim0, Dim1, Dim2, Dim3, Dim5, i, j, k, l, n>( TensorExpr( a, b ) );
}

////////////////////////// A(i, j, k, l, m, n, )*B(n, ) //////////////////////////
// Base class:
template <class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
class Tensor6_times_Tensor1_single6
{
  Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> iterA;
  Tensor1_Expr<B, U, Dim5, n> iterB;

  template <int Current_Dim>
  typename promote<T, U>::V
  eval( const int N1, const int N2, const int N3, const int N4, const int N5, const Number<Current_Dim> & ) const
  {
    return iterA( N1, N2, N3, N4, N5, Current_Dim - 1 ) * iterB( Current_Dim - 1 ) + eval( N1, N2, N3, N4, N5, Number<Current_Dim - 1>() );
  }
  typename promote<T, U>::V
  eval( const int N1, const int N2, const int N3, const int N4, const int N5, const Number<1> & ) const
  {
    return iterA( N1, N2, N3, N4, N5, 0 ) * iterB( 0 );
  }

  public:
  Tensor6_times_Tensor1_single6(
      const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
      const Tensor1_Expr<B, U, Dim5, n> &b ) :
      iterA( a ), iterB( b )
  {
  }

  typename promote<T, U>::V operator()( const int N1, const int N2, const int N3, const int N4, const int N5 ) const
  {
    return eval( N1, N2, N3, N4, N5, Number<Dim5>() );
  }
};

// operator for A*b:
template <class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
Tensor5_Expr<Tensor6_times_Tensor1_single6<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n>, typename promote<T, U>::V, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m>
operator*( const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
    const Tensor1_Expr<B, U, Dim5, n> &b )
{
  using TensorExpr
      = Tensor6_times_Tensor1_single6<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n>;
  return Tensor5_Expr<TensorExpr, typename promote<T, U>::V, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m>( TensorExpr( a, b ) );
}

// operator for b*A:
template <class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
Tensor5_Expr<Tensor6_times_Tensor1_single6<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n>, typename promote<T, U>::V, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m>
operator*( const Tensor1_Expr<B, U, Dim5, n> &b,
    const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a )
{
  using TensorExpr
      = Tensor6_times_Tensor1_single6<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n>;
  return Tensor5_Expr<TensorExpr, typename promote<T, U>::V, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m>( TensorExpr( a, b ) );
}


} //namespace FTensor
