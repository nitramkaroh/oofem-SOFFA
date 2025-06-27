#pragma once

namespace FTensor
{
////////////////////////// A(i, j, k, l, m, )*B(i, ) //////////////////////////
// Base class:
template <class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, char i, char j, char k, char l, char m>
class Tensor5_times_Tensor1_single1
{
  Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> iterA;
  Tensor1_Expr<B, U, Dim0, i> iterB;

  template <int Current_Dim>
  typename promote<T, U>::V
  eval( const int N1, const int N2, const int N3, const int N4, const Number<Current_Dim> & ) const
  {
    return iterA( Current_Dim - 1, N1, N2, N3, N4 ) * iterB( Current_Dim - 1 ) + eval( N1, N2, N3, N4, Number<Current_Dim - 1>() );
  }
  typename promote<T, U>::V
  eval( const int N1, const int N2, const int N3, const int N4, const Number<1> & ) const
  {
    return iterA( 0, N1, N2, N3, N4 ) * iterB( 0 );
  }

  public:
  Tensor5_times_Tensor1_single1(
      const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
      const Tensor1_Expr<B, U, Dim0, i> &b ) :
      iterA( a ), iterB( b )
  {
  }

  typename promote<T, U>::V operator( const int N1, const int N2, const int N3, const int N4 ) const
  {
    return eval( N1, N2, N3, N4, Number<Dim0>() );
  }
};

// operator for A*b:
template <class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, char i, char j, char k, char l, char m>
Tensor4_Expr<Tensor5_times_Tensor1_single1<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m>, typename promote<T, U>::V, Dim1, Dim2, Dim3, Dim4, j, k, l, m>
operator*( const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
    const Tensor1_Expr<B, U, Dim0, i> &b )
{
  using TensorExpr
      = Tensor5_times_Tensor1_single1<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m>;
  return Tensor4_Expr<TensorExpr, typename promote<T, U>::V, Dim1, Dim2, Dim3, Dim4, j, k, l, m>( TensorExpr( a, b ) );
}

// operator for b*A:
template <class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, char i, char j, char k, char l, char m>
Tensor4_Expr<Tensor5_times_Tensor1_single1<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m>, typename promote<T, U>::V, Dim1, Dim2, Dim3, Dim4, j, k, l, m>
operator*( const Tensor1_Expr<B, U, Dim0, i> &b,
    const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a )
{
  using TensorExpr
      = Tensor5_times_Tensor1_single1<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m>;
  return Tensor4_Expr<TensorExpr, typename promote<T, U>::V, Dim1, Dim2, Dim3, Dim4, j, k, l, m>( TensorExpr( a, b ) );
}

////////////////////////// A(i, j, k, l, m, )*B(j, ) //////////////////////////
// Base class:
template <class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, char i, char j, char k, char l, char m>
class Tensor5_times_Tensor1_single2
{
  Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> iterA;
  Tensor1_Expr<B, U, Dim1, j> iterB;

  template <int Current_Dim>
  typename promote<T, U>::V
  eval( const int N1, const int N2, const int N3, const int N4, const Number<Current_Dim> & ) const
  {
    return iterA( N1, Current_Dim - 1, N2, N3, N4 ) * iterB( Current_Dim - 1 ) + eval( N1, N2, N3, N4, Number<Current_Dim - 1>() );
  }
  typename promote<T, U>::V
  eval( const int N1, const int N2, const int N3, const int N4, const Number<1> & ) const
  {
    return iterA( N1, 0, N2, N3, N4 ) * iterB( 0 );
  }

  public:
  Tensor5_times_Tensor1_single2(
      const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
      const Tensor1_Expr<B, U, Dim1, j> &b ) :
      iterA( a ), iterB( b )
  {
  }

  typename promote<T, U>::V operator( const int N1, const int N2, const int N3, const int N4 ) const
  {
    return eval( N1, N2, N3, N4, Number<Dim1>() );
  }
};

// operator for A*b:
template <class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, char i, char j, char k, char l, char m>
Tensor4_Expr<Tensor5_times_Tensor1_single2<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m>, typename promote<T, U>::V, Dim0, Dim2, Dim3, Dim4, i, k, l, m>
operator*( const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
    const Tensor1_Expr<B, U, Dim1, j> &b )
{
  using TensorExpr
      = Tensor5_times_Tensor1_single2<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m>;
  return Tensor4_Expr<TensorExpr, typename promote<T, U>::V, Dim0, Dim2, Dim3, Dim4, i, k, l, m>( TensorExpr( a, b ) );
}

// operator for b*A:
template <class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, char i, char j, char k, char l, char m>
Tensor4_Expr<Tensor5_times_Tensor1_single2<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m>, typename promote<T, U>::V, Dim0, Dim2, Dim3, Dim4, i, k, l, m>
operator*( const Tensor1_Expr<B, U, Dim1, j> &b,
    const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a )
{
  using TensorExpr
      = Tensor5_times_Tensor1_single2<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m>;
  return Tensor4_Expr<TensorExpr, typename promote<T, U>::V, Dim0, Dim2, Dim3, Dim4, i, k, l, m>( TensorExpr( a, b ) );
}

////////////////////////// A(i, j, k, l, m, )*B(k, ) //////////////////////////
// Base class:
template <class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, char i, char j, char k, char l, char m>
class Tensor5_times_Tensor1_single3
{
  Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> iterA;
  Tensor1_Expr<B, U, Dim2, k> iterB;

  template <int Current_Dim>
  typename promote<T, U>::V
  eval( const int N1, const int N2, const int N3, const int N4, const Number<Current_Dim> & ) const
  {
    return iterA( N1, N2, Current_Dim - 1, N3, N4 ) * iterB( Current_Dim - 1 ) + eval( N1, N2, N3, N4, Number<Current_Dim - 1>() );
  }
  typename promote<T, U>::V
  eval( const int N1, const int N2, const int N3, const int N4, const Number<1> & ) const
  {
    return iterA( N1, N2, 0, N3, N4 ) * iterB( 0 );
  }

  public:
  Tensor5_times_Tensor1_single3(
      const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
      const Tensor1_Expr<B, U, Dim2, k> &b ) :
      iterA( a ), iterB( b )
  {
  }

  typename promote<T, U>::V operator( const int N1, const int N2, const int N3, const int N4 ) const
  {
    return eval( N1, N2, N3, N4, Number<Dim2>() );
  }
};

// operator for A*b:
template <class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, char i, char j, char k, char l, char m>
Tensor4_Expr<Tensor5_times_Tensor1_single3<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m>, typename promote<T, U>::V, Dim0, Dim1, Dim3, Dim4, i, j, l, m>
operator*( const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
    const Tensor1_Expr<B, U, Dim2, k> &b )
{
  using TensorExpr
      = Tensor5_times_Tensor1_single3<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m>;
  return Tensor4_Expr<TensorExpr, typename promote<T, U>::V, Dim0, Dim1, Dim3, Dim4, i, j, l, m>( TensorExpr( a, b ) );
}

// operator for b*A:
template <class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, char i, char j, char k, char l, char m>
Tensor4_Expr<Tensor5_times_Tensor1_single3<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m>, typename promote<T, U>::V, Dim0, Dim1, Dim3, Dim4, i, j, l, m>
operator*( const Tensor1_Expr<B, U, Dim2, k> &b,
    const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a )
{
  using TensorExpr
      = Tensor5_times_Tensor1_single3<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m>;
  return Tensor4_Expr<TensorExpr, typename promote<T, U>::V, Dim0, Dim1, Dim3, Dim4, i, j, l, m>( TensorExpr( a, b ) );
}

////////////////////////// A(i, j, k, l, m, )*B(l, ) //////////////////////////
// Base class:
template <class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, char i, char j, char k, char l, char m>
class Tensor5_times_Tensor1_single4
{
  Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> iterA;
  Tensor1_Expr<B, U, Dim3, l> iterB;

  template <int Current_Dim>
  typename promote<T, U>::V
  eval( const int N1, const int N2, const int N3, const int N4, const Number<Current_Dim> & ) const
  {
    return iterA( N1, N2, N3, Current_Dim - 1, N4 ) * iterB( Current_Dim - 1 ) + eval( N1, N2, N3, N4, Number<Current_Dim - 1>() );
  }
  typename promote<T, U>::V
  eval( const int N1, const int N2, const int N3, const int N4, const Number<1> & ) const
  {
    return iterA( N1, N2, N3, 0, N4 ) * iterB( 0 );
  }

  public:
  Tensor5_times_Tensor1_single4(
      const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
      const Tensor1_Expr<B, U, Dim3, l> &b ) :
      iterA( a ), iterB( b )
  {
  }

  typename promote<T, U>::V operator( const int N1, const int N2, const int N3, const int N4 ) const
  {
    return eval( N1, N2, N3, N4, Number<Dim3>() );
  }
};

// operator for A*b:
template <class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, char i, char j, char k, char l, char m>
Tensor4_Expr<Tensor5_times_Tensor1_single4<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m>, typename promote<T, U>::V, Dim0, Dim1, Dim2, Dim4, i, j, k, m>
operator*( const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
    const Tensor1_Expr<B, U, Dim3, l> &b )
{
  using TensorExpr
      = Tensor5_times_Tensor1_single4<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m>;
  return Tensor4_Expr<TensorExpr, typename promote<T, U>::V, Dim0, Dim1, Dim2, Dim4, i, j, k, m>( TensorExpr( a, b ) );
}

// operator for b*A:
template <class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, char i, char j, char k, char l, char m>
Tensor4_Expr<Tensor5_times_Tensor1_single4<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m>, typename promote<T, U>::V, Dim0, Dim1, Dim2, Dim4, i, j, k, m>
operator*( const Tensor1_Expr<B, U, Dim3, l> &b,
    const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a )
{
  using TensorExpr
      = Tensor5_times_Tensor1_single4<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m>;
  return Tensor4_Expr<TensorExpr, typename promote<T, U>::V, Dim0, Dim1, Dim2, Dim4, i, j, k, m>( TensorExpr( a, b ) );
}

////////////////////////// A(i, j, k, l, m, )*B(m, ) //////////////////////////
// Base class:
template <class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, char i, char j, char k, char l, char m>
class Tensor5_times_Tensor1_single5
{
  Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> iterA;
  Tensor1_Expr<B, U, Dim4, m> iterB;

  template <int Current_Dim>
  typename promote<T, U>::V
  eval( const int N1, const int N2, const int N3, const int N4, const Number<Current_Dim> & ) const
  {
    return iterA( N1, N2, N3, N4, Current_Dim - 1 ) * iterB( Current_Dim - 1 ) + eval( N1, N2, N3, N4, Number<Current_Dim - 1>() );
  }
  typename promote<T, U>::V
  eval( const int N1, const int N2, const int N3, const int N4, const Number<1> & ) const
  {
    return iterA( N1, N2, N3, N4, 0 ) * iterB( 0 );
  }

  public:
  Tensor5_times_Tensor1_single5(
      const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
      const Tensor1_Expr<B, U, Dim4, m> &b ) :
      iterA( a ), iterB( b )
  {
  }

  typename promote<T, U>::V operator( const int N1, const int N2, const int N3, const int N4 ) const
  {
    return eval( N1, N2, N3, N4, Number<Dim4>() );
  }
};

// operator for A*b:
template <class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, char i, char j, char k, char l, char m>
Tensor4_Expr<Tensor5_times_Tensor1_single5<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m>, typename promote<T, U>::V, Dim0, Dim1, Dim2, Dim3, i, j, k, l>
operator*( const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
    const Tensor1_Expr<B, U, Dim4, m> &b )
{
  using TensorExpr
      = Tensor5_times_Tensor1_single5<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m>;
  return Tensor4_Expr<TensorExpr, typename promote<T, U>::V, Dim0, Dim1, Dim2, Dim3, i, j, k, l>( TensorExpr( a, b ) );
}

// operator for b*A:
template <class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, char i, char j, char k, char l, char m>
Tensor4_Expr<Tensor5_times_Tensor1_single5<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m>, typename promote<T, U>::V, Dim0, Dim1, Dim2, Dim3, i, j, k, l>
operator*( const Tensor1_Expr<B, U, Dim4, m> &b,
    const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a )
{
  using TensorExpr
      = Tensor5_times_Tensor1_single5<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m>;
  return Tensor4_Expr<TensorExpr, typename promote<T, U>::V, Dim0, Dim1, Dim2, Dim3, i, j, k, l>( TensorExpr( a, b ) );
}



  
  
}
