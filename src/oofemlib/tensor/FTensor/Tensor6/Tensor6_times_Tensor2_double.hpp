#pragma once

namespace FTensor
{
// Universal A(i, j, k, l, m, n, )*B(o, p, ) double contraction
template <class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p, int DimA, int DimB, int DimC, int DimD, int DimX, int DimY, char a, char b, char c, char d, char x, char y>
class Tensor6_times_Tensor2_double
{
  Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> iterA;
  Tensor2_Expr<B, U, Dim6, Dim7, o, p> iterB;

  public:
  typename promote<T, U>::V operator()( const int N1, const int N2, const int N3, const int N4 ) const
  {
    typename promote<T, U>::V result( 0 );
    for ( int xx = 0; xx < DimX, ++xx )
      for ( int yy = 0; yy < DimY, ++yy ) {
        // Permutation is where the indices get checked.
        result += Permutation6<DimA, DimB, DimC, DimD, DimX, DimY, a, b, c, d, x, y>().eval(
                      iterA, N1, N2, N3, N4, xx, yy )
            * Permutation2<DimX, DimY, x, y>().eval( iterB, xx, yy );
      }
    return result;
  }

  Tensor6_times_Tensor2_double(
      const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &iter_a,
      const Tensor2_Expr<B, U, Dim6, Dim7, o, p> &iter_b ) :
      iterA( iter_a ), iterB( iter_b )
  {
  }
};

// A(i, j, k, l, m, n, )*B(j, i, )
template <class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*( const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
    const Tensor2_Expr<B, U, Dim1, Dim0, j, i> &b )
{
  using TensorExpr
      = Tensor6_times_Tensor2_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim1, Dim0, i, j, k, l, m, n, j, i, Dim2, Dim3, Dim4, Dim5, Dim1, Dim0, k, l, m, n, j, i>;
  return Tensor4_Expr<TensorExpr, typename promote<T, U>::V, Dim2, Dim3, Dim4, Dim5, k, l, m, n>( TensorExpr( a, b ) );
}

// B(j, i, )*A(i, j, k, l, m, n, )
template <class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, i, j, k, l, m, n>
auto operator*( const Tensor2_Expr<B, U, Dim1, Dim0, j, i> &b,
    const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a )
{
  return a * b;
}

// A(i, j, k, l, m, n, )*B(i, j, )
template <class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*( const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
    const Tensor2_Expr<B, U, Dim0, Dim1, i, j> &b )
{
  using TensorExpr
      = Tensor6_times_Tensor2_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim0, Dim1, i, j, k, l, m, n, i, j, Dim2, Dim3, Dim4, Dim5, Dim0, Dim1, k, l, m, n, i, j>;
  return Tensor4_Expr<TensorExpr, typename promote<T, U>::V, Dim2, Dim3, Dim4, Dim5, k, l, m, n>( TensorExpr( a, b ) );
}

// B(i, j, )*A(i, j, k, l, m, n, )
template <class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, i, j, k, l, m, n>
auto operator*( const Tensor2_Expr<B, U, Dim0, Dim1, i, j> &b,
    const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a )
{
  return a * b;
}

// A(i, j, k, l, m, n, )*B(k, i, )
template <class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*( const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
    const Tensor2_Expr<B, U, Dim2, Dim0, k, i> &b )
{
  using TensorExpr
      = Tensor6_times_Tensor2_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim2, Dim0, i, j, k, l, m, n, k, i, Dim1, Dim3, Dim4, Dim5, Dim2, Dim0, j, l, m, n, k, i>;
  return Tensor4_Expr<TensorExpr, typename promote<T, U>::V, Dim1, Dim3, Dim4, Dim5, j, l, m, n>( TensorExpr( a, b ) );
}

// B(k, i, )*A(i, j, k, l, m, n, )
template <class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, i, j, k, l, m, n>
auto operator*( const Tensor2_Expr<B, U, Dim2, Dim0, k, i> &b,
    const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a )
{
  return a * b;
}

// A(i, j, k, l, m, n, )*B(i, k, )
template <class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*( const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
    const Tensor2_Expr<B, U, Dim0, Dim2, i, k> &b )
{
  using TensorExpr
      = Tensor6_times_Tensor2_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim0, Dim2, i, j, k, l, m, n, i, k, Dim1, Dim3, Dim4, Dim5, Dim0, Dim2, j, l, m, n, i, k>;
  return Tensor4_Expr<TensorExpr, typename promote<T, U>::V, Dim1, Dim3, Dim4, Dim5, j, l, m, n>( TensorExpr( a, b ) );
}

// B(i, k, )*A(i, j, k, l, m, n, )
template <class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, i, j, k, l, m, n>
auto operator*( const Tensor2_Expr<B, U, Dim0, Dim2, i, k> &b,
    const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a )
{
  return a * b;
}

// A(i, j, k, l, m, n, )*B(l, i, )
template <class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*( const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
    const Tensor2_Expr<B, U, Dim3, Dim0, l, i> &b )
{
  using TensorExpr
      = Tensor6_times_Tensor2_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim3, Dim0, i, j, k, l, m, n, l, i, Dim1, Dim2, Dim4, Dim5, Dim3, Dim0, j, k, m, n, l, i>;
  return Tensor4_Expr<TensorExpr, typename promote<T, U>::V, Dim1, Dim2, Dim4, Dim5, j, k, m, n>( TensorExpr( a, b ) );
}

// B(l, i, )*A(i, j, k, l, m, n, )
template <class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, i, j, k, l, m, n>
auto operator*( const Tensor2_Expr<B, U, Dim3, Dim0, l, i> &b,
    const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a )
{
  return a * b;
}

// A(i, j, k, l, m, n, )*B(i, l, )
template <class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*( const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
    const Tensor2_Expr<B, U, Dim0, Dim3, i, l> &b )
{
  using TensorExpr
      = Tensor6_times_Tensor2_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim0, Dim3, i, j, k, l, m, n, i, l, Dim1, Dim2, Dim4, Dim5, Dim0, Dim3, j, k, m, n, i, l>;
  return Tensor4_Expr<TensorExpr, typename promote<T, U>::V, Dim1, Dim2, Dim4, Dim5, j, k, m, n>( TensorExpr( a, b ) );
}

// B(i, l, )*A(i, j, k, l, m, n, )
template <class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, i, j, k, l, m, n>
auto operator*( const Tensor2_Expr<B, U, Dim0, Dim3, i, l> &b,
    const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a )
{
  return a * b;
}

// A(i, j, k, l, m, n, )*B(m, i, )
template <class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*( const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
    const Tensor2_Expr<B, U, Dim4, Dim0, m, i> &b )
{
  using TensorExpr
      = Tensor6_times_Tensor2_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim4, Dim0, i, j, k, l, m, n, m, i, Dim1, Dim2, Dim3, Dim5, Dim4, Dim0, j, k, l, n, m, i>;
  return Tensor4_Expr<TensorExpr, typename promote<T, U>::V, Dim1, Dim2, Dim3, Dim5, j, k, l, n>( TensorExpr( a, b ) );
}

// B(m, i, )*A(i, j, k, l, m, n, )
template <class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, i, j, k, l, m, n>
auto operator*( const Tensor2_Expr<B, U, Dim4, Dim0, m, i> &b,
    const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a )
{
  return a * b;
}

// A(i, j, k, l, m, n, )*B(i, m, )
template <class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*( const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
    const Tensor2_Expr<B, U, Dim0, Dim4, i, m> &b )
{
  using TensorExpr
      = Tensor6_times_Tensor2_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim0, Dim4, i, j, k, l, m, n, i, m, Dim1, Dim2, Dim3, Dim5, Dim0, Dim4, j, k, l, n, i, m>;
  return Tensor4_Expr<TensorExpr, typename promote<T, U>::V, Dim1, Dim2, Dim3, Dim5, j, k, l, n>( TensorExpr( a, b ) );
}

// B(i, m, )*A(i, j, k, l, m, n, )
template <class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, i, j, k, l, m, n>
auto operator*( const Tensor2_Expr<B, U, Dim0, Dim4, i, m> &b,
    const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a )
{
  return a * b;
}

// A(i, j, k, l, m, n, )*B(n, i, )
template <class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*( const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
    const Tensor2_Expr<B, U, Dim5, Dim0, n, i> &b )
{
  using TensorExpr
      = Tensor6_times_Tensor2_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim5, Dim0, i, j, k, l, m, n, n, i, Dim1, Dim2, Dim3, Dim4, Dim5, Dim0, j, k, l, m, n, i>;
  return Tensor4_Expr<TensorExpr, typename promote<T, U>::V, Dim1, Dim2, Dim3, Dim4, j, k, l, m>( TensorExpr( a, b ) );
}

// B(n, i, )*A(i, j, k, l, m, n, )
template <class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, i, j, k, l, m, n>
auto operator*( const Tensor2_Expr<B, U, Dim5, Dim0, n, i> &b,
    const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a )
{
  return a * b;
}

// A(i, j, k, l, m, n, )*B(i, n, )
template <class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*( const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
    const Tensor2_Expr<B, U, Dim0, Dim5, i, n> &b )
{
  using TensorExpr
      = Tensor6_times_Tensor2_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim0, Dim5, i, j, k, l, m, n, i, n, Dim1, Dim2, Dim3, Dim4, Dim0, Dim5, j, k, l, m, i, n>;
  return Tensor4_Expr<TensorExpr, typename promote<T, U>::V, Dim1, Dim2, Dim3, Dim4, j, k, l, m>( TensorExpr( a, b ) );
}

// B(i, n, )*A(i, j, k, l, m, n, )
template <class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, i, j, k, l, m, n>
auto operator*( const Tensor2_Expr<B, U, Dim0, Dim5, i, n> &b,
    const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a )
{
  return a * b;
}

// A(i, j, k, l, m, n, )*B(k, j, )
template <class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*( const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
    const Tensor2_Expr<B, U, Dim2, Dim1, k, j> &b )
{
  using TensorExpr
      = Tensor6_times_Tensor2_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim2, Dim1, i, j, k, l, m, n, k, j, Dim0, Dim3, Dim4, Dim5, Dim2, Dim1, i, l, m, n, k, j>;
  return Tensor4_Expr<TensorExpr, typename promote<T, U>::V, Dim0, Dim3, Dim4, Dim5, i, l, m, n>( TensorExpr( a, b ) );
}

// B(k, j, )*A(i, j, k, l, m, n, )
template <class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, i, j, k, l, m, n>
auto operator*( const Tensor2_Expr<B, U, Dim2, Dim1, k, j> &b,
    const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a )
{
  return a * b;
}

// A(i, j, k, l, m, n, )*B(j, k, )
template <class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*( const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
    const Tensor2_Expr<B, U, Dim1, Dim2, j, k> &b )
{
  using TensorExpr
      = Tensor6_times_Tensor2_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim1, Dim2, i, j, k, l, m, n, j, k, Dim0, Dim3, Dim4, Dim5, Dim1, Dim2, i, l, m, n, j, k>;
  return Tensor4_Expr<TensorExpr, typename promote<T, U>::V, Dim0, Dim3, Dim4, Dim5, i, l, m, n>( TensorExpr( a, b ) );
}

// B(j, k, )*A(i, j, k, l, m, n, )
template <class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, i, j, k, l, m, n>
auto operator*( const Tensor2_Expr<B, U, Dim1, Dim2, j, k> &b,
    const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a )
{
  return a * b;
}

// A(i, j, k, l, m, n, )*B(l, j, )
template <class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*( const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
    const Tensor2_Expr<B, U, Dim3, Dim1, l, j> &b )
{
  using TensorExpr
      = Tensor6_times_Tensor2_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim3, Dim1, i, j, k, l, m, n, l, j, Dim0, Dim2, Dim4, Dim5, Dim3, Dim1, i, k, m, n, l, j>;
  return Tensor4_Expr<TensorExpr, typename promote<T, U>::V, Dim0, Dim2, Dim4, Dim5, i, k, m, n>( TensorExpr( a, b ) );
}

// B(l, j, )*A(i, j, k, l, m, n, )
template <class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, i, j, k, l, m, n>
auto operator*( const Tensor2_Expr<B, U, Dim3, Dim1, l, j> &b,
    const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a )
{
  return a * b;
}

// A(i, j, k, l, m, n, )*B(j, l, )
template <class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*( const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
    const Tensor2_Expr<B, U, Dim1, Dim3, j, l> &b )
{
  using TensorExpr
      = Tensor6_times_Tensor2_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim1, Dim3, i, j, k, l, m, n, j, l, Dim0, Dim2, Dim4, Dim5, Dim1, Dim3, i, k, m, n, j, l>;
  return Tensor4_Expr<TensorExpr, typename promote<T, U>::V, Dim0, Dim2, Dim4, Dim5, i, k, m, n>( TensorExpr( a, b ) );
}

// B(j, l, )*A(i, j, k, l, m, n, )
template <class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, i, j, k, l, m, n>
auto operator*( const Tensor2_Expr<B, U, Dim1, Dim3, j, l> &b,
    const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a )
{
  return a * b;
}

// A(i, j, k, l, m, n, )*B(m, j, )
template <class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*( const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
    const Tensor2_Expr<B, U, Dim4, Dim1, m, j> &b )
{
  using TensorExpr
      = Tensor6_times_Tensor2_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim4, Dim1, i, j, k, l, m, n, m, j, Dim0, Dim2, Dim3, Dim5, Dim4, Dim1, i, k, l, n, m, j>;
  return Tensor4_Expr<TensorExpr, typename promote<T, U>::V, Dim0, Dim2, Dim3, Dim5, i, k, l, n>( TensorExpr( a, b ) );
}

// B(m, j, )*A(i, j, k, l, m, n, )
template <class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, i, j, k, l, m, n>
auto operator*( const Tensor2_Expr<B, U, Dim4, Dim1, m, j> &b,
    const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a )
{
  return a * b;
}

// A(i, j, k, l, m, n, )*B(j, m, )
template <class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*( const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
    const Tensor2_Expr<B, U, Dim1, Dim4, j, m> &b )
{
  using TensorExpr
      = Tensor6_times_Tensor2_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim1, Dim4, i, j, k, l, m, n, j, m, Dim0, Dim2, Dim3, Dim5, Dim1, Dim4, i, k, l, n, j, m>;
  return Tensor4_Expr<TensorExpr, typename promote<T, U>::V, Dim0, Dim2, Dim3, Dim5, i, k, l, n>( TensorExpr( a, b ) );
}

// B(j, m, )*A(i, j, k, l, m, n, )
template <class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, i, j, k, l, m, n>
auto operator*( const Tensor2_Expr<B, U, Dim1, Dim4, j, m> &b,
    const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a )
{
  return a * b;
}

// A(i, j, k, l, m, n, )*B(n, j, )
template <class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*( const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
    const Tensor2_Expr<B, U, Dim5, Dim1, n, j> &b )
{
  using TensorExpr
      = Tensor6_times_Tensor2_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim5, Dim1, i, j, k, l, m, n, n, j, Dim0, Dim2, Dim3, Dim4, Dim5, Dim1, i, k, l, m, n, j>;
  return Tensor4_Expr<TensorExpr, typename promote<T, U>::V, Dim0, Dim2, Dim3, Dim4, i, k, l, m>( TensorExpr( a, b ) );
}

// B(n, j, )*A(i, j, k, l, m, n, )
template <class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, i, j, k, l, m, n>
auto operator*( const Tensor2_Expr<B, U, Dim5, Dim1, n, j> &b,
    const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a )
{
  return a * b;
}

// A(i, j, k, l, m, n, )*B(j, n, )
template <class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*( const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
    const Tensor2_Expr<B, U, Dim1, Dim5, j, n> &b )
{
  using TensorExpr
      = Tensor6_times_Tensor2_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim1, Dim5, i, j, k, l, m, n, j, n, Dim0, Dim2, Dim3, Dim4, Dim1, Dim5, i, k, l, m, j, n>;
  return Tensor4_Expr<TensorExpr, typename promote<T, U>::V, Dim0, Dim2, Dim3, Dim4, i, k, l, m>( TensorExpr( a, b ) );
}

// B(j, n, )*A(i, j, k, l, m, n, )
template <class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, i, j, k, l, m, n>
auto operator*( const Tensor2_Expr<B, U, Dim1, Dim5, j, n> &b,
    const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a )
{
  return a * b;
}

// A(i, j, k, l, m, n, )*B(l, k, )
template <class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*( const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
    const Tensor2_Expr<B, U, Dim3, Dim2, l, k> &b )
{
  using TensorExpr
      = Tensor6_times_Tensor2_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim3, Dim2, i, j, k, l, m, n, l, k, Dim0, Dim1, Dim4, Dim5, Dim3, Dim2, i, j, m, n, l, k>;
  return Tensor4_Expr<TensorExpr, typename promote<T, U>::V, Dim0, Dim1, Dim4, Dim5, i, j, m, n>( TensorExpr( a, b ) );
}

// B(l, k, )*A(i, j, k, l, m, n, )
template <class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, i, j, k, l, m, n>
auto operator*( const Tensor2_Expr<B, U, Dim3, Dim2, l, k> &b,
    const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a )
{
  return a * b;
}

// A(i, j, k, l, m, n, )*B(k, l, )
template <class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*( const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
    const Tensor2_Expr<B, U, Dim2, Dim3, k, l> &b )
{
  using TensorExpr
      = Tensor6_times_Tensor2_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim2, Dim3, i, j, k, l, m, n, k, l, Dim0, Dim1, Dim4, Dim5, Dim2, Dim3, i, j, m, n, k, l>;
  return Tensor4_Expr<TensorExpr, typename promote<T, U>::V, Dim0, Dim1, Dim4, Dim5, i, j, m, n>( TensorExpr( a, b ) );
}

// B(k, l, )*A(i, j, k, l, m, n, )
template <class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, i, j, k, l, m, n>
auto operator*( const Tensor2_Expr<B, U, Dim2, Dim3, k, l> &b,
    const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a )
{
  return a * b;
}

// A(i, j, k, l, m, n, )*B(m, k, )
template <class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*( const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
    const Tensor2_Expr<B, U, Dim4, Dim2, m, k> &b )
{
  using TensorExpr
      = Tensor6_times_Tensor2_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim4, Dim2, i, j, k, l, m, n, m, k, Dim0, Dim1, Dim3, Dim5, Dim4, Dim2, i, j, l, n, m, k>;
  return Tensor4_Expr<TensorExpr, typename promote<T, U>::V, Dim0, Dim1, Dim3, Dim5, i, j, l, n>( TensorExpr( a, b ) );
}

// B(m, k, )*A(i, j, k, l, m, n, )
template <class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, i, j, k, l, m, n>
auto operator*( const Tensor2_Expr<B, U, Dim4, Dim2, m, k> &b,
    const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a )
{
  return a * b;
}

// A(i, j, k, l, m, n, )*B(k, m, )
template <class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*( const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
    const Tensor2_Expr<B, U, Dim2, Dim4, k, m> &b )
{
  using TensorExpr
      = Tensor6_times_Tensor2_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim2, Dim4, i, j, k, l, m, n, k, m, Dim0, Dim1, Dim3, Dim5, Dim2, Dim4, i, j, l, n, k, m>;
  return Tensor4_Expr<TensorExpr, typename promote<T, U>::V, Dim0, Dim1, Dim3, Dim5, i, j, l, n>( TensorExpr( a, b ) );
}

// B(k, m, )*A(i, j, k, l, m, n, )
template <class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, i, j, k, l, m, n>
auto operator*( const Tensor2_Expr<B, U, Dim2, Dim4, k, m> &b,
    const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a )
{
  return a * b;
}

// A(i, j, k, l, m, n, )*B(n, k, )
template <class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*( const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
    const Tensor2_Expr<B, U, Dim5, Dim2, n, k> &b )
{
  using TensorExpr
      = Tensor6_times_Tensor2_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim5, Dim2, i, j, k, l, m, n, n, k, Dim0, Dim1, Dim3, Dim4, Dim5, Dim2, i, j, l, m, n, k>;
  return Tensor4_Expr<TensorExpr, typename promote<T, U>::V, Dim0, Dim1, Dim3, Dim4, i, j, l, m>( TensorExpr( a, b ) );
}

// B(n, k, )*A(i, j, k, l, m, n, )
template <class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, i, j, k, l, m, n>
auto operator*( const Tensor2_Expr<B, U, Dim5, Dim2, n, k> &b,
    const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a )
{
  return a * b;
}

// A(i, j, k, l, m, n, )*B(k, n, )
template <class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*( const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
    const Tensor2_Expr<B, U, Dim2, Dim5, k, n> &b )
{
  using TensorExpr
      = Tensor6_times_Tensor2_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim2, Dim5, i, j, k, l, m, n, k, n, Dim0, Dim1, Dim3, Dim4, Dim2, Dim5, i, j, l, m, k, n>;
  return Tensor4_Expr<TensorExpr, typename promote<T, U>::V, Dim0, Dim1, Dim3, Dim4, i, j, l, m>( TensorExpr( a, b ) );
}

// B(k, n, )*A(i, j, k, l, m, n, )
template <class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, i, j, k, l, m, n>
auto operator*( const Tensor2_Expr<B, U, Dim2, Dim5, k, n> &b,
    const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a )
{
  return a * b;
}

// A(i, j, k, l, m, n, )*B(m, l, )
template <class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*( const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
    const Tensor2_Expr<B, U, Dim4, Dim3, m, l> &b )
{
  using TensorExpr
      = Tensor6_times_Tensor2_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim4, Dim3, i, j, k, l, m, n, m, l, Dim0, Dim1, Dim2, Dim5, Dim4, Dim3, i, j, k, n, m, l>;
  return Tensor4_Expr<TensorExpr, typename promote<T, U>::V, Dim0, Dim1, Dim2, Dim5, i, j, k, n>( TensorExpr( a, b ) );
}

// B(m, l, )*A(i, j, k, l, m, n, )
template <class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, i, j, k, l, m, n>
auto operator*( const Tensor2_Expr<B, U, Dim4, Dim3, m, l> &b,
    const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a )
{
  return a * b;
}

// A(i, j, k, l, m, n, )*B(l, m, )
template <class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*( const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
    const Tensor2_Expr<B, U, Dim3, Dim4, l, m> &b )
{
  using TensorExpr
      = Tensor6_times_Tensor2_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim3, Dim4, i, j, k, l, m, n, l, m, Dim0, Dim1, Dim2, Dim5, Dim3, Dim4, i, j, k, n, l, m>;
  return Tensor4_Expr<TensorExpr, typename promote<T, U>::V, Dim0, Dim1, Dim2, Dim5, i, j, k, n>( TensorExpr( a, b ) );
}

// B(l, m, )*A(i, j, k, l, m, n, )
template <class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, i, j, k, l, m, n>
auto operator*( const Tensor2_Expr<B, U, Dim3, Dim4, l, m> &b,
    const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a )
{
  return a * b;
}

// A(i, j, k, l, m, n, )*B(n, l, )
template <class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*( const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
    const Tensor2_Expr<B, U, Dim5, Dim3, n, l> &b )
{
  using TensorExpr
      = Tensor6_times_Tensor2_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim5, Dim3, i, j, k, l, m, n, n, l, Dim0, Dim1, Dim2, Dim4, Dim5, Dim3, i, j, k, m, n, l>;
  return Tensor4_Expr<TensorExpr, typename promote<T, U>::V, Dim0, Dim1, Dim2, Dim4, i, j, k, m>( TensorExpr( a, b ) );
}

// B(n, l, )*A(i, j, k, l, m, n, )
template <class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, i, j, k, l, m, n>
auto operator*( const Tensor2_Expr<B, U, Dim5, Dim3, n, l> &b,
    const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a )
{
  return a * b;
}

// A(i, j, k, l, m, n, )*B(l, n, )
template <class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*( const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
    const Tensor2_Expr<B, U, Dim3, Dim5, l, n> &b )
{
  using TensorExpr
      = Tensor6_times_Tensor2_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim3, Dim5, i, j, k, l, m, n, l, n, Dim0, Dim1, Dim2, Dim4, Dim3, Dim5, i, j, k, m, l, n>;
  return Tensor4_Expr<TensorExpr, typename promote<T, U>::V, Dim0, Dim1, Dim2, Dim4, i, j, k, m>( TensorExpr( a, b ) );
}

// B(l, n, )*A(i, j, k, l, m, n, )
template <class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, i, j, k, l, m, n>
auto operator*( const Tensor2_Expr<B, U, Dim3, Dim5, l, n> &b,
    const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a )
{
  return a * b;
}

// A(i, j, k, l, m, n, )*B(n, m, )
template <class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*( const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
    const Tensor2_Expr<B, U, Dim5, Dim4, n, m> &b )
{
  using TensorExpr
      = Tensor6_times_Tensor2_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim5, Dim4, i, j, k, l, m, n, n, m, Dim0, Dim1, Dim2, Dim3, Dim5, Dim4, i, j, k, l, n, m>;
  return Tensor4_Expr<TensorExpr, typename promote<T, U>::V, Dim0, Dim1, Dim2, Dim3, i, j, k, l>( TensorExpr( a, b ) );
}

// B(n, m, )*A(i, j, k, l, m, n, )
template <class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, i, j, k, l, m, n>
auto operator*( const Tensor2_Expr<B, U, Dim5, Dim4, n, m> &b,
    const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a )
{
  return a * b;
}

// A(i, j, k, l, m, n, )*B(m, n, )
template <class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*( const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a,
    const Tensor2_Expr<B, U, Dim4, Dim5, m, n> &b )
{
  using TensorExpr
      = Tensor6_times_Tensor2_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim4, Dim5, i, j, k, l, m, n, m, n, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n>;
  return Tensor4_Expr<TensorExpr, typename promote<T, U>::V, Dim0, Dim1, Dim2, Dim3, i, j, k, l>( TensorExpr( a, b ) );
}

// B(m, n, )*A(i, j, k, l, m, n, )
template <class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, i, j, k, l, m, n>
auto operator*( const Tensor2_Expr<B, U, Dim4, Dim5, m, n> &b,
    const Tensor6_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n> &a )
{
  return a * b;
}


  }
