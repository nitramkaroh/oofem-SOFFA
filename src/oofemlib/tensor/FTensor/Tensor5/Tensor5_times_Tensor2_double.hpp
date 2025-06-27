#pragma once

namespace FTensor
{
  // A(i,j,k,l) * B(m,n) double contraction
  template <class A, class B, class T, class U, int Dim0, int Dim1, int Dim2,
            int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l,
            char m, char n, char o, int DimA, int DimB, int DimC, int DimX, int DimY, char a,
            char b, char c, char x, char y>
  class Tensor5_times_Tensor2_double
  {
    Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> iterA;
    Tensor2_Expr<B, U, Dim5, Dim6, n, o> iterB;

  public:
    typename promote<T, U>::V operator()(const int N1, const int N2, const int N3) const
    {
      typename promote<T, U>::V result(0);
      for(int xx = 0; xx < DimX; ++xx)
        for(int yy = 0; yy < DimY; ++yy)
          {
            // Permutation is where the indices get checked.
            result += Permutation5<DimA, DimB, DimC, DimX, DimY, a, b, c, x, y>().eval(
                        iterA, N1, N2, N3, xx, yy)
                      * Permutation2<DimX, DimY, x, y>().eval(iterB, xx, yy);
          }
      return result;
    }

    Tensor5_times_Tensor2_double(
      const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &iter_a,
      const Tensor2_Expr<B, U, Dim5, Dim6, n, o> &iter_b)
        : iterA(iter_a), iterB(iter_b)
    {}
  };

  // A(i, j, k, l, m, )*B(j, i, )
  template <class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, i, j, k, l, m>

  auto operator*( const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
      const Tensor2_Expr<B, U, Dim1, Dim0, j, i> &b )
  {
    using TensorExpr
        = Tensor5_times_Tensor2_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim1, Dim0, i, j, k, l, m, j, i, Dim2, Dim3, Dim4, Dim1, Dim0, k, l, m, j, i>;
    return Tensor3_Expr<TensorExpr, typename promote<T, U>::V, Dim2, Dim3, Dim4, k, l, m>( TensorExpr( a, b ) );
  }

  // B(j, i, )*A(i, j, k, l, m, )
  template <class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, i, j, k, l, m>

  auto operator*( const Tensor2_Expr<B, U, Dim1, Dim0, j, i> &b,
      const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a )
  {
    return a * b;
  }

  // A(i, j, k, l, m, )*B(i, j, )
  template <class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, i, j, k, l, m>

  auto operator*( const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
      const Tensor2_Expr<B, U, Dim0, Dim1, i, j> &b )
  {
    using TensorExpr
        = Tensor5_times_Tensor2_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim0, Dim1, i, j, k, l, m, i, j, Dim2, Dim3, Dim4, Dim0, Dim1, k, l, m, i, j>;
    return Tensor3_Expr<TensorExpr, typename promote<T, U>::V, Dim2, Dim3, Dim4, k, l, m>( TensorExpr( a, b ) );
  }

  // B(i, j, )*A(i, j, k, l, m, )
  template <class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, i, j, k, l, m>

  auto operator*( const Tensor2_Expr<B, U, Dim0, Dim1, i, j> &b,
      const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a )
  {
    return a * b;
  }

  // A(i, j, k, l, m, )*B(k, i, )
  template <class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, i, j, k, l, m>

  auto operator*( const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
      const Tensor2_Expr<B, U, Dim2, Dim0, k, i> &b )
  {
    using TensorExpr
        = Tensor5_times_Tensor2_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim2, Dim0, i, j, k, l, m, k, i, Dim1, Dim3, Dim4, Dim2, Dim0, j, l, m, k, i>;
    return Tensor3_Expr<TensorExpr, typename promote<T, U>::V, Dim1, Dim3, Dim4, j, l, m>( TensorExpr( a, b ) );
  }

  // B(k, i, )*A(i, j, k, l, m, )
  template <class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, i, j, k, l, m>

  auto operator*( const Tensor2_Expr<B, U, Dim2, Dim0, k, i> &b,
      const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a )
  {
    return a * b;
  }

  // A(i, j, k, l, m, )*B(i, k, )
  template <class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, i, j, k, l, m>

  auto operator*( const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
      const Tensor2_Expr<B, U, Dim0, Dim2, i, k> &b )
  {
    using TensorExpr
        = Tensor5_times_Tensor2_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim0, Dim2, i, j, k, l, m, i, k, Dim1, Dim3, Dim4, Dim0, Dim2, j, l, m, i, k>;
    return Tensor3_Expr<TensorExpr, typename promote<T, U>::V, Dim1, Dim3, Dim4, j, l, m>( TensorExpr( a, b ) );
  }

  // B(i, k, )*A(i, j, k, l, m, )
  template <class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, i, j, k, l, m>

  auto operator*( const Tensor2_Expr<B, U, Dim0, Dim2, i, k> &b,
      const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a )
  {
    return a * b;
  }

  // A(i, j, k, l, m, )*B(l, i, )
  template <class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, i, j, k, l, m>

  auto operator*( const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
      const Tensor2_Expr<B, U, Dim3, Dim0, l, i> &b )
  {
    using TensorExpr
        = Tensor5_times_Tensor2_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim3, Dim0, i, j, k, l, m, l, i, Dim1, Dim2, Dim4, Dim3, Dim0, j, k, m, l, i>;
    return Tensor3_Expr<TensorExpr, typename promote<T, U>::V, Dim1, Dim2, Dim4, j, k, m>( TensorExpr( a, b ) );
  }

  // B(l, i, )*A(i, j, k, l, m, )
  template <class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, i, j, k, l, m>

  auto operator*( const Tensor2_Expr<B, U, Dim3, Dim0, l, i> &b,
      const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a )
  {
    return a * b;
  }

  // A(i, j, k, l, m, )*B(i, l, )
  template <class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, i, j, k, l, m>

  auto operator*( const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
      const Tensor2_Expr<B, U, Dim0, Dim3, i, l> &b )
  {
    using TensorExpr
        = Tensor5_times_Tensor2_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim0, Dim3, i, j, k, l, m, i, l, Dim1, Dim2, Dim4, Dim0, Dim3, j, k, m, i, l>;
    return Tensor3_Expr<TensorExpr, typename promote<T, U>::V, Dim1, Dim2, Dim4, j, k, m>( TensorExpr( a, b ) );
  }

  // B(i, l, )*A(i, j, k, l, m, )
  template <class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, i, j, k, l, m>

  auto operator*( const Tensor2_Expr<B, U, Dim0, Dim3, i, l> &b,
      const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a )
  {
    return a * b;
  }

  // A(i, j, k, l, m, )*B(m, i, )
  template <class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, i, j, k, l, m>

  auto operator*( const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
      const Tensor2_Expr<B, U, Dim4, Dim0, m, i> &b )
  {
    using TensorExpr
        = Tensor5_times_Tensor2_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim4, Dim0, i, j, k, l, m, m, i, Dim1, Dim2, Dim3, Dim4, Dim0, j, k, l, m, i>;
    return Tensor3_Expr<TensorExpr, typename promote<T, U>::V, Dim1, Dim2, Dim3, j, k, l>( TensorExpr( a, b ) );
  }

  // B(m, i, )*A(i, j, k, l, m, )
  template <class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, i, j, k, l, m>

  auto operator*( const Tensor2_Expr<B, U, Dim4, Dim0, m, i> &b,
      const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a )
  {
    return a * b;
  }

  // A(i, j, k, l, m, )*B(i, m, )
  template <class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, i, j, k, l, m>

  auto operator*( const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
      const Tensor2_Expr<B, U, Dim0, Dim4, i, m> &b )
  {
    using TensorExpr
        = Tensor5_times_Tensor2_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim0, Dim4, i, j, k, l, m, i, m, Dim1, Dim2, Dim3, Dim0, Dim4, j, k, l, i, m>;
    return Tensor3_Expr<TensorExpr, typename promote<T, U>::V, Dim1, Dim2, Dim3, j, k, l>( TensorExpr( a, b ) );
  }

  // B(i, m, )*A(i, j, k, l, m, )
  template <class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, i, j, k, l, m>

  auto operator*( const Tensor2_Expr<B, U, Dim0, Dim4, i, m> &b,
      const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a )
  {
    return a * b;
  }

  // A(i, j, k, l, m, )*B(k, j, )
  template <class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, i, j, k, l, m>

  auto operator*( const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
      const Tensor2_Expr<B, U, Dim2, Dim1, k, j> &b )
  {
    using TensorExpr
        = Tensor5_times_Tensor2_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim2, Dim1, i, j, k, l, m, k, j, Dim0, Dim3, Dim4, Dim2, Dim1, i, l, m, k, j>;
    return Tensor3_Expr<TensorExpr, typename promote<T, U>::V, Dim0, Dim3, Dim4, i, l, m>( TensorExpr( a, b ) );
  }

  // B(k, j, )*A(i, j, k, l, m, )
  template <class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, i, j, k, l, m>

  auto operator*( const Tensor2_Expr<B, U, Dim2, Dim1, k, j> &b,
      const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a )
  {
    return a * b;
  }

  // A(i, j, k, l, m, )*B(j, k, )
  template <class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, i, j, k, l, m>

  auto operator*( const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
      const Tensor2_Expr<B, U, Dim1, Dim2, j, k> &b )
  {
    using TensorExpr
        = Tensor5_times_Tensor2_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim1, Dim2, i, j, k, l, m, j, k, Dim0, Dim3, Dim4, Dim1, Dim2, i, l, m, j, k>;
    return Tensor3_Expr<TensorExpr, typename promote<T, U>::V, Dim0, Dim3, Dim4, i, l, m>( TensorExpr( a, b ) );
  }

  // B(j, k, )*A(i, j, k, l, m, )
  template <class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, i, j, k, l, m>

  auto operator*( const Tensor2_Expr<B, U, Dim1, Dim2, j, k> &b,
      const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a )
  {
    return a * b;
  }

  // A(i, j, k, l, m, )*B(l, j, )
  template <class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, i, j, k, l, m>

  auto operator*( const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
      const Tensor2_Expr<B, U, Dim3, Dim1, l, j> &b )
  {
    using TensorExpr
        = Tensor5_times_Tensor2_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim3, Dim1, i, j, k, l, m, l, j, Dim0, Dim2, Dim4, Dim3, Dim1, i, k, m, l, j>;
    return Tensor3_Expr<TensorExpr, typename promote<T, U>::V, Dim0, Dim2, Dim4, i, k, m>( TensorExpr( a, b ) );
  }

  // B(l, j, )*A(i, j, k, l, m, )
  template <class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, i, j, k, l, m>

  auto operator*( const Tensor2_Expr<B, U, Dim3, Dim1, l, j> &b,
      const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a )
  {
    return a * b;
  }

  // A(i, j, k, l, m, )*B(j, l, )
  template <class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, i, j, k, l, m>

  auto operator*( const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
      const Tensor2_Expr<B, U, Dim1, Dim3, j, l> &b )
  {
    using TensorExpr
        = Tensor5_times_Tensor2_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim1, Dim3, i, j, k, l, m, j, l, Dim0, Dim2, Dim4, Dim1, Dim3, i, k, m, j, l>;
    return Tensor3_Expr<TensorExpr, typename promote<T, U>::V, Dim0, Dim2, Dim4, i, k, m>( TensorExpr( a, b ) );
  }

  // B(j, l, )*A(i, j, k, l, m, )
  template <class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, i, j, k, l, m>

  auto operator*( const Tensor2_Expr<B, U, Dim1, Dim3, j, l> &b,
      const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a )
  {
    return a * b;
  }

  // A(i, j, k, l, m, )*B(m, j, )
  template <class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, i, j, k, l, m>

  auto operator*( const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
      const Tensor2_Expr<B, U, Dim4, Dim1, m, j> &b )
  {
    using TensorExpr
        = Tensor5_times_Tensor2_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim4, Dim1, i, j, k, l, m, m, j, Dim0, Dim2, Dim3, Dim4, Dim1, i, k, l, m, j>;
    return Tensor3_Expr<TensorExpr, typename promote<T, U>::V, Dim0, Dim2, Dim3, i, k, l>( TensorExpr( a, b ) );
  }

  // B(m, j, )*A(i, j, k, l, m, )
  template <class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, i, j, k, l, m>

  auto operator*( const Tensor2_Expr<B, U, Dim4, Dim1, m, j> &b,
      const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a )
  {
    return a * b;
  }

  // A(i, j, k, l, m, )*B(j, m, )
  template <class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, i, j, k, l, m>

  auto operator*( const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
      const Tensor2_Expr<B, U, Dim1, Dim4, j, m> &b )
  {
    using TensorExpr
        = Tensor5_times_Tensor2_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim1, Dim4, i, j, k, l, m, j, m, Dim0, Dim2, Dim3, Dim1, Dim4, i, k, l, j, m>;
    return Tensor3_Expr<TensorExpr, typename promote<T, U>::V, Dim0, Dim2, Dim3, i, k, l>( TensorExpr( a, b ) );
  }

  // B(j, m, )*A(i, j, k, l, m, )
  template <class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, i, j, k, l, m>

  auto operator*( const Tensor2_Expr<B, U, Dim1, Dim4, j, m> &b,
      const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a )
  {
    return a * b;
  }

  // A(i, j, k, l, m, )*B(l, k, )
  template <class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, i, j, k, l, m>

  auto operator*( const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
      const Tensor2_Expr<B, U, Dim3, Dim2, l, k> &b )
  {
    using TensorExpr
        = Tensor5_times_Tensor2_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim3, Dim2, i, j, k, l, m, l, k, Dim0, Dim1, Dim4, Dim3, Dim2, i, j, m, l, k>;
    return Tensor3_Expr<TensorExpr, typename promote<T, U>::V, Dim0, Dim1, Dim4, i, j, m>( TensorExpr( a, b ) );
  }

  // B(l, k, )*A(i, j, k, l, m, )
  template <class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, i, j, k, l, m>

  auto operator*( const Tensor2_Expr<B, U, Dim3, Dim2, l, k> &b,
      const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a )
  {
    return a * b;
  }

  // A(i, j, k, l, m, )*B(k, l, )
  template <class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, i, j, k, l, m>

  auto operator*( const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
      const Tensor2_Expr<B, U, Dim2, Dim3, k, l> &b )
  {
    using TensorExpr
        = Tensor5_times_Tensor2_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim2, Dim3, i, j, k, l, m, k, l, Dim0, Dim1, Dim4, Dim2, Dim3, i, j, m, k, l>;
    return Tensor3_Expr<TensorExpr, typename promote<T, U>::V, Dim0, Dim1, Dim4, i, j, m>( TensorExpr( a, b ) );
  }

  // B(k, l, )*A(i, j, k, l, m, )
  template <class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, i, j, k, l, m>

  auto operator*( const Tensor2_Expr<B, U, Dim2, Dim3, k, l> &b,
      const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a )
  {
    return a * b;
  }

  // A(i, j, k, l, m, )*B(m, k, )
  template <class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, i, j, k, l, m>

  auto operator*( const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
      const Tensor2_Expr<B, U, Dim4, Dim2, m, k> &b )
  {
    using TensorExpr
        = Tensor5_times_Tensor2_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim4, Dim2, i, j, k, l, m, m, k, Dim0, Dim1, Dim3, Dim4, Dim2, i, j, l, m, k>;
    return Tensor3_Expr<TensorExpr, typename promote<T, U>::V, Dim0, Dim1, Dim3, i, j, l>( TensorExpr( a, b ) );
  }

  // B(m, k, )*A(i, j, k, l, m, )
  template <class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, i, j, k, l, m>

  auto operator*( const Tensor2_Expr<B, U, Dim4, Dim2, m, k> &b,
      const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a )
  {
    return a * b;
  }

  // A(i, j, k, l, m, )*B(k, m, )
  template <class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, i, j, k, l, m>

  auto operator*( const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
      const Tensor2_Expr<B, U, Dim2, Dim4, k, m> &b )
  {
    using TensorExpr
        = Tensor5_times_Tensor2_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim2, Dim4, i, j, k, l, m, k, m, Dim0, Dim1, Dim3, Dim2, Dim4, i, j, l, k, m>;
    return Tensor3_Expr<TensorExpr, typename promote<T, U>::V, Dim0, Dim1, Dim3, i, j, l>( TensorExpr( a, b ) );
  }

  // B(k, m, )*A(i, j, k, l, m, )
  template <class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, i, j, k, l, m>

  auto operator*( const Tensor2_Expr<B, U, Dim2, Dim4, k, m> &b,
      const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a )
  {
    return a * b;
  }

  // A(i, j, k, l, m, )*B(m, l, )
  template <class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, i, j, k, l, m>

  auto operator*( const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
      const Tensor2_Expr<B, U, Dim4, Dim3, m, l> &b )
  {
    using TensorExpr
        = Tensor5_times_Tensor2_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim4, Dim3, i, j, k, l, m, m, l, Dim0, Dim1, Dim2, Dim4, Dim3, i, j, k, m, l>;
    return Tensor3_Expr<TensorExpr, typename promote<T, U>::V, Dim0, Dim1, Dim2, i, j, k>( TensorExpr( a, b ) );
  }

  // B(m, l, )*A(i, j, k, l, m, )
  template <class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, i, j, k, l, m>

  auto operator*( const Tensor2_Expr<B, U, Dim4, Dim3, m, l> &b,
      const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a )
  {
    return a * b;
  }

  // A(i, j, k, l, m, )*B(l, m, )
  template <class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, i, j, k, l, m>

  auto operator*( const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
      const Tensor2_Expr<B, U, Dim3, Dim4, l, m> &b )
  {
    using TensorExpr
        = Tensor5_times_Tensor2_double<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim3, Dim4, i, j, k, l, m, l, m, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m>;
    return Tensor3_Expr<TensorExpr, typename promote<T, U>::V, Dim0, Dim1, Dim2, i, j, k>( TensorExpr( a, b ) );
  }

  // B(l, m, )*A(i, j, k, l, m, )
  template <class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, i, j, k, l, m>

  auto operator*( const Tensor2_Expr<B, U, Dim3, Dim4, l, m> &b,
      const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a )
  {
    return a * b;
  }
  }
