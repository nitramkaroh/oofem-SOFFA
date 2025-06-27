#pragma once

namespace FTensor
{
  // A(i,j,k,l,m) * B(n,o,p) triple contraction
  template <class A, class B, class T, class U, int Dim0, int Dim1, int Dim2,
            int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k,
            char l, char m, char n, char o, char p, int DimA, int DimB, int DimX, int DimY,
            int DimZ, char a, char b, char x, char y, char z>
  class Tensor5_times_Tensor3_triple
  {
    Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> iterA;
    Tensor3_Expr<B, U, Dim5, Dim6, Dim7, n, o, p> iterB;

  public:
    typename promote<T, U>::V operator()(const int N1, const int N2) const
    {
      typename promote<T, U>::V result(0);
      for(int xx = 0; xx < DimX; ++xx)
        for(int yy = 0; yy < DimY; ++yy)
          for(int zz = 0; zz < DimZ; ++zz)
            {
              // Permutation is where the indices get checked.
              result
                += Permutation5<DimA, DimB, DimX, DimY, DimZ, a, b, x, y, z>().eval(
                     iterA, N1, N2, xx, yy, zz)
                   * Permutation3<DimX, DimY, DimZ, x, y, z>().eval(iterB, xx,
                                                                    yy, zz);
            }
      return result;
    }

    Tensor5_times_Tensor3_triple(
      const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &iter_a,
      const Tensor3_Expr<B, U, Dim5, Dim6, Dim7, n, o, p> &iter_b)
        : iterA(iter_a), iterB(iter_b)
    {}
  };

  // A(i, j, k, l, m, )*B(k, j, i, )
  template <class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, i, j, k, l, m>

  auto operator*( const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
      const Tensor3_Expr<B, U, Dim2, Dim1, Dim0, k, j, i> &b )
  {
    using TensorExpr
        = Tensor5_times_Tensor3_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim2, Dim1, Dim0, i, j, k, l, m, k, j, i, Dim3, Dim4, Dim2, Dim1, Dim0, l, m, k, j, i>;
    return Tensor2_Expr<TensorExpr, typename promote<T, U>::V, Dim3, Dim4, l, m>( TensorExpr( a, b ) );
  }

  // B(k, j, i, )*A(i, j, k, l, m, )
  template <class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, i, j, k, l, m>

  auto operator*( const Tensor3_Expr<B, U, Dim2, Dim1, Dim0, k, j, i> &b,
      const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a )
  {
    return a * b;
  }

  // A(i, j, k, l, m, )*B(j, k, i, )
  template <class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, i, j, k, l, m>

  auto operator*( const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
      const Tensor3_Expr<B, U, Dim1, Dim2, Dim0, j, k, i> &b )
  {
    using TensorExpr
        = Tensor5_times_Tensor3_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim1, Dim2, Dim0, i, j, k, l, m, j, k, i, Dim3, Dim4, Dim1, Dim2, Dim0, l, m, j, k, i>;
    return Tensor2_Expr<TensorExpr, typename promote<T, U>::V, Dim3, Dim4, l, m>( TensorExpr( a, b ) );
  }

  // B(j, k, i, )*A(i, j, k, l, m, )
  template <class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, i, j, k, l, m>

  auto operator*( const Tensor3_Expr<B, U, Dim1, Dim2, Dim0, j, k, i> &b,
      const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a )
  {
    return a * b;
  }

  // A(i, j, k, l, m, )*B(k, i, j, )
  template <class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, i, j, k, l, m>

  auto operator*( const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
      const Tensor3_Expr<B, U, Dim2, Dim0, Dim1, k, i, j> &b )
  {
    using TensorExpr
        = Tensor5_times_Tensor3_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim2, Dim0, Dim1, i, j, k, l, m, k, i, j, Dim3, Dim4, Dim2, Dim0, Dim1, l, m, k, i, j>;
    return Tensor2_Expr<TensorExpr, typename promote<T, U>::V, Dim3, Dim4, l, m>( TensorExpr( a, b ) );
  }

  // B(k, i, j, )*A(i, j, k, l, m, )
  template <class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, i, j, k, l, m>

  auto operator*( const Tensor3_Expr<B, U, Dim2, Dim0, Dim1, k, i, j> &b,
      const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a )
  {
    return a * b;
  }

  // A(i, j, k, l, m, )*B(j, i, k, )
  template <class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, i, j, k, l, m>

  auto operator*( const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
      const Tensor3_Expr<B, U, Dim1, Dim0, Dim2, j, i, k> &b )
  {
    using TensorExpr
        = Tensor5_times_Tensor3_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim1, Dim0, Dim2, i, j, k, l, m, j, i, k, Dim3, Dim4, Dim1, Dim0, Dim2, l, m, j, i, k>;
    return Tensor2_Expr<TensorExpr, typename promote<T, U>::V, Dim3, Dim4, l, m>( TensorExpr( a, b ) );
  }

  // B(j, i, k, )*A(i, j, k, l, m, )
  template <class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, i, j, k, l, m>

  auto operator*( const Tensor3_Expr<B, U, Dim1, Dim0, Dim2, j, i, k> &b,
      const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a )
  {
    return a * b;
  }

  // A(i, j, k, l, m, )*B(i, k, j, )
  template <class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, i, j, k, l, m>

  auto operator*( const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
      const Tensor3_Expr<B, U, Dim0, Dim2, Dim1, i, k, j> &b )
  {
    using TensorExpr
        = Tensor5_times_Tensor3_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim0, Dim2, Dim1, i, j, k, l, m, i, k, j, Dim3, Dim4, Dim0, Dim2, Dim1, l, m, i, k, j>;
    return Tensor2_Expr<TensorExpr, typename promote<T, U>::V, Dim3, Dim4, l, m>( TensorExpr( a, b ) );
  }

  // B(i, k, j, )*A(i, j, k, l, m, )
  template <class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, i, j, k, l, m>

  auto operator*( const Tensor3_Expr<B, U, Dim0, Dim2, Dim1, i, k, j> &b,
      const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a )
  {
    return a * b;
  }

  // A(i, j, k, l, m, )*B(i, j, k, )
  template <class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, i, j, k, l, m>

  auto operator*( const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
      const Tensor3_Expr<B, U, Dim0, Dim1, Dim2, i, j, k> &b )
  {
    using TensorExpr
        = Tensor5_times_Tensor3_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim0, Dim1, Dim2, i, j, k, l, m, i, j, k, Dim3, Dim4, Dim0, Dim1, Dim2, l, m, i, j, k>;
    return Tensor2_Expr<TensorExpr, typename promote<T, U>::V, Dim3, Dim4, l, m>( TensorExpr( a, b ) );
  }

  // B(i, j, k, )*A(i, j, k, l, m, )
  template <class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, i, j, k, l, m>

  auto operator*( const Tensor3_Expr<B, U, Dim0, Dim1, Dim2, i, j, k> &b,
      const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a )
  {
    return a * b;
  }

  // A(i, j, k, l, m, )*B(l, j, i, )
  template <class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, i, j, k, l, m>

  auto operator*( const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
      const Tensor3_Expr<B, U, Dim3, Dim1, Dim0, l, j, i> &b )
  {
    using TensorExpr
        = Tensor5_times_Tensor3_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim3, Dim1, Dim0, i, j, k, l, m, l, j, i, Dim2, Dim4, Dim3, Dim1, Dim0, k, m, l, j, i>;
    return Tensor2_Expr<TensorExpr, typename promote<T, U>::V, Dim2, Dim4, k, m>( TensorExpr( a, b ) );
  }

  // B(l, j, i, )*A(i, j, k, l, m, )
  template <class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, i, j, k, l, m>

  auto operator*( const Tensor3_Expr<B, U, Dim3, Dim1, Dim0, l, j, i> &b,
      const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a )
  {
    return a * b;
  }

  // A(i, j, k, l, m, )*B(j, l, i, )
  template <class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, i, j, k, l, m>

  auto operator*( const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
      const Tensor3_Expr<B, U, Dim1, Dim3, Dim0, j, l, i> &b )
  {
    using TensorExpr
        = Tensor5_times_Tensor3_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim1, Dim3, Dim0, i, j, k, l, m, j, l, i, Dim2, Dim4, Dim1, Dim3, Dim0, k, m, j, l, i>;
    return Tensor2_Expr<TensorExpr, typename promote<T, U>::V, Dim2, Dim4, k, m>( TensorExpr( a, b ) );
  }

  // B(j, l, i, )*A(i, j, k, l, m, )
  template <class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, i, j, k, l, m>

  auto operator*( const Tensor3_Expr<B, U, Dim1, Dim3, Dim0, j, l, i> &b,
      const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a )
  {
    return a * b;
  }

  // A(i, j, k, l, m, )*B(l, i, j, )
  template <class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, i, j, k, l, m>

  auto operator*( const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
      const Tensor3_Expr<B, U, Dim3, Dim0, Dim1, l, i, j> &b )
  {
    using TensorExpr
        = Tensor5_times_Tensor3_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim3, Dim0, Dim1, i, j, k, l, m, l, i, j, Dim2, Dim4, Dim3, Dim0, Dim1, k, m, l, i, j>;
    return Tensor2_Expr<TensorExpr, typename promote<T, U>::V, Dim2, Dim4, k, m>( TensorExpr( a, b ) );
  }

  // B(l, i, j, )*A(i, j, k, l, m, )
  template <class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, i, j, k, l, m>

  auto operator*( const Tensor3_Expr<B, U, Dim3, Dim0, Dim1, l, i, j> &b,
      const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a )
  {
    return a * b;
  }

  // A(i, j, k, l, m, )*B(j, i, l, )
  template <class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, i, j, k, l, m>

  auto operator*( const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
      const Tensor3_Expr<B, U, Dim1, Dim0, Dim3, j, i, l> &b )
  {
    using TensorExpr
        = Tensor5_times_Tensor3_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim1, Dim0, Dim3, i, j, k, l, m, j, i, l, Dim2, Dim4, Dim1, Dim0, Dim3, k, m, j, i, l>;
    return Tensor2_Expr<TensorExpr, typename promote<T, U>::V, Dim2, Dim4, k, m>( TensorExpr( a, b ) );
  }

  // B(j, i, l, )*A(i, j, k, l, m, )
  template <class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, i, j, k, l, m>

  auto operator*( const Tensor3_Expr<B, U, Dim1, Dim0, Dim3, j, i, l> &b,
      const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a )
  {
    return a * b;
  }

  // A(i, j, k, l, m, )*B(i, l, j, )
  template <class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, i, j, k, l, m>

  auto operator*( const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
      const Tensor3_Expr<B, U, Dim0, Dim3, Dim1, i, l, j> &b )
  {
    using TensorExpr
        = Tensor5_times_Tensor3_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim0, Dim3, Dim1, i, j, k, l, m, i, l, j, Dim2, Dim4, Dim0, Dim3, Dim1, k, m, i, l, j>;
    return Tensor2_Expr<TensorExpr, typename promote<T, U>::V, Dim2, Dim4, k, m>( TensorExpr( a, b ) );
  }

  // B(i, l, j, )*A(i, j, k, l, m, )
  template <class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, i, j, k, l, m>

  auto operator*( const Tensor3_Expr<B, U, Dim0, Dim3, Dim1, i, l, j> &b,
      const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a )
  {
    return a * b;
  }

  // A(i, j, k, l, m, )*B(i, j, l, )
  template <class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, i, j, k, l, m>

  auto operator*( const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
      const Tensor3_Expr<B, U, Dim0, Dim1, Dim3, i, j, l> &b )
  {
    using TensorExpr
        = Tensor5_times_Tensor3_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim0, Dim1, Dim3, i, j, k, l, m, i, j, l, Dim2, Dim4, Dim0, Dim1, Dim3, k, m, i, j, l>;
    return Tensor2_Expr<TensorExpr, typename promote<T, U>::V, Dim2, Dim4, k, m>( TensorExpr( a, b ) );
  }

  // B(i, j, l, )*A(i, j, k, l, m, )
  template <class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, i, j, k, l, m>

  auto operator*( const Tensor3_Expr<B, U, Dim0, Dim1, Dim3, i, j, l> &b,
      const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a )
  {
    return a * b;
  }

  // A(i, j, k, l, m, )*B(m, j, i, )
  template <class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, i, j, k, l, m>

  auto operator*( const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
      const Tensor3_Expr<B, U, Dim4, Dim1, Dim0, m, j, i> &b )
  {
    using TensorExpr
        = Tensor5_times_Tensor3_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim4, Dim1, Dim0, i, j, k, l, m, m, j, i, Dim2, Dim3, Dim4, Dim1, Dim0, k, l, m, j, i>;
    return Tensor2_Expr<TensorExpr, typename promote<T, U>::V, Dim2, Dim3, k, l>( TensorExpr( a, b ) );
  }

  // B(m, j, i, )*A(i, j, k, l, m, )
  template <class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, i, j, k, l, m>

  auto operator*( const Tensor3_Expr<B, U, Dim4, Dim1, Dim0, m, j, i> &b,
      const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a )
  {
    return a * b;
  }

  // A(i, j, k, l, m, )*B(j, m, i, )
  template <class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, i, j, k, l, m>

  auto operator*( const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
      const Tensor3_Expr<B, U, Dim1, Dim4, Dim0, j, m, i> &b )
  {
    using TensorExpr
        = Tensor5_times_Tensor3_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim1, Dim4, Dim0, i, j, k, l, m, j, m, i, Dim2, Dim3, Dim1, Dim4, Dim0, k, l, j, m, i>;
    return Tensor2_Expr<TensorExpr, typename promote<T, U>::V, Dim2, Dim3, k, l>( TensorExpr( a, b ) );
  }

  // B(j, m, i, )*A(i, j, k, l, m, )
  template <class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, i, j, k, l, m>

  auto operator*( const Tensor3_Expr<B, U, Dim1, Dim4, Dim0, j, m, i> &b,
      const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a )
  {
    return a * b;
  }

  // A(i, j, k, l, m, )*B(m, i, j, )
  template <class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, i, j, k, l, m>

  auto operator*( const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
      const Tensor3_Expr<B, U, Dim4, Dim0, Dim1, m, i, j> &b )
  {
    using TensorExpr
        = Tensor5_times_Tensor3_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim4, Dim0, Dim1, i, j, k, l, m, m, i, j, Dim2, Dim3, Dim4, Dim0, Dim1, k, l, m, i, j>;
    return Tensor2_Expr<TensorExpr, typename promote<T, U>::V, Dim2, Dim3, k, l>( TensorExpr( a, b ) );
  }

  // B(m, i, j, )*A(i, j, k, l, m, )
  template <class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, i, j, k, l, m>

  auto operator*( const Tensor3_Expr<B, U, Dim4, Dim0, Dim1, m, i, j> &b,
      const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a )
  {
    return a * b;
  }

  // A(i, j, k, l, m, )*B(j, i, m, )
  template <class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, i, j, k, l, m>

  auto operator*( const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
      const Tensor3_Expr<B, U, Dim1, Dim0, Dim4, j, i, m> &b )
  {
    using TensorExpr
        = Tensor5_times_Tensor3_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim1, Dim0, Dim4, i, j, k, l, m, j, i, m, Dim2, Dim3, Dim1, Dim0, Dim4, k, l, j, i, m>;
    return Tensor2_Expr<TensorExpr, typename promote<T, U>::V, Dim2, Dim3, k, l>( TensorExpr( a, b ) );
  }

  // B(j, i, m, )*A(i, j, k, l, m, )
  template <class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, i, j, k, l, m>

  auto operator*( const Tensor3_Expr<B, U, Dim1, Dim0, Dim4, j, i, m> &b,
      const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a )
  {
    return a * b;
  }

  // A(i, j, k, l, m, )*B(i, m, j, )
  template <class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, i, j, k, l, m>

  auto operator*( const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
      const Tensor3_Expr<B, U, Dim0, Dim4, Dim1, i, m, j> &b )
  {
    using TensorExpr
        = Tensor5_times_Tensor3_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim0, Dim4, Dim1, i, j, k, l, m, i, m, j, Dim2, Dim3, Dim0, Dim4, Dim1, k, l, i, m, j>;
    return Tensor2_Expr<TensorExpr, typename promote<T, U>::V, Dim2, Dim3, k, l>( TensorExpr( a, b ) );
  }

  // B(i, m, j, )*A(i, j, k, l, m, )
  template <class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, i, j, k, l, m>

  auto operator*( const Tensor3_Expr<B, U, Dim0, Dim4, Dim1, i, m, j> &b,
      const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a )
  {
    return a * b;
  }

  // A(i, j, k, l, m, )*B(i, j, m, )
  template <class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, i, j, k, l, m>

  auto operator*( const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
      const Tensor3_Expr<B, U, Dim0, Dim1, Dim4, i, j, m> &b )
  {
    using TensorExpr
        = Tensor5_times_Tensor3_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim0, Dim1, Dim4, i, j, k, l, m, i, j, m, Dim2, Dim3, Dim0, Dim1, Dim4, k, l, i, j, m>;
    return Tensor2_Expr<TensorExpr, typename promote<T, U>::V, Dim2, Dim3, k, l>( TensorExpr( a, b ) );
  }

  // B(i, j, m, )*A(i, j, k, l, m, )
  template <class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, i, j, k, l, m>

  auto operator*( const Tensor3_Expr<B, U, Dim0, Dim1, Dim4, i, j, m> &b,
      const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a )
  {
    return a * b;
  }

  // A(i, j, k, l, m, )*B(l, k, i, )
  template <class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, i, j, k, l, m>

  auto operator*( const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
      const Tensor3_Expr<B, U, Dim3, Dim2, Dim0, l, k, i> &b )
  {
    using TensorExpr
        = Tensor5_times_Tensor3_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim3, Dim2, Dim0, i, j, k, l, m, l, k, i, Dim1, Dim4, Dim3, Dim2, Dim0, j, m, l, k, i>;
    return Tensor2_Expr<TensorExpr, typename promote<T, U>::V, Dim1, Dim4, j, m>( TensorExpr( a, b ) );
  }

  // B(l, k, i, )*A(i, j, k, l, m, )
  template <class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, i, j, k, l, m>

  auto operator*( const Tensor3_Expr<B, U, Dim3, Dim2, Dim0, l, k, i> &b,
      const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a )
  {
    return a * b;
  }

  // A(i, j, k, l, m, )*B(k, l, i, )
  template <class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, i, j, k, l, m>

  auto operator*( const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
      const Tensor3_Expr<B, U, Dim2, Dim3, Dim0, k, l, i> &b )
  {
    using TensorExpr
        = Tensor5_times_Tensor3_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim2, Dim3, Dim0, i, j, k, l, m, k, l, i, Dim1, Dim4, Dim2, Dim3, Dim0, j, m, k, l, i>;
    return Tensor2_Expr<TensorExpr, typename promote<T, U>::V, Dim1, Dim4, j, m>( TensorExpr( a, b ) );
  }

  // B(k, l, i, )*A(i, j, k, l, m, )
  template <class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, i, j, k, l, m>

  auto operator*( const Tensor3_Expr<B, U, Dim2, Dim3, Dim0, k, l, i> &b,
      const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a )
  {
    return a * b;
  }

  // A(i, j, k, l, m, )*B(l, i, k, )
  template <class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, i, j, k, l, m>

  auto operator*( const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
      const Tensor3_Expr<B, U, Dim3, Dim0, Dim2, l, i, k> &b )
  {
    using TensorExpr
        = Tensor5_times_Tensor3_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim3, Dim0, Dim2, i, j, k, l, m, l, i, k, Dim1, Dim4, Dim3, Dim0, Dim2, j, m, l, i, k>;
    return Tensor2_Expr<TensorExpr, typename promote<T, U>::V, Dim1, Dim4, j, m>( TensorExpr( a, b ) );
  }

  // B(l, i, k, )*A(i, j, k, l, m, )
  template <class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, i, j, k, l, m>

  auto operator*( const Tensor3_Expr<B, U, Dim3, Dim0, Dim2, l, i, k> &b,
      const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a )
  {
    return a * b;
  }

  // A(i, j, k, l, m, )*B(k, i, l, )
  template <class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, i, j, k, l, m>

  auto operator*( const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
      const Tensor3_Expr<B, U, Dim2, Dim0, Dim3, k, i, l> &b )
  {
    using TensorExpr
        = Tensor5_times_Tensor3_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim2, Dim0, Dim3, i, j, k, l, m, k, i, l, Dim1, Dim4, Dim2, Dim0, Dim3, j, m, k, i, l>;
    return Tensor2_Expr<TensorExpr, typename promote<T, U>::V, Dim1, Dim4, j, m>( TensorExpr( a, b ) );
  }

  // B(k, i, l, )*A(i, j, k, l, m, )
  template <class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, i, j, k, l, m>

  auto operator*( const Tensor3_Expr<B, U, Dim2, Dim0, Dim3, k, i, l> &b,
      const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a )
  {
    return a * b;
  }

  // A(i, j, k, l, m, )*B(i, l, k, )
  template <class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, i, j, k, l, m>

  auto operator*( const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
      const Tensor3_Expr<B, U, Dim0, Dim3, Dim2, i, l, k> &b )
  {
    using TensorExpr
        = Tensor5_times_Tensor3_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim0, Dim3, Dim2, i, j, k, l, m, i, l, k, Dim1, Dim4, Dim0, Dim3, Dim2, j, m, i, l, k>;
    return Tensor2_Expr<TensorExpr, typename promote<T, U>::V, Dim1, Dim4, j, m>( TensorExpr( a, b ) );
  }

  // B(i, l, k, )*A(i, j, k, l, m, )
  template <class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, i, j, k, l, m>

  auto operator*( const Tensor3_Expr<B, U, Dim0, Dim3, Dim2, i, l, k> &b,
      const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a )
  {
    return a * b;
  }

  // A(i, j, k, l, m, )*B(i, k, l, )
  template <class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, i, j, k, l, m>

  auto operator*( const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
      const Tensor3_Expr<B, U, Dim0, Dim2, Dim3, i, k, l> &b )
  {
    using TensorExpr
        = Tensor5_times_Tensor3_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim0, Dim2, Dim3, i, j, k, l, m, i, k, l, Dim1, Dim4, Dim0, Dim2, Dim3, j, m, i, k, l>;
    return Tensor2_Expr<TensorExpr, typename promote<T, U>::V, Dim1, Dim4, j, m>( TensorExpr( a, b ) );
  }

  // B(i, k, l, )*A(i, j, k, l, m, )
  template <class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, i, j, k, l, m>

  auto operator*( const Tensor3_Expr<B, U, Dim0, Dim2, Dim3, i, k, l> &b,
      const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a )
  {
    return a * b;
  }

  // A(i, j, k, l, m, )*B(m, k, i, )
  template <class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, i, j, k, l, m>

  auto operator*( const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
      const Tensor3_Expr<B, U, Dim4, Dim2, Dim0, m, k, i> &b )
  {
    using TensorExpr
        = Tensor5_times_Tensor3_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim4, Dim2, Dim0, i, j, k, l, m, m, k, i, Dim1, Dim3, Dim4, Dim2, Dim0, j, l, m, k, i>;
    return Tensor2_Expr<TensorExpr, typename promote<T, U>::V, Dim1, Dim3, j, l>( TensorExpr( a, b ) );
  }

  // B(m, k, i, )*A(i, j, k, l, m, )
  template <class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, i, j, k, l, m>

  auto operator*( const Tensor3_Expr<B, U, Dim4, Dim2, Dim0, m, k, i> &b,
      const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a )
  {
    return a * b;
  }

  // A(i, j, k, l, m, )*B(k, m, i, )
  template <class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, i, j, k, l, m>

  auto operator*( const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
      const Tensor3_Expr<B, U, Dim2, Dim4, Dim0, k, m, i> &b )
  {
    using TensorExpr
        = Tensor5_times_Tensor3_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim2, Dim4, Dim0, i, j, k, l, m, k, m, i, Dim1, Dim3, Dim2, Dim4, Dim0, j, l, k, m, i>;
    return Tensor2_Expr<TensorExpr, typename promote<T, U>::V, Dim1, Dim3, j, l>( TensorExpr( a, b ) );
  }

  // B(k, m, i, )*A(i, j, k, l, m, )
  template <class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, i, j, k, l, m>

  auto operator*( const Tensor3_Expr<B, U, Dim2, Dim4, Dim0, k, m, i> &b,
      const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a )
  {
    return a * b;
  }

  // A(i, j, k, l, m, )*B(m, i, k, )
  template <class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, i, j, k, l, m>

  auto operator*( const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
      const Tensor3_Expr<B, U, Dim4, Dim0, Dim2, m, i, k> &b )
  {
    using TensorExpr
        = Tensor5_times_Tensor3_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim4, Dim0, Dim2, i, j, k, l, m, m, i, k, Dim1, Dim3, Dim4, Dim0, Dim2, j, l, m, i, k>;
    return Tensor2_Expr<TensorExpr, typename promote<T, U>::V, Dim1, Dim3, j, l>( TensorExpr( a, b ) );
  }

  // B(m, i, k, )*A(i, j, k, l, m, )
  template <class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, i, j, k, l, m>

  auto operator*( const Tensor3_Expr<B, U, Dim4, Dim0, Dim2, m, i, k> &b,
      const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a )
  {
    return a * b;
  }

  // A(i, j, k, l, m, )*B(k, i, m, )
  template <class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, i, j, k, l, m>

  auto operator*( const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
      const Tensor3_Expr<B, U, Dim2, Dim0, Dim4, k, i, m> &b )
  {
    using TensorExpr
        = Tensor5_times_Tensor3_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim2, Dim0, Dim4, i, j, k, l, m, k, i, m, Dim1, Dim3, Dim2, Dim0, Dim4, j, l, k, i, m>;
    return Tensor2_Expr<TensorExpr, typename promote<T, U>::V, Dim1, Dim3, j, l>( TensorExpr( a, b ) );
  }

  // B(k, i, m, )*A(i, j, k, l, m, )
  template <class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, i, j, k, l, m>

  auto operator*( const Tensor3_Expr<B, U, Dim2, Dim0, Dim4, k, i, m> &b,
      const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a )
  {
    return a * b;
  }

  // A(i, j, k, l, m, )*B(i, m, k, )
  template <class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, i, j, k, l, m>

  auto operator*( const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
      const Tensor3_Expr<B, U, Dim0, Dim4, Dim2, i, m, k> &b )
  {
    using TensorExpr
        = Tensor5_times_Tensor3_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim0, Dim4, Dim2, i, j, k, l, m, i, m, k, Dim1, Dim3, Dim0, Dim4, Dim2, j, l, i, m, k>;
    return Tensor2_Expr<TensorExpr, typename promote<T, U>::V, Dim1, Dim3, j, l>( TensorExpr( a, b ) );
  }

  // B(i, m, k, )*A(i, j, k, l, m, )
  template <class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, i, j, k, l, m>

  auto operator*( const Tensor3_Expr<B, U, Dim0, Dim4, Dim2, i, m, k> &b,
      const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a )
  {
    return a * b;
  }

  // A(i, j, k, l, m, )*B(i, k, m, )
  template <class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, i, j, k, l, m>

  auto operator*( const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
      const Tensor3_Expr<B, U, Dim0, Dim2, Dim4, i, k, m> &b )
  {
    using TensorExpr
        = Tensor5_times_Tensor3_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim0, Dim2, Dim4, i, j, k, l, m, i, k, m, Dim1, Dim3, Dim0, Dim2, Dim4, j, l, i, k, m>;
    return Tensor2_Expr<TensorExpr, typename promote<T, U>::V, Dim1, Dim3, j, l>( TensorExpr( a, b ) );
  }

  // B(i, k, m, )*A(i, j, k, l, m, )
  template <class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, i, j, k, l, m>

  auto operator*( const Tensor3_Expr<B, U, Dim0, Dim2, Dim4, i, k, m> &b,
      const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a )
  {
    return a * b;
  }

  // A(i, j, k, l, m, )*B(m, l, i, )
  template <class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, i, j, k, l, m>

  auto operator*( const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
      const Tensor3_Expr<B, U, Dim4, Dim3, Dim0, m, l, i> &b )
  {
    using TensorExpr
        = Tensor5_times_Tensor3_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim4, Dim3, Dim0, i, j, k, l, m, m, l, i, Dim1, Dim2, Dim4, Dim3, Dim0, j, k, m, l, i>;
    return Tensor2_Expr<TensorExpr, typename promote<T, U>::V, Dim1, Dim2, j, k>( TensorExpr( a, b ) );
  }

  // B(m, l, i, )*A(i, j, k, l, m, )
  template <class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, i, j, k, l, m>

  auto operator*( const Tensor3_Expr<B, U, Dim4, Dim3, Dim0, m, l, i> &b,
      const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a )
  {
    return a * b;
  }

  // A(i, j, k, l, m, )*B(l, m, i, )
  template <class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, i, j, k, l, m>

  auto operator*( const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
      const Tensor3_Expr<B, U, Dim3, Dim4, Dim0, l, m, i> &b )
  {
    using TensorExpr
        = Tensor5_times_Tensor3_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim3, Dim4, Dim0, i, j, k, l, m, l, m, i, Dim1, Dim2, Dim3, Dim4, Dim0, j, k, l, m, i>;
    return Tensor2_Expr<TensorExpr, typename promote<T, U>::V, Dim1, Dim2, j, k>( TensorExpr( a, b ) );
  }

  // B(l, m, i, )*A(i, j, k, l, m, )
  template <class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, i, j, k, l, m>

  auto operator*( const Tensor3_Expr<B, U, Dim3, Dim4, Dim0, l, m, i> &b,
      const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a )
  {
    return a * b;
  }

  // A(i, j, k, l, m, )*B(m, i, l, )
  template <class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, i, j, k, l, m>

  auto operator*( const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
      const Tensor3_Expr<B, U, Dim4, Dim0, Dim3, m, i, l> &b )
  {
    using TensorExpr
        = Tensor5_times_Tensor3_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim4, Dim0, Dim3, i, j, k, l, m, m, i, l, Dim1, Dim2, Dim4, Dim0, Dim3, j, k, m, i, l>;
    return Tensor2_Expr<TensorExpr, typename promote<T, U>::V, Dim1, Dim2, j, k>( TensorExpr( a, b ) );
  }

  // B(m, i, l, )*A(i, j, k, l, m, )
  template <class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, i, j, k, l, m>

  auto operator*( const Tensor3_Expr<B, U, Dim4, Dim0, Dim3, m, i, l> &b,
      const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a )
  {
    return a * b;
  }

  // A(i, j, k, l, m, )*B(l, i, m, )
  template <class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, i, j, k, l, m>

  auto operator*( const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
      const Tensor3_Expr<B, U, Dim3, Dim0, Dim4, l, i, m> &b )
  {
    using TensorExpr
        = Tensor5_times_Tensor3_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim3, Dim0, Dim4, i, j, k, l, m, l, i, m, Dim1, Dim2, Dim3, Dim0, Dim4, j, k, l, i, m>;
    return Tensor2_Expr<TensorExpr, typename promote<T, U>::V, Dim1, Dim2, j, k>( TensorExpr( a, b ) );
  }

  // B(l, i, m, )*A(i, j, k, l, m, )
  template <class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, i, j, k, l, m>

  auto operator*( const Tensor3_Expr<B, U, Dim3, Dim0, Dim4, l, i, m> &b,
      const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a )
  {
    return a * b;
  }

  // A(i, j, k, l, m, )*B(i, m, l, )
  template <class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, i, j, k, l, m>

  auto operator*( const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
      const Tensor3_Expr<B, U, Dim0, Dim4, Dim3, i, m, l> &b )
  {
    using TensorExpr
        = Tensor5_times_Tensor3_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim0, Dim4, Dim3, i, j, k, l, m, i, m, l, Dim1, Dim2, Dim0, Dim4, Dim3, j, k, i, m, l>;
    return Tensor2_Expr<TensorExpr, typename promote<T, U>::V, Dim1, Dim2, j, k>( TensorExpr( a, b ) );
  }

  // B(i, m, l, )*A(i, j, k, l, m, )
  template <class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, i, j, k, l, m>

  auto operator*( const Tensor3_Expr<B, U, Dim0, Dim4, Dim3, i, m, l> &b,
      const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a )
  {
    return a * b;
  }

  // A(i, j, k, l, m, )*B(i, l, m, )
  template <class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, i, j, k, l, m>

  auto operator*( const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
      const Tensor3_Expr<B, U, Dim0, Dim3, Dim4, i, l, m> &b )
  {
    using TensorExpr
        = Tensor5_times_Tensor3_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim0, Dim3, Dim4, i, j, k, l, m, i, l, m, Dim1, Dim2, Dim0, Dim3, Dim4, j, k, i, l, m>;
    return Tensor2_Expr<TensorExpr, typename promote<T, U>::V, Dim1, Dim2, j, k>( TensorExpr( a, b ) );
  }

  // B(i, l, m, )*A(i, j, k, l, m, )
  template <class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, i, j, k, l, m>

  auto operator*( const Tensor3_Expr<B, U, Dim0, Dim3, Dim4, i, l, m> &b,
      const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a )
  {
    return a * b;
  }

  // A(i, j, k, l, m, )*B(l, k, j, )
  template <class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, i, j, k, l, m>

  auto operator*( const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
      const Tensor3_Expr<B, U, Dim3, Dim2, Dim1, l, k, j> &b )
  {
    using TensorExpr
        = Tensor5_times_Tensor3_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim3, Dim2, Dim1, i, j, k, l, m, l, k, j, Dim0, Dim4, Dim3, Dim2, Dim1, i, m, l, k, j>;
    return Tensor2_Expr<TensorExpr, typename promote<T, U>::V, Dim0, Dim4, i, m>( TensorExpr( a, b ) );
  }

  // B(l, k, j, )*A(i, j, k, l, m, )
  template <class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, i, j, k, l, m>

  auto operator*( const Tensor3_Expr<B, U, Dim3, Dim2, Dim1, l, k, j> &b,
      const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a )
  {
    return a * b;
  }

  // A(i, j, k, l, m, )*B(k, l, j, )
  template <class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, i, j, k, l, m>

  auto operator*( const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
      const Tensor3_Expr<B, U, Dim2, Dim3, Dim1, k, l, j> &b )
  {
    using TensorExpr
        = Tensor5_times_Tensor3_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim2, Dim3, Dim1, i, j, k, l, m, k, l, j, Dim0, Dim4, Dim2, Dim3, Dim1, i, m, k, l, j>;
    return Tensor2_Expr<TensorExpr, typename promote<T, U>::V, Dim0, Dim4, i, m>( TensorExpr( a, b ) );
  }

  // B(k, l, j, )*A(i, j, k, l, m, )
  template <class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, i, j, k, l, m>

  auto operator*( const Tensor3_Expr<B, U, Dim2, Dim3, Dim1, k, l, j> &b,
      const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a )
  {
    return a * b;
  }

  // A(i, j, k, l, m, )*B(l, j, k, )
  template <class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, i, j, k, l, m>

  auto operator*( const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
      const Tensor3_Expr<B, U, Dim3, Dim1, Dim2, l, j, k> &b )
  {
    using TensorExpr
        = Tensor5_times_Tensor3_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim3, Dim1, Dim2, i, j, k, l, m, l, j, k, Dim0, Dim4, Dim3, Dim1, Dim2, i, m, l, j, k>;
    return Tensor2_Expr<TensorExpr, typename promote<T, U>::V, Dim0, Dim4, i, m>( TensorExpr( a, b ) );
  }

  // B(l, j, k, )*A(i, j, k, l, m, )
  template <class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, i, j, k, l, m>

  auto operator*( const Tensor3_Expr<B, U, Dim3, Dim1, Dim2, l, j, k> &b,
      const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a )
  {
    return a * b;
  }

  // A(i, j, k, l, m, )*B(k, j, l, )
  template <class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, i, j, k, l, m>

  auto operator*( const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
      const Tensor3_Expr<B, U, Dim2, Dim1, Dim3, k, j, l> &b )
  {
    using TensorExpr
        = Tensor5_times_Tensor3_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim2, Dim1, Dim3, i, j, k, l, m, k, j, l, Dim0, Dim4, Dim2, Dim1, Dim3, i, m, k, j, l>;
    return Tensor2_Expr<TensorExpr, typename promote<T, U>::V, Dim0, Dim4, i, m>( TensorExpr( a, b ) );
  }

  // B(k, j, l, )*A(i, j, k, l, m, )
  template <class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, i, j, k, l, m>

  auto operator*( const Tensor3_Expr<B, U, Dim2, Dim1, Dim3, k, j, l> &b,
      const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a )
  {
    return a * b;
  }

  // A(i, j, k, l, m, )*B(j, l, k, )
  template <class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, i, j, k, l, m>

  auto operator*( const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
      const Tensor3_Expr<B, U, Dim1, Dim3, Dim2, j, l, k> &b )
  {
    using TensorExpr
        = Tensor5_times_Tensor3_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim1, Dim3, Dim2, i, j, k, l, m, j, l, k, Dim0, Dim4, Dim1, Dim3, Dim2, i, m, j, l, k>;
    return Tensor2_Expr<TensorExpr, typename promote<T, U>::V, Dim0, Dim4, i, m>( TensorExpr( a, b ) );
  }

  // B(j, l, k, )*A(i, j, k, l, m, )
  template <class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, i, j, k, l, m>

  auto operator*( const Tensor3_Expr<B, U, Dim1, Dim3, Dim2, j, l, k> &b,
      const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a )
  {
    return a * b;
  }

  // A(i, j, k, l, m, )*B(j, k, l, )
  template <class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, i, j, k, l, m>

  auto operator*( const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
      const Tensor3_Expr<B, U, Dim1, Dim2, Dim3, j, k, l> &b )
  {
    using TensorExpr
        = Tensor5_times_Tensor3_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim1, Dim2, Dim3, i, j, k, l, m, j, k, l, Dim0, Dim4, Dim1, Dim2, Dim3, i, m, j, k, l>;
    return Tensor2_Expr<TensorExpr, typename promote<T, U>::V, Dim0, Dim4, i, m>( TensorExpr( a, b ) );
  }

  // B(j, k, l, )*A(i, j, k, l, m, )
  template <class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, i, j, k, l, m>

  auto operator*( const Tensor3_Expr<B, U, Dim1, Dim2, Dim3, j, k, l> &b,
      const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a )
  {
    return a * b;
  }

  // A(i, j, k, l, m, )*B(m, k, j, )
  template <class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, i, j, k, l, m>

  auto operator*( const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
      const Tensor3_Expr<B, U, Dim4, Dim2, Dim1, m, k, j> &b )
  {
    using TensorExpr
        = Tensor5_times_Tensor3_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim4, Dim2, Dim1, i, j, k, l, m, m, k, j, Dim0, Dim3, Dim4, Dim2, Dim1, i, l, m, k, j>;
    return Tensor2_Expr<TensorExpr, typename promote<T, U>::V, Dim0, Dim3, i, l>( TensorExpr( a, b ) );
  }

  // B(m, k, j, )*A(i, j, k, l, m, )
  template <class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, i, j, k, l, m>

  auto operator*( const Tensor3_Expr<B, U, Dim4, Dim2, Dim1, m, k, j> &b,
      const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a )
  {
    return a * b;
  }

  // A(i, j, k, l, m, )*B(k, m, j, )
  template <class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, i, j, k, l, m>

  auto operator*( const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
      const Tensor3_Expr<B, U, Dim2, Dim4, Dim1, k, m, j> &b )
  {
    using TensorExpr
        = Tensor5_times_Tensor3_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim2, Dim4, Dim1, i, j, k, l, m, k, m, j, Dim0, Dim3, Dim2, Dim4, Dim1, i, l, k, m, j>;
    return Tensor2_Expr<TensorExpr, typename promote<T, U>::V, Dim0, Dim3, i, l>( TensorExpr( a, b ) );
  }

  // B(k, m, j, )*A(i, j, k, l, m, )
  template <class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, i, j, k, l, m>

  auto operator*( const Tensor3_Expr<B, U, Dim2, Dim4, Dim1, k, m, j> &b,
      const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a )
  {
    return a * b;
  }

  // A(i, j, k, l, m, )*B(m, j, k, )
  template <class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, i, j, k, l, m>

  auto operator*( const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
      const Tensor3_Expr<B, U, Dim4, Dim1, Dim2, m, j, k> &b )
  {
    using TensorExpr
        = Tensor5_times_Tensor3_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim4, Dim1, Dim2, i, j, k, l, m, m, j, k, Dim0, Dim3, Dim4, Dim1, Dim2, i, l, m, j, k>;
    return Tensor2_Expr<TensorExpr, typename promote<T, U>::V, Dim0, Dim3, i, l>( TensorExpr( a, b ) );
  }

  // B(m, j, k, )*A(i, j, k, l, m, )
  template <class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, i, j, k, l, m>

  auto operator*( const Tensor3_Expr<B, U, Dim4, Dim1, Dim2, m, j, k> &b,
      const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a )
  {
    return a * b;
  }

  // A(i, j, k, l, m, )*B(k, j, m, )
  template <class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, i, j, k, l, m>

  auto operator*( const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
      const Tensor3_Expr<B, U, Dim2, Dim1, Dim4, k, j, m> &b )
  {
    using TensorExpr
        = Tensor5_times_Tensor3_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim2, Dim1, Dim4, i, j, k, l, m, k, j, m, Dim0, Dim3, Dim2, Dim1, Dim4, i, l, k, j, m>;
    return Tensor2_Expr<TensorExpr, typename promote<T, U>::V, Dim0, Dim3, i, l>( TensorExpr( a, b ) );
  }

  // B(k, j, m, )*A(i, j, k, l, m, )
  template <class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, i, j, k, l, m>

  auto operator*( const Tensor3_Expr<B, U, Dim2, Dim1, Dim4, k, j, m> &b,
      const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a )
  {
    return a * b;
  }

  // A(i, j, k, l, m, )*B(j, m, k, )
  template <class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, i, j, k, l, m>

  auto operator*( const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
      const Tensor3_Expr<B, U, Dim1, Dim4, Dim2, j, m, k> &b )
  {
    using TensorExpr
        = Tensor5_times_Tensor3_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim1, Dim4, Dim2, i, j, k, l, m, j, m, k, Dim0, Dim3, Dim1, Dim4, Dim2, i, l, j, m, k>;
    return Tensor2_Expr<TensorExpr, typename promote<T, U>::V, Dim0, Dim3, i, l>( TensorExpr( a, b ) );
  }

  // B(j, m, k, )*A(i, j, k, l, m, )
  template <class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, i, j, k, l, m>

  auto operator*( const Tensor3_Expr<B, U, Dim1, Dim4, Dim2, j, m, k> &b,
      const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a )
  {
    return a * b;
  }

  // A(i, j, k, l, m, )*B(j, k, m, )
  template <class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, i, j, k, l, m>

  auto operator*( const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
      const Tensor3_Expr<B, U, Dim1, Dim2, Dim4, j, k, m> &b )
  {
    using TensorExpr
        = Tensor5_times_Tensor3_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim1, Dim2, Dim4, i, j, k, l, m, j, k, m, Dim0, Dim3, Dim1, Dim2, Dim4, i, l, j, k, m>;
    return Tensor2_Expr<TensorExpr, typename promote<T, U>::V, Dim0, Dim3, i, l>( TensorExpr( a, b ) );
  }

  // B(j, k, m, )*A(i, j, k, l, m, )
  template <class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, i, j, k, l, m>

  auto operator*( const Tensor3_Expr<B, U, Dim1, Dim2, Dim4, j, k, m> &b,
      const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a )
  {
    return a * b;
  }

  // A(i, j, k, l, m, )*B(m, l, j, )
  template <class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, i, j, k, l, m>

  auto operator*( const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
      const Tensor3_Expr<B, U, Dim4, Dim3, Dim1, m, l, j> &b )
  {
    using TensorExpr
        = Tensor5_times_Tensor3_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim4, Dim3, Dim1, i, j, k, l, m, m, l, j, Dim0, Dim2, Dim4, Dim3, Dim1, i, k, m, l, j>;
    return Tensor2_Expr<TensorExpr, typename promote<T, U>::V, Dim0, Dim2, i, k>( TensorExpr( a, b ) );
  }

  // B(m, l, j, )*A(i, j, k, l, m, )
  template <class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, i, j, k, l, m>

  auto operator*( const Tensor3_Expr<B, U, Dim4, Dim3, Dim1, m, l, j> &b,
      const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a )
  {
    return a * b;
  }

  // A(i, j, k, l, m, )*B(l, m, j, )
  template <class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, i, j, k, l, m>

  auto operator*( const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
      const Tensor3_Expr<B, U, Dim3, Dim4, Dim1, l, m, j> &b )
  {
    using TensorExpr
        = Tensor5_times_Tensor3_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim3, Dim4, Dim1, i, j, k, l, m, l, m, j, Dim0, Dim2, Dim3, Dim4, Dim1, i, k, l, m, j>;
    return Tensor2_Expr<TensorExpr, typename promote<T, U>::V, Dim0, Dim2, i, k>( TensorExpr( a, b ) );
  }

  // B(l, m, j, )*A(i, j, k, l, m, )
  template <class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, i, j, k, l, m>

  auto operator*( const Tensor3_Expr<B, U, Dim3, Dim4, Dim1, l, m, j> &b,
      const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a )
  {
    return a * b;
  }

  // A(i, j, k, l, m, )*B(m, j, l, )
  template <class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, i, j, k, l, m>

  auto operator*( const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
      const Tensor3_Expr<B, U, Dim4, Dim1, Dim3, m, j, l> &b )
  {
    using TensorExpr
        = Tensor5_times_Tensor3_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim4, Dim1, Dim3, i, j, k, l, m, m, j, l, Dim0, Dim2, Dim4, Dim1, Dim3, i, k, m, j, l>;
    return Tensor2_Expr<TensorExpr, typename promote<T, U>::V, Dim0, Dim2, i, k>( TensorExpr( a, b ) );
  }

  // B(m, j, l, )*A(i, j, k, l, m, )
  template <class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, i, j, k, l, m>

  auto operator*( const Tensor3_Expr<B, U, Dim4, Dim1, Dim3, m, j, l> &b,
      const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a )
  {
    return a * b;
  }

  // A(i, j, k, l, m, )*B(l, j, m, )
  template <class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, i, j, k, l, m>

  auto operator*( const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
      const Tensor3_Expr<B, U, Dim3, Dim1, Dim4, l, j, m> &b )
  {
    using TensorExpr
        = Tensor5_times_Tensor3_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim3, Dim1, Dim4, i, j, k, l, m, l, j, m, Dim0, Dim2, Dim3, Dim1, Dim4, i, k, l, j, m>;
    return Tensor2_Expr<TensorExpr, typename promote<T, U>::V, Dim0, Dim2, i, k>( TensorExpr( a, b ) );
  }

  // B(l, j, m, )*A(i, j, k, l, m, )
  template <class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, i, j, k, l, m>

  auto operator*( const Tensor3_Expr<B, U, Dim3, Dim1, Dim4, l, j, m> &b,
      const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a )
  {
    return a * b;
  }

  // A(i, j, k, l, m, )*B(j, m, l, )
  template <class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, i, j, k, l, m>

  auto operator*( const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
      const Tensor3_Expr<B, U, Dim1, Dim4, Dim3, j, m, l> &b )
  {
    using TensorExpr
        = Tensor5_times_Tensor3_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim1, Dim4, Dim3, i, j, k, l, m, j, m, l, Dim0, Dim2, Dim1, Dim4, Dim3, i, k, j, m, l>;
    return Tensor2_Expr<TensorExpr, typename promote<T, U>::V, Dim0, Dim2, i, k>( TensorExpr( a, b ) );
  }

  // B(j, m, l, )*A(i, j, k, l, m, )
  template <class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, i, j, k, l, m>

  auto operator*( const Tensor3_Expr<B, U, Dim1, Dim4, Dim3, j, m, l> &b,
      const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a )
  {
    return a * b;
  }

  // A(i, j, k, l, m, )*B(j, l, m, )
  template <class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, i, j, k, l, m>

  auto operator*( const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
      const Tensor3_Expr<B, U, Dim1, Dim3, Dim4, j, l, m> &b )
  {
    using TensorExpr
        = Tensor5_times_Tensor3_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim1, Dim3, Dim4, i, j, k, l, m, j, l, m, Dim0, Dim2, Dim1, Dim3, Dim4, i, k, j, l, m>;
    return Tensor2_Expr<TensorExpr, typename promote<T, U>::V, Dim0, Dim2, i, k>( TensorExpr( a, b ) );
  }

  // B(j, l, m, )*A(i, j, k, l, m, )
  template <class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, i, j, k, l, m>

  auto operator*( const Tensor3_Expr<B, U, Dim1, Dim3, Dim4, j, l, m> &b,
      const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a )
  {
    return a * b;
  }

  // A(i, j, k, l, m, )*B(m, l, k, )
  template <class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, i, j, k, l, m>

  auto operator*( const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
      const Tensor3_Expr<B, U, Dim4, Dim3, Dim2, m, l, k> &b )
  {
    using TensorExpr
        = Tensor5_times_Tensor3_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim4, Dim3, Dim2, i, j, k, l, m, m, l, k, Dim0, Dim1, Dim4, Dim3, Dim2, i, j, m, l, k>;
    return Tensor2_Expr<TensorExpr, typename promote<T, U>::V, Dim0, Dim1, i, j>( TensorExpr( a, b ) );
  }

  // B(m, l, k, )*A(i, j, k, l, m, )
  template <class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, i, j, k, l, m>

  auto operator*( const Tensor3_Expr<B, U, Dim4, Dim3, Dim2, m, l, k> &b,
      const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a )
  {
    return a * b;
  }

  // A(i, j, k, l, m, )*B(l, m, k, )
  template <class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, i, j, k, l, m>

  auto operator*( const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
      const Tensor3_Expr<B, U, Dim3, Dim4, Dim2, l, m, k> &b )
  {
    using TensorExpr
        = Tensor5_times_Tensor3_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim3, Dim4, Dim2, i, j, k, l, m, l, m, k, Dim0, Dim1, Dim3, Dim4, Dim2, i, j, l, m, k>;
    return Tensor2_Expr<TensorExpr, typename promote<T, U>::V, Dim0, Dim1, i, j>( TensorExpr( a, b ) );
  }

  // B(l, m, k, )*A(i, j, k, l, m, )
  template <class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, i, j, k, l, m>

  auto operator*( const Tensor3_Expr<B, U, Dim3, Dim4, Dim2, l, m, k> &b,
      const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a )
  {
    return a * b;
  }

  // A(i, j, k, l, m, )*B(m, k, l, )
  template <class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, i, j, k, l, m>

  auto operator*( const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
      const Tensor3_Expr<B, U, Dim4, Dim2, Dim3, m, k, l> &b )
  {
    using TensorExpr
        = Tensor5_times_Tensor3_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim4, Dim2, Dim3, i, j, k, l, m, m, k, l, Dim0, Dim1, Dim4, Dim2, Dim3, i, j, m, k, l>;
    return Tensor2_Expr<TensorExpr, typename promote<T, U>::V, Dim0, Dim1, i, j>( TensorExpr( a, b ) );
  }

  // B(m, k, l, )*A(i, j, k, l, m, )
  template <class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, i, j, k, l, m>

  auto operator*( const Tensor3_Expr<B, U, Dim4, Dim2, Dim3, m, k, l> &b,
      const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a )
  {
    return a * b;
  }

  // A(i, j, k, l, m, )*B(l, k, m, )
  template <class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, i, j, k, l, m>

  auto operator*( const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
      const Tensor3_Expr<B, U, Dim3, Dim2, Dim4, l, k, m> &b )
  {
    using TensorExpr
        = Tensor5_times_Tensor3_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim3, Dim2, Dim4, i, j, k, l, m, l, k, m, Dim0, Dim1, Dim3, Dim2, Dim4, i, j, l, k, m>;
    return Tensor2_Expr<TensorExpr, typename promote<T, U>::V, Dim0, Dim1, i, j>( TensorExpr( a, b ) );
  }

  // B(l, k, m, )*A(i, j, k, l, m, )
  template <class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, i, j, k, l, m>

  auto operator*( const Tensor3_Expr<B, U, Dim3, Dim2, Dim4, l, k, m> &b,
      const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a )
  {
    return a * b;
  }

  // A(i, j, k, l, m, )*B(k, m, l, )
  template <class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, i, j, k, l, m>

  auto operator*( const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
      const Tensor3_Expr<B, U, Dim2, Dim4, Dim3, k, m, l> &b )
  {
    using TensorExpr
        = Tensor5_times_Tensor3_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim2, Dim4, Dim3, i, j, k, l, m, k, m, l, Dim0, Dim1, Dim2, Dim4, Dim3, i, j, k, m, l>;
    return Tensor2_Expr<TensorExpr, typename promote<T, U>::V, Dim0, Dim1, i, j>( TensorExpr( a, b ) );
  }

  // B(k, m, l, )*A(i, j, k, l, m, )
  template <class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, i, j, k, l, m>

  auto operator*( const Tensor3_Expr<B, U, Dim2, Dim4, Dim3, k, m, l> &b,
      const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a )
  {
    return a * b;
  }

  // A(i, j, k, l, m, )*B(k, l, m, )
  template <class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, i, j, k, l, m>

  auto operator*( const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
      const Tensor3_Expr<B, U, Dim2, Dim3, Dim4, k, l, m> &b )
  {
    using TensorExpr
        = Tensor5_times_Tensor3_triple<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim2, Dim3, Dim4, i, j, k, l, m, k, l, m, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m>;
    return Tensor2_Expr<TensorExpr, typename promote<T, U>::V, Dim0, Dim1, i, j>( TensorExpr( a, b ) );
  }

  // B(k, l, m, )*A(i, j, k, l, m, )
  template <class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, i, j, k, l, m>

  auto operator*( const Tensor3_Expr<B, U, Dim2, Dim3, Dim4, k, l, m> &b,
      const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a )
  {
    return a * b;
  }


  }
