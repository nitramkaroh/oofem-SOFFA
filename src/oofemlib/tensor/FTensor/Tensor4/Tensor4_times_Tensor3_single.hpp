#pragma once

namespace FTensor
{
//Universal A(i, j, k, l, )*B(m, n, o, ) single contraction
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o, int DimA, int DimB, int DimC, int DimX, int DimY, int DimZ, char a, char b, char c, char x, char y, char z>
class Tensor4_times_Tensor3_single
{
  Tensor4_Expr<A, T, Dim0, Dim1, Dim2, Dim3, i, j, k, l> iterA;
  Tensor3_Expr<B, U, Dim4, Dim5, Dim6, m, n, o> iterB;

public:
  typename promote<T, U>::V operator()(const int N1, const int N2, const int N3, const int N4, const int N5) const
  {
    typename promote<T, U>::V result(0);
    for(int xx = 0; xx < DimX; ++xx)
      {
        // Permutation is where the indices get checked.
        result += Permutation4<DimA, DimB, DimC, DimX, a, b, c, x>().eval(
          iterA, N1, N2, N3, xx) * Permutation3<DimX, DimY, DimZ, x, y, z>().eval(iterB, xx, N4, N5);
      }
    return result;
  }

  Tensor4_times_Tensor3_single(
    const Tensor4_Expr<A, T, Dim0, Dim1, Dim2, Dim3, i, j, k, l> &iter_a,
    const Tensor3_Expr<B, U, Dim4, Dim5, Dim6, m, n, o> &iter_b)
      : iterA(iter_a), iterB(iter_b)
  {}
};

//A(i, j, k, l, )*B(i, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<A, T, Dim0, Dim1, Dim2, Dim3, i, j, k, l> &a,
               const Tensor3_Expr<B, U, Dim0, Dim4, Dim5, i, m, n> &b)
{
  using TensorExpr
    = Tensor4_times_Tensor3_single<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim0, Dim4, Dim5, i, j, k, l, i, m, n, Dim1, Dim2, Dim3, Dim0, Dim4, Dim5, j, k, l, i, m, n>;
  return Tensor5_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim2, Dim3, Dim4, Dim5, j, k, l, m, n>(TensorExpr(a, b));
}

//B(i, m, n, )*A(i, j, k, l, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor3_Expr<B, U, Dim0, Dim4, Dim5, i, m, n> &b,
               const Tensor4_Expr<A, T, Dim0, Dim1, Dim2, Dim3, i, j, k, l> &a)
{
  return a * b;
}

//A(i, j, k, l, )*B(m, i, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<A, T, Dim0, Dim1, Dim2, Dim3, i, j, k, l> &a,
               const Tensor3_Expr<B, U, Dim4, Dim0, Dim5, m, i, n> &b)
{
  using TensorExpr
    = Tensor4_times_Tensor3_single<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim0, Dim5, i, j, k, l, m, i, n, Dim1, Dim2, Dim3, Dim4, Dim0, Dim5, j, k, l, m, i, n>;
  return Tensor5_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim2, Dim3, Dim4, Dim5, j, k, l, m, n>(TensorExpr(a, b));
}

//B(m, i, n, )*A(i, j, k, l, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor3_Expr<B, U, Dim4, Dim0, Dim5, m, i, n> &b,
               const Tensor4_Expr<A, T, Dim0, Dim1, Dim2, Dim3, i, j, k, l> &a)
{
  return a * b;
}

//A(i, j, k, l, )*B(m, n, i, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<A, T, Dim0, Dim1, Dim2, Dim3, i, j, k, l> &a,
               const Tensor3_Expr<B, U, Dim4, Dim5, Dim0, m, n, i> &b)
{
  using TensorExpr
    = Tensor4_times_Tensor3_single<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim0, i, j, k, l, m, n, i, Dim1, Dim2, Dim3, Dim4, Dim5, Dim0, j, k, l, m, n, i>;
  return Tensor5_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim2, Dim3, Dim4, Dim5, j, k, l, m, n>(TensorExpr(a, b));
}

//B(m, n, i, )*A(i, j, k, l, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor3_Expr<B, U, Dim4, Dim5, Dim0, m, n, i> &b,
               const Tensor4_Expr<A, T, Dim0, Dim1, Dim2, Dim3, i, j, k, l> &a)
{
  return a * b;
}

//A(i, j, k, l, )*B(j, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<A, T, Dim0, Dim1, Dim2, Dim3, i, j, k, l> &a,
               const Tensor3_Expr<B, U, Dim1, Dim4, Dim5, j, m, n> &b)
{
  using TensorExpr
    = Tensor4_times_Tensor3_single<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim1, Dim4, Dim5, i, j, k, l, j, m, n, Dim0, Dim2, Dim3, Dim1, Dim4, Dim5, i, k, l, j, m, n>;
  return Tensor5_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim2, Dim3, Dim4, Dim5, i, k, l, m, n>(TensorExpr(a, b));
}

//B(j, m, n, )*A(i, j, k, l, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor3_Expr<B, U, Dim1, Dim4, Dim5, j, m, n> &b,
               const Tensor4_Expr<A, T, Dim0, Dim1, Dim2, Dim3, i, j, k, l> &a)
{
  return a * b;
}

//A(i, j, k, l, )*B(m, j, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<A, T, Dim0, Dim1, Dim2, Dim3, i, j, k, l> &a,
               const Tensor3_Expr<B, U, Dim4, Dim1, Dim5, m, j, n> &b)
{
  using TensorExpr
    = Tensor4_times_Tensor3_single<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim1, Dim5, i, j, k, l, m, j, n, Dim0, Dim2, Dim3, Dim4, Dim1, Dim5, i, k, l, m, j, n>;
  return Tensor5_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim2, Dim3, Dim4, Dim5, i, k, l, m, n>(TensorExpr(a, b));
}

//B(m, j, n, )*A(i, j, k, l, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor3_Expr<B, U, Dim4, Dim1, Dim5, m, j, n> &b,
               const Tensor4_Expr<A, T, Dim0, Dim1, Dim2, Dim3, i, j, k, l> &a)
{
  return a * b;
}

//A(i, j, k, l, )*B(m, n, j, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<A, T, Dim0, Dim1, Dim2, Dim3, i, j, k, l> &a,
               const Tensor3_Expr<B, U, Dim4, Dim5, Dim1, m, n, j> &b)
{
  using TensorExpr
    = Tensor4_times_Tensor3_single<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim1, i, j, k, l, m, n, j, Dim0, Dim2, Dim3, Dim4, Dim5, Dim1, i, k, l, m, n, j>;
  return Tensor5_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim2, Dim3, Dim4, Dim5, i, k, l, m, n>(TensorExpr(a, b));
}

//B(m, n, j, )*A(i, j, k, l, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor3_Expr<B, U, Dim4, Dim5, Dim1, m, n, j> &b,
               const Tensor4_Expr<A, T, Dim0, Dim1, Dim2, Dim3, i, j, k, l> &a)
{
  return a * b;
}

//A(i, j, k, l, )*B(k, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<A, T, Dim0, Dim1, Dim2, Dim3, i, j, k, l> &a,
               const Tensor3_Expr<B, U, Dim2, Dim4, Dim5, k, m, n> &b)
{
  using TensorExpr
    = Tensor4_times_Tensor3_single<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim2, Dim4, Dim5, i, j, k, l, k, m, n, Dim0, Dim1, Dim3, Dim2, Dim4, Dim5, i, j, l, k, m, n>;
  return Tensor5_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim3, Dim4, Dim5, i, j, l, m, n>(TensorExpr(a, b));
}

//B(k, m, n, )*A(i, j, k, l, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor3_Expr<B, U, Dim2, Dim4, Dim5, k, m, n> &b,
               const Tensor4_Expr<A, T, Dim0, Dim1, Dim2, Dim3, i, j, k, l> &a)
{
  return a * b;
}

//A(i, j, k, l, )*B(m, k, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<A, T, Dim0, Dim1, Dim2, Dim3, i, j, k, l> &a,
               const Tensor3_Expr<B, U, Dim4, Dim2, Dim5, m, k, n> &b)
{
  using TensorExpr
    = Tensor4_times_Tensor3_single<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim2, Dim5, i, j, k, l, m, k, n, Dim0, Dim1, Dim3, Dim4, Dim2, Dim5, i, j, l, m, k, n>;
  return Tensor5_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim3, Dim4, Dim5, i, j, l, m, n>(TensorExpr(a, b));
}

//B(m, k, n, )*A(i, j, k, l, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor3_Expr<B, U, Dim4, Dim2, Dim5, m, k, n> &b,
               const Tensor4_Expr<A, T, Dim0, Dim1, Dim2, Dim3, i, j, k, l> &a)
{
  return a * b;
}

//A(i, j, k, l, )*B(m, n, k, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<A, T, Dim0, Dim1, Dim2, Dim3, i, j, k, l> &a,
               const Tensor3_Expr<B, U, Dim4, Dim5, Dim2, m, n, k> &b)
{
  using TensorExpr
    = Tensor4_times_Tensor3_single<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim2, i, j, k, l, m, n, k, Dim0, Dim1, Dim3, Dim4, Dim5, Dim2, i, j, l, m, n, k>;
  return Tensor5_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim3, Dim4, Dim5, i, j, l, m, n>(TensorExpr(a, b));
}

//B(m, n, k, )*A(i, j, k, l, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor3_Expr<B, U, Dim4, Dim5, Dim2, m, n, k> &b,
               const Tensor4_Expr<A, T, Dim0, Dim1, Dim2, Dim3, i, j, k, l> &a)
{
  return a * b;
}

//A(i, j, k, l, )*B(l, m, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<A, T, Dim0, Dim1, Dim2, Dim3, i, j, k, l> &a,
               const Tensor3_Expr<B, U, Dim3, Dim4, Dim5, l, m, n> &b)
{
  using TensorExpr
    = Tensor4_times_Tensor3_single<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim3, Dim4, Dim5, i, j, k, l, l, m, n, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, i, j, k, l, m, n>;
  return Tensor5_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim2, Dim4, Dim5, i, j, k, m, n>(TensorExpr(a, b));
}

//B(l, m, n, )*A(i, j, k, l, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor3_Expr<B, U, Dim3, Dim4, Dim5, l, m, n> &b,
               const Tensor4_Expr<A, T, Dim0, Dim1, Dim2, Dim3, i, j, k, l> &a)
{
  return a * b;
}

//A(i, j, k, l, )*B(m, l, n, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<A, T, Dim0, Dim1, Dim2, Dim3, i, j, k, l> &a,
               const Tensor3_Expr<B, U, Dim4, Dim3, Dim5, m, l, n> &b)
{
  using TensorExpr
    = Tensor4_times_Tensor3_single<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim3, Dim5, i, j, k, l, m, l, n, Dim0, Dim1, Dim2, Dim4, Dim3, Dim5, i, j, k, m, l, n>;
  return Tensor5_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim2, Dim4, Dim5, i, j, k, m, n>(TensorExpr(a, b));
}

//B(m, l, n, )*A(i, j, k, l, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor3_Expr<B, U, Dim4, Dim3, Dim5, m, l, n> &b,
               const Tensor4_Expr<A, T, Dim0, Dim1, Dim2, Dim3, i, j, k, l> &a)
{
  return a * b;
}

//A(i, j, k, l, )*B(m, n, l, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor4_Expr<A, T, Dim0, Dim1, Dim2, Dim3, i, j, k, l> &a,
               const Tensor3_Expr<B, U, Dim4, Dim5, Dim3, m, n, l> &b)
{
  using TensorExpr
    = Tensor4_times_Tensor3_single<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim3, i, j, k, l, m, n, l, Dim0, Dim1, Dim2, Dim4, Dim5, Dim3, i, j, k, m, n, l>;
  return Tensor5_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim2, Dim4, Dim5, i, j, k, m, n>(TensorExpr(a, b));
}

//B(m, n, l, )*A(i, j, k, l, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, char i, char j, char k, char l, char m, char n>
auto operator*(const Tensor3_Expr<B, U, Dim4, Dim5, Dim3, m, n, l> &b,
               const Tensor4_Expr<A, T, Dim0, Dim1, Dim2, Dim3, i, j, k, l> &a)
{
  return a * b;
}

} //namespace FTensor
