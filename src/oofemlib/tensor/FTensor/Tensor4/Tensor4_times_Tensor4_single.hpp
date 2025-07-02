#pragma once

namespace FTensor
{
//Universal A(i, j, k, l, )*B(m, n, o, p, ) single contraction
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p, int DimA, int DimB, int DimC, int DimX, int DimY, int DimZ, int DimW, char a, char b, char c, char x, char y, char z, char w>
class Tensor4_times_Tensor4_single
{
  Tensor4_Expr<A, T, Dim0, Dim1, Dim2, Dim3, i, j, k, l> iterA;
  Tensor4_Expr<B, U, Dim4, Dim5, Dim6, Dim7, m, n, o, p> iterB;

public:
  typename promote<T, U>::V operator()(const int N1, const int N2, const int N3, const int N4, const int N5, const int N6) const
  {
    typename promote<T, U>::V result(0);
    for(int xx = 0; xx < DimX, ++xx)
      {
        // Permutation is where the indices get checked.
        result += Permutation4<DimA, DimB, DimC, DimX, a, b, c, x>().eval(
          iterA, N1, N2, N3, xx) * Permutation4<DimX, DimY, DimZ, DimW, x, y, z, w>().eval(iterB, xx, N4, N5, N6);
      }
    return result;
  }

  Tensor4_times_Tensor4_single(
    const Tensor4_Expr<A, T, Dim0, Dim1, Dim2, Dim3, i, j, k, l> &iter_a,
    const Tensor4_Expr<B, U, Dim4, Dim5, Dim6, Dim7, m, n, o, p> &iter_b)
      : iterA(iter_a), iterB(iter_b)
  {}
};

//A(i, j, k, l, )*B(i, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<A, T, Dim0, Dim1, Dim2, Dim3, i, j, k, l> &a,
               const Tensor4_Expr<B, U, Dim0, Dim4, Dim5, Dim6, i, m, n, o> &b)
{
  using TensorExpr
    = Tensor4_times_Tensor4_single<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim0, Dim4, Dim5, Dim6, i, j, k, l, i, m, n, o, Dim1, Dim2, Dim3, Dim0, Dim4, Dim5, Dim6, j, k, l, i, m, n, o>;
  return Tensor6_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, j, k, l, m, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, )*B(m, i, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<A, T, Dim0, Dim1, Dim2, Dim3, i, j, k, l> &a,
               const Tensor4_Expr<B, U, Dim4, Dim0, Dim5, Dim6, m, i, n, o> &b)
{
  using TensorExpr
    = Tensor4_times_Tensor4_single<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim0, Dim5, Dim6, i, j, k, l, m, i, n, o, Dim1, Dim2, Dim3, Dim4, Dim0, Dim5, Dim6, j, k, l, m, i, n, o>;
  return Tensor6_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, j, k, l, m, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, )*B(m, n, i, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<A, T, Dim0, Dim1, Dim2, Dim3, i, j, k, l> &a,
               const Tensor4_Expr<B, U, Dim4, Dim5, Dim0, Dim6, m, n, i, o> &b)
{
  using TensorExpr
    = Tensor4_times_Tensor4_single<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim0, Dim6, i, j, k, l, m, n, i, o, Dim1, Dim2, Dim3, Dim4, Dim5, Dim0, Dim6, j, k, l, m, n, i, o>;
  return Tensor6_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, j, k, l, m, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, )*B(m, n, o, i, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<A, T, Dim0, Dim1, Dim2, Dim3, i, j, k, l> &a,
               const Tensor4_Expr<B, U, Dim4, Dim5, Dim6, Dim0, m, n, o, i> &b)
{
  using TensorExpr
    = Tensor4_times_Tensor4_single<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim0, i, j, k, l, m, n, o, i, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim0, j, k, l, m, n, o, i>;
  return Tensor6_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, j, k, l, m, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, )*B(j, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<A, T, Dim0, Dim1, Dim2, Dim3, i, j, k, l> &a,
               const Tensor4_Expr<B, U, Dim1, Dim4, Dim5, Dim6, j, m, n, o> &b)
{
  using TensorExpr
    = Tensor4_times_Tensor4_single<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim1, Dim4, Dim5, Dim6, i, j, k, l, j, m, n, o, Dim0, Dim2, Dim3, Dim1, Dim4, Dim5, Dim6, i, k, l, j, m, n, o>;
  return Tensor6_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim2, Dim3, Dim4, Dim5, Dim6, i, k, l, m, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, )*B(m, j, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<A, T, Dim0, Dim1, Dim2, Dim3, i, j, k, l> &a,
               const Tensor4_Expr<B, U, Dim4, Dim1, Dim5, Dim6, m, j, n, o> &b)
{
  using TensorExpr
    = Tensor4_times_Tensor4_single<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim1, Dim5, Dim6, i, j, k, l, m, j, n, o, Dim0, Dim2, Dim3, Dim4, Dim1, Dim5, Dim6, i, k, l, m, j, n, o>;
  return Tensor6_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim2, Dim3, Dim4, Dim5, Dim6, i, k, l, m, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, )*B(m, n, j, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<A, T, Dim0, Dim1, Dim2, Dim3, i, j, k, l> &a,
               const Tensor4_Expr<B, U, Dim4, Dim5, Dim1, Dim6, m, n, j, o> &b)
{
  using TensorExpr
    = Tensor4_times_Tensor4_single<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim1, Dim6, i, j, k, l, m, n, j, o, Dim0, Dim2, Dim3, Dim4, Dim5, Dim1, Dim6, i, k, l, m, n, j, o>;
  return Tensor6_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim2, Dim3, Dim4, Dim5, Dim6, i, k, l, m, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, )*B(m, n, o, j, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<A, T, Dim0, Dim1, Dim2, Dim3, i, j, k, l> &a,
               const Tensor4_Expr<B, U, Dim4, Dim5, Dim6, Dim1, m, n, o, j> &b)
{
  using TensorExpr
    = Tensor4_times_Tensor4_single<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim1, i, j, k, l, m, n, o, j, Dim0, Dim2, Dim3, Dim4, Dim5, Dim6, Dim1, i, k, l, m, n, o, j>;
  return Tensor6_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim2, Dim3, Dim4, Dim5, Dim6, i, k, l, m, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, )*B(k, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<A, T, Dim0, Dim1, Dim2, Dim3, i, j, k, l> &a,
               const Tensor4_Expr<B, U, Dim2, Dim4, Dim5, Dim6, k, m, n, o> &b)
{
  using TensorExpr
    = Tensor4_times_Tensor4_single<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim2, Dim4, Dim5, Dim6, i, j, k, l, k, m, n, o, Dim0, Dim1, Dim3, Dim2, Dim4, Dim5, Dim6, i, j, l, k, m, n, o>;
  return Tensor6_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim3, Dim4, Dim5, Dim6, i, j, l, m, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, )*B(m, k, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<A, T, Dim0, Dim1, Dim2, Dim3, i, j, k, l> &a,
               const Tensor4_Expr<B, U, Dim4, Dim2, Dim5, Dim6, m, k, n, o> &b)
{
  using TensorExpr
    = Tensor4_times_Tensor4_single<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim2, Dim5, Dim6, i, j, k, l, m, k, n, o, Dim0, Dim1, Dim3, Dim4, Dim2, Dim5, Dim6, i, j, l, m, k, n, o>;
  return Tensor6_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim3, Dim4, Dim5, Dim6, i, j, l, m, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, )*B(m, n, k, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<A, T, Dim0, Dim1, Dim2, Dim3, i, j, k, l> &a,
               const Tensor4_Expr<B, U, Dim4, Dim5, Dim2, Dim6, m, n, k, o> &b)
{
  using TensorExpr
    = Tensor4_times_Tensor4_single<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim2, Dim6, i, j, k, l, m, n, k, o, Dim0, Dim1, Dim3, Dim4, Dim5, Dim2, Dim6, i, j, l, m, n, k, o>;
  return Tensor6_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim3, Dim4, Dim5, Dim6, i, j, l, m, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, )*B(m, n, o, k, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<A, T, Dim0, Dim1, Dim2, Dim3, i, j, k, l> &a,
               const Tensor4_Expr<B, U, Dim4, Dim5, Dim6, Dim2, m, n, o, k> &b)
{
  using TensorExpr
    = Tensor4_times_Tensor4_single<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim2, i, j, k, l, m, n, o, k, Dim0, Dim1, Dim3, Dim4, Dim5, Dim6, Dim2, i, j, l, m, n, o, k>;
  return Tensor6_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim3, Dim4, Dim5, Dim6, i, j, l, m, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, )*B(l, m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<A, T, Dim0, Dim1, Dim2, Dim3, i, j, k, l> &a,
               const Tensor4_Expr<B, U, Dim3, Dim4, Dim5, Dim6, l, m, n, o> &b)
{
  using TensorExpr
    = Tensor4_times_Tensor4_single<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim3, Dim4, Dim5, Dim6, i, j, k, l, l, m, n, o, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o>;
  return Tensor6_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim2, Dim4, Dim5, Dim6, i, j, k, m, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, )*B(m, l, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<A, T, Dim0, Dim1, Dim2, Dim3, i, j, k, l> &a,
               const Tensor4_Expr<B, U, Dim4, Dim3, Dim5, Dim6, m, l, n, o> &b)
{
  using TensorExpr
    = Tensor4_times_Tensor4_single<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim3, Dim5, Dim6, i, j, k, l, m, l, n, o, Dim0, Dim1, Dim2, Dim4, Dim3, Dim5, Dim6, i, j, k, m, l, n, o>;
  return Tensor6_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim2, Dim4, Dim5, Dim6, i, j, k, m, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, )*B(m, n, l, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<A, T, Dim0, Dim1, Dim2, Dim3, i, j, k, l> &a,
               const Tensor4_Expr<B, U, Dim4, Dim5, Dim3, Dim6, m, n, l, o> &b)
{
  using TensorExpr
    = Tensor4_times_Tensor4_single<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim3, Dim6, i, j, k, l, m, n, l, o, Dim0, Dim1, Dim2, Dim4, Dim5, Dim3, Dim6, i, j, k, m, n, l, o>;
  return Tensor6_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim2, Dim4, Dim5, Dim6, i, j, k, m, n, o>(TensorExpr(a, b));
}

//A(i, j, k, l, )*B(m, n, o, l, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor4_Expr<A, T, Dim0, Dim1, Dim2, Dim3, i, j, k, l> &a,
               const Tensor4_Expr<B, U, Dim4, Dim5, Dim6, Dim3, m, n, o, l> &b)
{
  using TensorExpr
    = Tensor4_times_Tensor4_single<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim3, i, j, k, l, m, n, o, l, Dim0, Dim1, Dim2, Dim4, Dim5, Dim6, Dim3, i, j, k, m, n, o, l>;
  return Tensor6_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim2, Dim4, Dim5, Dim6, i, j, k, m, n, o>(TensorExpr(a, b));
}

} //namespace FTensor
