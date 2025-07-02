#pragma once

namespace FTensor
{
//Universal A(i, j, k, l, m, )*B(n, o, p, ) single contraction
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, int Dim7, char i, char j, char k, char l, char m, char n, char o, char p, int DimA, int DimB, int DimC, int DimD, int DimX, int DimY, int DimZ, char a, char b, char c, char d, char x, char y, char z>
class Tensor5_times_Tensor3_single
{
  Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> iterA;
  Tensor3_Expr<B, U, Dim5, Dim6, Dim7, n, o, p> iterB;

public:
  typename promote<T, U>::V operator()(const int N1, const int N2, const int N3, const int N4, const int N5, const int N6) const
  {
    typename promote<T, U>::V result(0);
    for(int xx = 0; xx < DimX; ++xx)
      {
        // Permutation is where the indices get checked.
        result += Permutation5<DimA, DimB, DimC, DimD, DimX, a, b, c, d, x>().eval(
          iterA, N1, N2, N3, N4, xx) * Permutation3<DimX, DimY, DimZ, x, y, z>().eval(iterB, xx, N5, N6);
      }
    return result;
  }

  Tensor5_times_Tensor3_single(
    const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &iter_a,
    const Tensor3_Expr<B, U, Dim5, Dim6, Dim7, n, o, p> &iter_b)
      : iterA(iter_a), iterB(iter_b)
  {}
};

//A(i, j, k, l, m, )*B(i, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor3_Expr<B, U, Dim0, Dim5, Dim6, i, n, o> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor3_single<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim0, Dim5, Dim6, i, j, k, l, m, i, n, o, Dim1, Dim2, Dim3, Dim4, Dim0, Dim5, Dim6, j, k, l, m, i, n, o>;
  return Tensor6_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, j, k, l, m, n, o>(TensorExpr(a, b));
}

//B(i, n, o, )*A(i, j, k, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor3_Expr<B, U, Dim0, Dim5, Dim6, i, n, o> &b,
               const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a)
{
  return a * b;
}

//A(i, j, k, l, m, )*B(n, i, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor3_Expr<B, U, Dim5, Dim0, Dim6, n, i, o> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor3_single<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim0, Dim6, i, j, k, l, m, n, i, o, Dim1, Dim2, Dim3, Dim4, Dim5, Dim0, Dim6, j, k, l, m, n, i, o>;
  return Tensor6_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, j, k, l, m, n, o>(TensorExpr(a, b));
}

//B(n, i, o, )*A(i, j, k, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor3_Expr<B, U, Dim5, Dim0, Dim6, n, i, o> &b,
               const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a)
{
  return a * b;
}

//A(i, j, k, l, m, )*B(n, o, i, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor3_Expr<B, U, Dim5, Dim6, Dim0, n, o, i> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor3_single<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim0, i, j, k, l, m, n, o, i, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim0, j, k, l, m, n, o, i>;
  return Tensor6_Expr<TensorExpr, typename promote<T,U>::V, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, j, k, l, m, n, o>(TensorExpr(a, b));
}

//B(n, o, i, )*A(i, j, k, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor3_Expr<B, U, Dim5, Dim6, Dim0, n, o, i> &b,
               const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a)
{
  return a * b;
}

//A(i, j, k, l, m, )*B(j, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor3_Expr<B, U, Dim1, Dim5, Dim6, j, n, o> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor3_single<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim1, Dim5, Dim6, i, j, k, l, m, j, n, o, Dim0, Dim2, Dim3, Dim4, Dim1, Dim5, Dim6, i, k, l, m, j, n, o>;
  return Tensor6_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim2, Dim3, Dim4, Dim5, Dim6, i, k, l, m, n, o>(TensorExpr(a, b));
}

//B(j, n, o, )*A(i, j, k, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor3_Expr<B, U, Dim1, Dim5, Dim6, j, n, o> &b,
               const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a)
{
  return a * b;
}

//A(i, j, k, l, m, )*B(n, j, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor3_Expr<B, U, Dim5, Dim1, Dim6, n, j, o> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor3_single<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim1, Dim6, i, j, k, l, m, n, j, o, Dim0, Dim2, Dim3, Dim4, Dim5, Dim1, Dim6, i, k, l, m, n, j, o>;
  return Tensor6_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim2, Dim3, Dim4, Dim5, Dim6, i, k, l, m, n, o>(TensorExpr(a, b));
}

//B(n, j, o, )*A(i, j, k, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor3_Expr<B, U, Dim5, Dim1, Dim6, n, j, o> &b,
               const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a)
{
  return a * b;
}

//A(i, j, k, l, m, )*B(n, o, j, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor3_Expr<B, U, Dim5, Dim6, Dim1, n, o, j> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor3_single<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim1, i, j, k, l, m, n, o, j, Dim0, Dim2, Dim3, Dim4, Dim5, Dim6, Dim1, i, k, l, m, n, o, j>;
  return Tensor6_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim2, Dim3, Dim4, Dim5, Dim6, i, k, l, m, n, o>(TensorExpr(a, b));
}

//B(n, o, j, )*A(i, j, k, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor3_Expr<B, U, Dim5, Dim6, Dim1, n, o, j> &b,
               const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a)
{
  return a * b;
}

//A(i, j, k, l, m, )*B(k, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor3_Expr<B, U, Dim2, Dim5, Dim6, k, n, o> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor3_single<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim2, Dim5, Dim6, i, j, k, l, m, k, n, o, Dim0, Dim1, Dim3, Dim4, Dim2, Dim5, Dim6, i, j, l, m, k, n, o>;
  return Tensor6_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim3, Dim4, Dim5, Dim6, i, j, l, m, n, o>(TensorExpr(a, b));
}

//B(k, n, o, )*A(i, j, k, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor3_Expr<B, U, Dim2, Dim5, Dim6, k, n, o> &b,
               const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a)
{
  return a * b;
}

//A(i, j, k, l, m, )*B(n, k, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor3_Expr<B, U, Dim5, Dim2, Dim6, n, k, o> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor3_single<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim2, Dim6, i, j, k, l, m, n, k, o, Dim0, Dim1, Dim3, Dim4, Dim5, Dim2, Dim6, i, j, l, m, n, k, o>;
  return Tensor6_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim3, Dim4, Dim5, Dim6, i, j, l, m, n, o>(TensorExpr(a, b));
}

//B(n, k, o, )*A(i, j, k, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor3_Expr<B, U, Dim5, Dim2, Dim6, n, k, o> &b,
               const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a)
{
  return a * b;
}

//A(i, j, k, l, m, )*B(n, o, k, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor3_Expr<B, U, Dim5, Dim6, Dim2, n, o, k> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor3_single<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim2, i, j, k, l, m, n, o, k, Dim0, Dim1, Dim3, Dim4, Dim5, Dim6, Dim2, i, j, l, m, n, o, k>;
  return Tensor6_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim3, Dim4, Dim5, Dim6, i, j, l, m, n, o>(TensorExpr(a, b));
}

//B(n, o, k, )*A(i, j, k, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor3_Expr<B, U, Dim5, Dim6, Dim2, n, o, k> &b,
               const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a)
{
  return a * b;
}

//A(i, j, k, l, m, )*B(l, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor3_Expr<B, U, Dim3, Dim5, Dim6, l, n, o> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor3_single<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim3, Dim5, Dim6, i, j, k, l, m, l, n, o, Dim0, Dim1, Dim2, Dim4, Dim3, Dim5, Dim6, i, j, k, m, l, n, o>;
  return Tensor6_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim2, Dim4, Dim5, Dim6, i, j, k, m, n, o>(TensorExpr(a, b));
}

//B(l, n, o, )*A(i, j, k, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor3_Expr<B, U, Dim3, Dim5, Dim6, l, n, o> &b,
               const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a)
{
  return a * b;
}

//A(i, j, k, l, m, )*B(n, l, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor3_Expr<B, U, Dim5, Dim3, Dim6, n, l, o> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor3_single<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim3, Dim6, i, j, k, l, m, n, l, o, Dim0, Dim1, Dim2, Dim4, Dim5, Dim3, Dim6, i, j, k, m, n, l, o>;
  return Tensor6_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim2, Dim4, Dim5, Dim6, i, j, k, m, n, o>(TensorExpr(a, b));
}

//B(n, l, o, )*A(i, j, k, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor3_Expr<B, U, Dim5, Dim3, Dim6, n, l, o> &b,
               const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a)
{
  return a * b;
}

//A(i, j, k, l, m, )*B(n, o, l, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor3_Expr<B, U, Dim5, Dim6, Dim3, n, o, l> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor3_single<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim3, i, j, k, l, m, n, o, l, Dim0, Dim1, Dim2, Dim4, Dim5, Dim6, Dim3, i, j, k, m, n, o, l>;
  return Tensor6_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim2, Dim4, Dim5, Dim6, i, j, k, m, n, o>(TensorExpr(a, b));
}

//B(n, o, l, )*A(i, j, k, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor3_Expr<B, U, Dim5, Dim6, Dim3, n, o, l> &b,
               const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a)
{
  return a * b;
}

//A(i, j, k, l, m, )*B(m, n, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor3_Expr<B, U, Dim4, Dim5, Dim6, m, n, o> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor3_single<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim4, Dim5, Dim6, i, j, k, l, m, m, n, o, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, i, j, k, l, m, n, o>;
  return Tensor6_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim2, Dim3, Dim5, Dim6, i, j, k, l, n, o>(TensorExpr(a, b));
}

//B(m, n, o, )*A(i, j, k, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor3_Expr<B, U, Dim4, Dim5, Dim6, m, n, o> &b,
               const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a)
{
  return a * b;
}

//A(i, j, k, l, m, )*B(n, m, o, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor3_Expr<B, U, Dim5, Dim4, Dim6, n, m, o> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor3_single<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim4, Dim6, i, j, k, l, m, n, m, o, Dim0, Dim1, Dim2, Dim3, Dim5, Dim4, Dim6, i, j, k, l, n, m, o>;
  return Tensor6_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim2, Dim3, Dim5, Dim6, i, j, k, l, n, o>(TensorExpr(a, b));
}

//B(n, m, o, )*A(i, j, k, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor3_Expr<B, U, Dim5, Dim4, Dim6, n, m, o> &b,
               const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a)
{
  return a * b;
}

//A(i, j, k, l, m, )*B(n, o, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a,
               const Tensor3_Expr<B, U, Dim5, Dim6, Dim4, n, o, m> &b)
{
  using TensorExpr
    = Tensor5_times_Tensor3_single<A, B, T, U, Dim0, Dim1, Dim2, Dim3, Dim4, Dim5, Dim6, Dim4, i, j, k, l, m, n, o, m, Dim0, Dim1, Dim2, Dim3, Dim5, Dim6, Dim4, i, j, k, l, n, o, m>;
  return Tensor6_Expr<TensorExpr, typename promote<T,U>::V, Dim0, Dim1, Dim2, Dim3, Dim5, Dim6, i, j, k, l, n, o>(TensorExpr(a, b));
}

//B(n, o, m, )*A(i, j, k, l, m, )
template<class A, class B, class T, class U, int Dim0, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5, int Dim6, char i, char j, char k, char l, char m, char n, char o>
auto operator*(const Tensor3_Expr<B, U, Dim5, Dim6, Dim4, n, o, m> &b,
               const Tensor5_Expr<A, T, Dim0, Dim1, Dim2, Dim3, Dim4, i, j, k, l, m> &a)
{
  return a * b;
}

} //namespace FTensor
